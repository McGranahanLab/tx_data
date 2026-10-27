"""Build wgs_battenberg_inventory → canonical Parquet.

Source is a per-sample TSV listing each WGS tumour sample's chosen
Battenberg solution, purity, ploidy and pointers to the Battenberg results,
the multi-sample MSAI file and the mosdepth summary.

Identity comes from the source `PATIENT` column, which despite its name
holds a WGS *tumour* label: bare `LTX0001` for a patient's only tumour,
`LTX0002_tumour1` / `LTX0002_tumour2` for multi-tumour patients (also the
Battenberg output folder names). The build keeps it verbatim as
`wgs_tumour_label` and derives WES-consistent `patient_id` and `tumour_id`
(`LTX0001-Tumour1`) via `tx_data.ids.wgs_tumour_id`. An unrecognised label
fails the build.

`region` is parsed best-effort from the sample name: name tokens are
unreliable with no closed catalogue, so it is NULL (and logged) when the
name doesn't follow `<patient>_SU_T<n>-R<n>--<hash>`. Names whose tokens
disagree with the label are logged but keep the label-derived identity.

`optimal_solution` is kept verbatim; its psi / rho tokens are not parsed
into columns because `purity` / `ploidy` are the values to use.

Path columns stay as absolute `/nemo/...` strings (normalised); callers
resolve them with `tx_data.paths.nemo_path()`.
"""

from __future__ import annotations

import posixpath
import re

import polars as pl

from tx_data.builds._base import canonical_output_path, log
from tx_data.ids import parse_wgs_tumour_label, wgs_tumour_id
from tx_data.sources import resolve_source

TABLE = "wgs_battenberg_inventory"

SOURCE_RENAMES = {
    "PATIENT": "wgs_tumour_label",
    "TUMOUR_SAMPLE": "sample_name_hash",
    "results_dir": "results_dir",
    "optimal_solution": "optimal_solution",
    "subclones": "subclones_file",
    "purity": "purity",
    "ploidy": "ploidy",
    "msai": "msai_path",
    "mosdepth_summary": "mosdepth_summary_path",
    "coverage": "coverage",
}

SAMPLE_RE = re.compile(
    r"^(?P<patient>(?:LTX|CRUK)\d+)_SU_T(?P<ordinal>\d+)-(?P<region>R\d+)--(?P<hash>.+)$"
)

OUTPUT_COLUMNS = [
    "patient_id", "tumour_id", "wgs_tumour_label", "region", "sample_name_hash",
    "results_dir", "optimal_solution", "subclones_file", "subclones_path",
    "purity", "ploidy", "msai_path", "mosdepth_summary_path", "coverage",
]


def _norm_path(p: str | None) -> str | None:
    """Collapse doubled slashes and drop a trailing slash; keep it absolute."""
    if p is None:
        return None
    return posixpath.normpath(p)


def _parse_sample(name: str) -> tuple[str, int, str]:
    """Return (patient_id, tumour_ordinal, region) from a WGS sample name."""
    m = SAMPLE_RE.match(name)
    if m is None:
        raise ValueError(name)
    return m.group("patient"), int(m.group("ordinal")), m.group("region")


def build() -> pl.DataFrame:
    src = resolve_source(TABLE)
    log(f"read {src}")

    df = pl.read_csv(
        src,
        separator="\t",
        null_values=["NA", "", "nan"],
        infer_schema_length=10000,
        schema_overrides={"PATIENT": pl.String, "TUMOUR_SAMPLE": pl.String},
    )

    missing = set(SOURCE_RENAMES) - set(df.columns)
    if missing:
        raise ValueError(
            f"{TABLE} source missing expected columns {sorted(missing)}; got {df.columns!r}"
        )
    df = df.select(list(SOURCE_RENAMES)).rename(SOURCE_RENAMES)

    # ---- identity: tumour label → patient_id / tumour_id (authoritative) ----
    labels = df.get_column("wgs_tumour_label").unique().to_list()
    bad_labels = []
    for label in labels:
        try:
            parse_wgs_tumour_label(label)
        except (ValueError, TypeError):
            bad_labels.append(label)
    if bad_labels:
        raise ValueError(
            f"{TABLE}: {len(bad_labels)} PATIENT (tumour label) value(s) not of the form "
            f"<LTX|CRUK><digits>[_tumour<n>]; first: {bad_labels[:5]}"
        )
    df = df.with_columns(
        pl.col("wgs_tumour_label")
        .replace_strict({lb: parse_wgs_tumour_label(lb)[0] for lb in labels},
                        return_dtype=pl.String)
        .alias("patient_id"),
        pl.col("wgs_tumour_label")
        .replace_strict({lb: wgs_tumour_id(lb) for lb in labels}, return_dtype=pl.String)
        .alias("tumour_id"),
    )

    # ---- region: best-effort from the sample name ----------------------------
    regions: dict[str, str | None] = {}
    unparsed: list[str] = []
    disagree: list[str] = []
    for label, name in df.select("wgs_tumour_label", "sample_name_hash").iter_rows():
        try:
            name_patient, name_ordinal, region = _parse_sample(name)
        except ValueError:
            unparsed.append(name)
            regions[name] = None
            continue
        regions[name] = region
        if (name_patient, name_ordinal) != parse_wgs_tumour_label(label):
            disagree.append(f"{label} / {name}")
    if unparsed:
        log(
            f"{len(unparsed)} sample name(s) don't match <patient>_SU_T<n>-R<n>--<hash>; "
            f"region left NULL. First: {unparsed[:5]}"
        )
    if disagree:
        log(
            f"warning: {len(disagree)} sample name(s) whose tokens disagree with the "
            f"tumour label (label wins). First: {disagree[:5]}"
        )
    df = df.with_columns(
        pl.col("sample_name_hash")
        .replace_strict(regions, return_dtype=pl.String)
        .alias("region"),
    )

    dupes = df.group_by("sample_name_hash").len().filter(pl.col("len") > 1)
    if dupes.height:
        raise ValueError(
            f"{TABLE}: {dupes.height} duplicate sample rows; "
            f"first: {dupes.head(5).to_dicts()}"
        )

    # ---- paths ---------------------------------------------------------------
    path_cols = ["results_dir", "msai_path", "mosdepth_summary_path"]
    df = df.with_columns(
        [pl.col(c).map_elements(_norm_path, return_dtype=pl.String) for c in path_cols]
    )
    df = df.with_columns(
        pl.when(pl.col("results_dir").is_not_null() & pl.col("subclones_file").is_not_null())
        .then(pl.col("results_dir") + "/" + pl.col("subclones_file"))
        .alias("subclones_path"),
        pl.col("purity").cast(pl.Float64),
        pl.col("ploidy").cast(pl.Float64),
        pl.col("coverage").cast(pl.Float64),
    )

    out_of_range = df.filter(~pl.col("purity").is_between(0.0, 1.0, closed="right"))
    if out_of_range.height:
        log(f"warning: {out_of_range.height} sample(s) with purity outside (0, 1]")

    out_df = df.select(OUTPUT_COLUMNS).sort(["tumour_id", "region", "sample_name_hash"], nulls_last=True)
    out = canonical_output_path(TABLE)
    out_df.write_parquet(out)
    log(f"wrote {out} ({out_df.height} rows × {out_df.width} cols)")
    return out_df


if __name__ == "__main__":
    build()
