"""Build wgs_battenberg_inventory → canonical Parquet.

Source is a per-sample TSV listing each WGS tumour sample's chosen
Battenberg solution, purity, ploidy and pointers to the Battenberg results,
the multi-sample MSAI file and the mosdepth summary.

The WES `sample_hash_map` does not cover WGS samples, so `tumour_id` and
`region` are parsed best-effort from the sample name. Name tokens are
unreliable with no closed catalogue, so a name that doesn't follow
`<patient>_SU_T<n>-R<n>--<hash>` (or whose patient disagrees with
`PATIENT`) keeps its row with NULL `tumour_id` / `region` and is logged.
`patient_id` and `sample_name_hash` are always the source values.

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
from tx_data.ids import apply_convention_scalar
from tx_data.sources import resolve_source

TABLE = "wgs_battenberg_inventory"

SOURCE_RENAMES = {
    "PATIENT": "patient_id",
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
    "patient_id", "tumour_id", "region", "sample_name_hash",
    "results_dir", "optimal_solution", "subclones_file", "subclones_path",
    "purity", "ploidy", "msai_path", "mosdepth_summary_path", "coverage",
]


def _norm_path(p: str | None) -> str | None:
    """Collapse doubled slashes and drop a trailing slash; keep it absolute."""
    if p is None:
        return None
    return posixpath.normpath(p)


def _parse_sample(name: str) -> tuple[str, str, str]:
    """Return (patient_id, tumour_id, region) from a WGS sample name."""
    m = SAMPLE_RE.match(name)
    if m is None:
        raise ValueError(name)
    patient = m.group("patient")
    tumour_id = apply_convention_scalar(f"{patient}_{int(m.group('ordinal'))}")
    return patient, tumour_id, m.group("region")


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

    # ---- identity: best-effort parse; unparsed samples keep NULL tumour_id ----
    tumour_ids: dict[str, str | None] = {}
    regions: dict[str, str | None] = {}
    unparsed: list[str] = []
    mismatched: list[str] = []
    for patient, name in df.select("patient_id", "sample_name_hash").iter_rows():
        try:
            parsed_patient, tumour_id, region = _parse_sample(name)
        except ValueError:
            unparsed.append(name)
            tumour_ids[name] = regions[name] = None
            continue
        if parsed_patient != patient:
            mismatched.append(name)
            tumour_ids[name] = regions[name] = None
            continue
        tumour_ids[name], regions[name] = tumour_id, region
    if unparsed:
        log(
            f"{len(unparsed)} sample name(s) don't match <patient>_SU_T<n>-R<n>--<hash>; "
            f"tumour_id/region left NULL. First: {unparsed[:5]}"
        )
    if mismatched:
        log(
            f"{len(mismatched)} sample name(s) disagree with PATIENT; "
            f"tumour_id/region left NULL. First: {mismatched[:5]}"
        )

    df = df.with_columns(
        pl.col("sample_name_hash")
        .replace_strict(tumour_ids, return_dtype=pl.String)
        .alias("tumour_id"),
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

    out_df = df.select(OUTPUT_COLUMNS).sort(["patient_id", "tumour_id", "region", "sample_name_hash"], nulls_last=True)
    out = canonical_output_path(TABLE)
    out_df.write_parquet(out)
    log(f"wrote {out} ({out_df.height} rows × {out_df.width} cols)")
    return out_df


if __name__ == "__main__":
    build()
