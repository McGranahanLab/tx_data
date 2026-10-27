"""Build wgs_battenberg_inventory → canonical Parquet.

Source is a per-sample TSV listing each WGS tumour sample's chosen
Battenberg solution, purity, ploidy and pointers to the Battenberg results,
the multi-sample MSAI file and the mosdepth summary.

The WES `sample_hash_map` does not cover WGS samples, so `tumour_id` and
`region` are parsed from the sample name itself. Any name that does not
follow `<patient>_SU_T<n>-R<n>--<hash>` fails the build loudly rather than
being dropped: an unparsed sample would silently vanish from every
downstream WGS analysis.

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
SOLUTION_RE = re.compile(r"_psi(?P<psi>\d+(?:\.\d+)?)_rho(?P<rho>\d+(?:\.\d+)?)$")

OUTPUT_COLUMNS = [
    "patient_id", "tumour_id", "region", "sample_name_hash",
    "results_dir", "optimal_solution", "solution_psi", "solution_rho",
    "subclones_file", "subclones_path",
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


def _parse_solution(label: str | None) -> tuple[float | None, float | None]:
    if label is None:
        return None, None
    m = SOLUTION_RE.search(label)
    if m is None:
        return None, None
    return float(m.group("psi")), float(m.group("rho"))


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

    # ---- identity: parse every sample name, fail on anything unexpected ----
    parsed: dict[str, tuple[str, str, str]] = {}
    bad: list[str] = []
    for name in df.get_column("sample_name_hash").to_list():
        try:
            parsed[name] = _parse_sample(name)
        except ValueError:
            bad.append(name)
    if bad:
        raise ValueError(
            f"{TABLE}: {len(bad)} sample name(s) don't match "
            f"<patient>_SU_T<n>-R<n>--<hash>; first: {bad[:5]}"
        )

    df = df.with_columns(
        pl.col("sample_name_hash")
        .replace_strict({k: v[0] for k, v in parsed.items()}, return_dtype=pl.String)
        .alias("_parsed_patient"),
        pl.col("sample_name_hash")
        .replace_strict({k: v[1] for k, v in parsed.items()}, return_dtype=pl.String)
        .alias("tumour_id"),
        pl.col("sample_name_hash")
        .replace_strict({k: v[2] for k, v in parsed.items()}, return_dtype=pl.String)
        .alias("region"),
    )
    mismatch = df.filter(pl.col("patient_id") != pl.col("_parsed_patient"))
    if mismatch.height:
        raise ValueError(
            f"{TABLE}: PATIENT disagrees with the sample name for {mismatch.height} row(s); "
            f"first: {mismatch.select('patient_id', 'sample_name_hash').head(5).to_dicts()}"
        )

    dupes = df.group_by("sample_name_hash").len().filter(pl.col("len") > 1)
    if dupes.height:
        raise ValueError(
            f"{TABLE}: {dupes.height} duplicate sample rows; "
            f"first: {dupes.head(5).to_dicts()}"
        )

    # ---- paths + solution label ---------------------------------------------
    path_cols = ["results_dir", "msai_path", "mosdepth_summary_path"]
    df = df.with_columns(
        [pl.col(c).map_elements(_norm_path, return_dtype=pl.String) for c in path_cols]
    )
    df = df.with_columns(
        pl.when(pl.col("results_dir").is_not_null() & pl.col("subclones_file").is_not_null())
        .then(pl.col("results_dir") + "/" + pl.col("subclones_file"))
        .alias("subclones_path"),
        pl.col("optimal_solution")
        .map_elements(lambda s: _parse_solution(s)[0], return_dtype=pl.Float64)
        .alias("solution_psi"),
        pl.col("optimal_solution")
        .map_elements(lambda s: _parse_solution(s)[1], return_dtype=pl.Float64)
        .alias("solution_rho"),
        pl.col("purity").cast(pl.Float64),
        pl.col("ploidy").cast(pl.Float64),
        pl.col("coverage").cast(pl.Float64),
    )

    unlabelled = df.filter(
        pl.col("optimal_solution").is_not_null() & pl.col("solution_psi").is_null()
    )
    if unlabelled.height:
        log(
            f"{unlabelled.height} optimal_solution label(s) without _psi/_rho tokens; "
            f"first: {unlabelled.get_column('optimal_solution').head(3).to_list()}"
        )

    out_of_range = df.filter(~pl.col("purity").is_between(0.0, 1.0, closed="right"))
    if out_of_range.height:
        log(f"warning: {out_of_range.height} sample(s) with purity outside (0, 1]")

    out_df = df.select(OUTPUT_COLUMNS).sort(["tumour_id", "region", "sample_name_hash"])
    out = canonical_output_path(TABLE)
    out_df.write_parquet(out)
    log(f"wrote {out} ({out_df.height} rows × {out_df.width} cols)")
    return out_df


if __name__ == "__main__":
    build()
