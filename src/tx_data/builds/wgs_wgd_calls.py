"""Build wgs_wgd_calls → canonical Parquet.

Per-clone WGD calls for the WGS cohort, from darwins_analysis `WGD_WGS_pipeline` (step 05
`wgd_calls.tsv`: every call, with `quality_pass` for the quality cut-off; Piotr, 2026-11-19). Same
derived columns as the WES `wgd_calls` (`tx_data.builds.wgd_calls.add_derived`); they are null on rows
without calls (status `qc_failed` / `no_alpaca_tree`).

Where the pipeline output lives is a build parameter, as for `alpaca_corrected` (nothing hardcoded):
`build(source_dir=...)`, `scripts/build_wgs_wgd_calls.py --source-dir <run>/wgd/05_table`, or the
environment variable TX_WGS_WGD_CALLS_DIR (used by build_all). Without one, the `_sources.yml` entry
(mock layout in nemo_mock/) is used.
"""

from __future__ import annotations

import os
from pathlib import Path

import polars as pl

from tx_data.builds._base import canonical_output_path, log
from tx_data.builds.wgd_calls import add_derived
from tx_data.sources import resolve_source

TABLE = "wgs_wgd_calls"
ENV_DIR = "TX_WGS_WGD_CALLS_DIR"
NO_CALLS = ("qc_failed", "no_alpaca_tree")

SOURCE_SCHEMA = {
    "tumour_id": pl.String, "clone": pl.String, "parent": pl.String, "number_of_gds_at_clone": pl.Int64,
    "total_gds_per_tumour": pl.Int64, "quality": pl.Float64, "status": pl.String, "class": pl.String,
    "sub_class": pl.String, "quality_code": pl.String, "worst_event_quality": pl.Float64,
    "worst_event_quality_code": pl.String, "n_unplaced_events": pl.Int64, "n_rejected_events": pl.Int64,
    "ploidy_prior": pl.String, "quality_pass": pl.Boolean, "upstream_qc": pl.String, "ccf_qc": pl.String,
    "truncal_min_ccf": pl.Float64, "qc_pass": pl.Boolean, "qc_reason": pl.String,
}
DERIVED = {"is_root": pl.Boolean, "is_wgd_event_clone": pl.Boolean, "cumulative_gds_to_clone": pl.Int64,
           "has_wgd_ancestor": pl.Boolean}
SCHEMA = {**SOURCE_SCHEMA, **DERIVED}


def _source(source_dir: Path | None) -> Path:
    d = source_dir or (Path(os.environ[ENV_DIR]) if os.environ.get(ENV_DIR) else None)
    return Path(d) / "wgd_calls.tsv" if d is not None else resolve_source(TABLE)


def read_calls(path: Path) -> pl.DataFrame:
    df = pl.read_csv(path, separator="\t", null_values=["NA", ""], infer_schema_length=None,
                     schema_overrides={c: t for c, t in SOURCE_SCHEMA.items() if t == pl.String})
    missing = [c for c in SOURCE_SCHEMA if c not in df.columns]
    if missing:
        raise ValueError(f"{path}: missing columns {missing} (WGD_WGS_pipeline step 05 from 2026-11-19 on)")
    df = df.select(pl.col(c).cast(t, strict=True) for c, t in SOURCE_SCHEMA.items())
    dup = df.filter(pl.col("clone").is_not_null() & pl.struct("tumour_id", "clone").is_duplicated())
    if dup.height:
        raise ValueError(f"{TABLE}: duplicate (tumour_id, clone):\n{dup.head()}")
    return df


def build(source_dir: Path | None = None) -> pl.DataFrame:
    src = _source(source_dir)
    log(f"read {src}")
    df = read_calls(src)
    calls = df.filter(~pl.col("status").is_in(list(NO_CALLS)) & pl.col("clone").is_not_null())
    rest = df.join(calls.select("tumour_id", "clone"), on=["tumour_id", "clone"], how="anti", nulls_equal=True)
    df = (pl.concat([add_derived(calls) if calls.height else calls.with_columns(
                         *[pl.lit(None, t).alias(c) for c, t in DERIVED.items()]),
                     rest.with_columns(*[pl.lit(None, t).alias(c) for c, t in DERIVED.items()])],
                    how="diagonal_relaxed")
            .select([pl.col(c).cast(t) for c, t in SCHEMA.items()])
            .sort("tumour_id", "clone", nulls_last=True))
    log(f"{df['tumour_id'].n_unique()} tumours; status: "
        + ", ".join(f"{s}={n}" for s, n in df.unique("tumour_id").group_by("status").len().sort("status").iter_rows()))
    out = canonical_output_path(TABLE)
    df.write_parquet(out)
    log(f"wrote {out} ({df.height} rows × {df.width} cols)")
    return df


if __name__ == "__main__":
    build()
