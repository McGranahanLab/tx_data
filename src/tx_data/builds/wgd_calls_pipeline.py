"""Build wgd_calls_pipeline (WES) → canonical Parquet.

Per-clone WGD calls for the WES cohort from darwins_analysis `WGD_pipeline` with `config.wes.toml` (step 05
`wgd_calls.tsv`), same columns and derived fields as `wgs_wgd_calls` (Piotr, 2026-11-26). The previous
release table stays `wgd_calls` until consumers are switched.

Where the pipeline output lives is a build parameter: `build(source_dir=...)`,
`scripts/build_wgd_calls_pipeline.py --source-dir <run>/wgd/05_table`, or TX_WES_WGD_CALLS_DIR (used by
build_all), e.g. `<darwins_analysis>/WGD_pipeline/output/wes/wgd/05_table`. Without one, the `_sources.yml`
entry (mock layout in nemo_mock/ only) is used. Built after the pipeline, which reads the tx_data WES tables.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl

from tx_data.builds.wgs_wgd_calls import _source, build_table

TABLE = "wgd_calls_pipeline"
ENV_DIR = "TX_WES_WGD_CALLS_DIR"


def build(source_dir: Path | None = None) -> pl.DataFrame:
    return build_table(TABLE, _source(source_dir, TABLE, ENV_DIR, "wes"))


if __name__ == "__main__":
    build()
