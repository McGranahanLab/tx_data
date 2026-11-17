"""Build wgs_alpaca_wgd_ratios → canonical Parquet.

ALPACA per-clone WGD ratios (`alpaca wgd`) for the WGS cohort, from
`wgs_resources/alpaca/output/tumour_results/<label>/wgd_ratio_scores.parquet` (Piotr,
2026-11-17: same format as WES). Same 12 columns + `is_mrca` as `alpaca_wgd_ratios`, plus the
WGS id columns; tumour_id from the inventory.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl

from tx_data.builds._wgs_tumour_files import ID_COLS, build_per_tumour, require
from tx_data.builds.alpaca_wgd_ratios import SCHEMA as WES_SCHEMA

TABLE = "wgs_alpaca_wgd_ratios"
SCHEMA = {**ID_COLS, **{k: v for k, v in WES_SCHEMA.items() if k != "tumour_id"}, "is_mrca": pl.Boolean}


def read_one(label: str, path: Path) -> pl.DataFrame:
    df = pl.read_parquet(path)
    cols = [c for c in WES_SCHEMA if c != "tumour_id"]
    require(df, cols, path)
    df = df.select(pl.col(c).cast(WES_SCHEMA[c]) for c in cols).with_columns(pl.col(pl.Float64).fill_nan(None))
    if df.select(pl.struct("clone").is_duplicated().any()).item():
        raise ValueError(f"{path}: duplicate clones")
    return df.with_columns((pl.col("parent") == "diploid").alias("is_mrca"))


def build() -> pl.DataFrame:
    return build_per_tumour(TABLE, SCHEMA, read_one, ("tumour_id", "clone"))


if __name__ == "__main__":
    build()
