"""Build wgs_alpaca → canonical Parquet.

ALPACA clone copy number for the WGS cohort, from
`wgs_resources/alpaca/output/tumour_results/<label>/ALPACA_output_<label>.parquet` (Piotr,
2026-11-17: same format as the WES output). Same columns and derived fields as `alpaca`
(`tx_data.builds.alpaca.derive`), plus the WGS id columns; tumour_id from the inventory.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl

from tx_data.builds._wgs_tumour_files import ID_COLS, build_per_tumour, require
from tx_data.builds.alpaca import derive

TABLE = "wgs_alpaca"
SCHEMA = {**ID_COLS, "segment": pl.String, "clone": pl.String, "pred_CN_A": pl.Int64, "pred_CN_B": pl.Int64,
          "segment_chr": pl.String, "segment_start": pl.Int64, "segment_end": pl.Int64, "total_cn": pl.Int64,
          "minor_cn": pl.Int64, "major_cn": pl.Int64, "is_homdel": pl.Boolean, "is_loh": pl.Boolean}


def read_one(label: str, path: Path) -> pl.DataFrame:
    df = pl.read_parquet(path)
    require(df, ("segment", "clone", "pred_CN_A", "pred_CN_B"), path)
    df = df.select(pl.col("segment").cast(pl.String), pl.col("clone").cast(pl.String),
                   pl.col("pred_CN_A").cast(pl.Int64), pl.col("pred_CN_B").cast(pl.Int64))
    return derive(df)


def build() -> pl.DataFrame:
    return build_per_tumour(TABLE, SCHEMA, read_one, ("tumour_id", "clone", "segment_chr", "segment_start"))


if __name__ == "__main__":
    build()
