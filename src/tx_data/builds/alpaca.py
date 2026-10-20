"""Build alpaca → canonical Parquet.

ALPACA cohort-wide CSV → long Parquet with parsed segment coordinates and
derived CN summaries.

Adds derived columns per `catalog/alpaca.yml`:
  segment_chr, segment_start, segment_end, total_cn, minor_cn, major_cn,
  is_homdel, is_loh.
"""

from __future__ import annotations

import polars as pl

from tx_data.builds._base import canonical_output_path, log
from tx_data.normalize import canonical_chr_expr
from tx_data.sources import resolve_source

TABLE = "alpaca"


def build() -> pl.DataFrame:
    src = resolve_source(TABLE)
    log(f"read {src}")

    df = pl.read_csv(
        src,
        separator=",",
        null_values=["NA", ""],
        infer_schema_length=10000,
        schema_overrides={
            "tumour_id": pl.String,
            "segment": pl.String,
            "clone": pl.String,
            "pred_CN_A": pl.Int64,
            "pred_CN_B": pl.Int64,
        },
    )

    # Parse `segment` = "<chr>_<start>_<end>".
    seg_parts = pl.col("segment").str.split_exact("_", 2).struct.rename_fields(
        ["segment_chr", "segment_start_s", "segment_end_s"]
    )
    df = df.with_columns(seg_parts.alias("_seg")).unnest("_seg")
    df = df.with_columns(
        [
            pl.col("segment_start_s").cast(pl.Int64).alias("segment_start"),
            pl.col("segment_end_s").cast(pl.Int64).alias("segment_end"),
            # Strip any `chr` prefix on segment_chr so it matches the
            # canonical no-prefix form used everywhere else.
            canonical_chr_expr("segment_chr"),
        ]
    ).drop(["segment_start_s", "segment_end_s"])

    # Derived CN summaries.
    df = df.with_columns(
        [
            (pl.col("pred_CN_A") + pl.col("pred_CN_B")).alias("total_cn"),
            pl.min_horizontal("pred_CN_A", "pred_CN_B").alias("minor_cn"),
            pl.max_horizontal("pred_CN_A", "pred_CN_B").alias("major_cn"),
        ]
    ).with_columns(
        [
            ((pl.col("pred_CN_A") == 0) & (pl.col("pred_CN_B") == 0)).alias("is_homdel"),
            ((pl.col("minor_cn") == 0) & (pl.col("total_cn") > 0)).alias("is_loh"),
        ]
    )

    out = canonical_output_path(TABLE)
    df.write_parquet(out)
    log(f"wrote {out} ({df.height} rows × {df.width} cols)")
    return df


if __name__ == "__main__":
    build()
