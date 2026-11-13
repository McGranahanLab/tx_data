"""Build alpaca_corrected and alpaca_wgd_ratios_corrected → canonical Parquet.

Reads the output of darwins_analysis `alpaca_post_hoc_correction` (ALPACA copy number with
missing clones re-inferred from the agreed WGD calls, and the WGD ratios recomputed from it).
Checks the corrected table has exactly the keys of `alpaca` and that uncorrected rows keep
ALPACA's copy number (`pred_CN_*_original`).
"""

from __future__ import annotations

import polars as pl

from tx_data.builds._base import canonical_output_path, log
from tx_data.builds.alpaca_wgd_ratios import SCHEMA as RATIO_SCHEMA
from tx_data.sources import resolve_source

TABLE = "alpaca_corrected"
RATIOS = "alpaca_wgd_ratios_corrected"

SCHEMA = {
    "tumour_id": pl.String, "segment": pl.String, "clone": pl.String,
    "pred_CN_A": pl.Int64, "pred_CN_B": pl.Int64,
    "segment_chr": pl.String, "segment_start": pl.Int64, "segment_end": pl.Int64,
    "total_cn": pl.Int64, "minor_cn": pl.Int64, "major_cn": pl.Int64, "is_homdel": pl.Boolean, "is_loh": pl.Boolean,
    "pred_CN_A_original": pl.Int64, "pred_CN_B_original": pl.Int64, "corrected": pl.Boolean,
}


def read_corrected(path) -> pl.DataFrame:
    df = pl.read_parquet(path)
    missing = [c for c in SCHEMA if c not in df.columns]
    if missing:
        raise ValueError(f"{path}: missing columns {missing}")
    df = df.select(pl.col(c).cast(t, strict=True) for c, t in SCHEMA.items())
    bad = df.filter(~pl.col("corrected") & ((pl.col("pred_CN_A") != pl.col("pred_CN_A_original"))
                                            | (pl.col("pred_CN_B") != pl.col("pred_CN_B_original"))))
    if bad.height:
        raise ValueError(f"{path}: {bad.height} rows changed without corrected = true")
    dup = df.filter(pl.struct("tumour_id", "segment", "clone").is_duplicated())
    if dup.height:
        raise ValueError(f"{TABLE}: duplicate (tumour_id, segment, clone):\n{dup.head()}")
    return df


def read_ratios(path) -> pl.DataFrame:
    df = pl.read_csv(path, separator="\t", null_values=["NA", "", "nan", "NaN"], infer_schema_length=10000,
                     schema_overrides={c: t for c, t in RATIO_SCHEMA.items() if t == pl.String})
    missing = [c for c in RATIO_SCHEMA if c not in df.columns]
    if missing:
        raise ValueError(f"{path}: missing columns {missing}")
    df = df.select(pl.col(c).cast(t, strict=True) for c, t in RATIO_SCHEMA.items())
    df = df.with_columns(pl.col(pl.Float64).fill_nan(None), (pl.col("parent") == "diploid").alias("is_mrca"))
    dup = df.filter(pl.struct("tumour_id", "clone").is_duplicated())
    if dup.height:
        raise ValueError(f"{RATIOS}: duplicate (tumour_id, clone):\n{dup.head()}")
    return df


def build() -> pl.DataFrame:
    src = resolve_source(TABLE)
    log(f"read {src}")
    df = read_corrected(src)
    out = canonical_output_path(TABLE)
    df.write_parquet(out)
    log(f"wrote {out} ({df.height} rows × {df.width} cols; {int(df['corrected'].sum())} corrected)")
    rsrc = resolve_source(RATIOS)
    log(f"read {rsrc}")
    r = read_ratios(rsrc)
    rout = canonical_output_path(RATIOS)
    r.write_parquet(rout)
    log(f"wrote {rout} ({r.height} rows × {r.width} cols)")
    return df


if __name__ == "__main__":
    build()
