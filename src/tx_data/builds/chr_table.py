"""Build chr_table_hg19 and chr_table_hg38 → canonical Parquet: chromosome lengths, columns chr, len.

From UCSC's `<build>.chrom.sizes` (committed under references/<build>/, downloaded 2026-11-19 from
https://hgdownload.soe.ucsc.edu/goldenPath/<build>/bigZips/<build>.chrom.sizes). Primary chromosomes
only, 1-22, X, Y, in that order; `chr` without the `chr` prefix (tx_data's canonical form).
"""

from __future__ import annotations

import polars as pl

from tx_data.builds._base import canonical_output_path, log
from tx_data.sources import resolve_source

PRIMARY = [str(i) for i in range(1, 23)] + ["X", "Y"]
SCHEMA = {"chr": pl.String, "len": pl.Int64}


def build_one(build: str) -> pl.DataFrame:
    table = f"chr_table_{build}"
    src = resolve_source(table)
    raw = pl.read_csv(src, separator="\t", has_header=False, new_columns=["chrom", "len"],
                      schema_overrides={"chrom": pl.String, "len": pl.Int64})
    df = (raw.with_columns(pl.col("chrom").str.replace(r"^chr", "").alias("chr"))
             .filter(pl.col("chr").is_in(PRIMARY)).select("chr", "len"))
    if sorted(df["chr"]) != sorted(PRIMARY) or df["chr"].is_duplicated().any():
        raise ValueError(f"{src}: expected exactly the primary chromosomes {PRIMARY}")
    df = df.with_columns(pl.col("chr").replace_strict({c: i for i, c in enumerate(PRIMARY)}).alias("_o")).sort("_o").drop("_o")
    out = canonical_output_path(table)
    df.write_parquet(out)
    log(f"wrote {out} ({df.height} rows; genome {df['len'].sum():,} bp)")
    return df


def build() -> dict[str, pl.DataFrame]:
    return {b: build_one(b) for b in ("hg19", "hg38")}


if __name__ == "__main__":
    build()
