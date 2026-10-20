"""Build driver_list → canonical Parquet.

Boolean columns are read as TRUE/FALSE strings from the source; we coerce
to real bools. Coordinates are kept as-is (both hg19 and hg38 present).
"""

from __future__ import annotations

import polars as pl

from tx_data.builds._base import canonical_output_path, log
from tx_data.normalize import canonical_chr_expr
from tx_data.sources import resolve_source

TABLE = "driver_list"

_BOOL_COLS = [
    "driver_gene", "mut_driver", "lung_mut_driver",
    "CN_driver", "lung_CN_driver",
]


def build() -> pl.DataFrame:
    src = resolve_source(TABLE)
    log(f"read {src}")

    # Force chromosome columns to String — otherwise polars may infer i64
    # from the first N rows (all autosomes) and fail on X / Y / MT.
    df = pl.read_csv(
        src,
        separator="\t",
        null_values=["NA", ""],
        infer_schema_length=10000,
        schema_overrides={"chr_hg19": pl.String, "chr_hg38": pl.String},
    )

    # Coerce booleans (source has "TRUE"/"FALSE" strings).
    coercions = []
    for c in _BOOL_COLS:
        if c in df.columns:
            coercions.append(
                pl.when(pl.col(c).cast(pl.Utf8).str.to_uppercase() == "TRUE")
                .then(True)
                .when(pl.col(c).cast(pl.Utf8).str.to_uppercase() == "FALSE")
                .then(False)
                .otherwise(None)
                .alias(c)
            )
    if coercions:
        df = df.with_columns(coercions)

    # Normalise chromosome columns (strip `chr` prefix) so joins line up with
    # muttable.chr, alphamissense.chr, and alpaca.segment_chr.
    df = df.with_columns(
        [canonical_chr_expr("chr_hg19"), canonical_chr_expr("chr_hg38")]
    )

    out = canonical_output_path(TABLE)
    df.write_parquet(out)
    log(f"wrote {out} ({df.height} rows × {df.width} cols)")
    return df


if __name__ == "__main__":
    build()
