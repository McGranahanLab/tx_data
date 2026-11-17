"""Build wgs_conipher_mutations → canonical Parquet.

CONIPHER mutations of each WGS tumour, from the 2026_08_mean_max_absent run (Piotr, 2026-11-17):
`<label>/<label>.SCoutput.CLEAN.tsv` (one row per mutation × sample) joined with
`<label>/treeTable.tsv` on (mutation_id, SAMPLE) for the tree cluster. Same shape as the WES
CONIPHER output. Input of the WGD pipeline's step 00 (PGDD input), via WGD_WGS_pipeline.

Kept: the columns step 00 reads (renamed to snake case) plus PLOIDY when present.
`tree_cluster` is null for mutations treeTable leaves out (clusters not on the tree).
"""

from __future__ import annotations

from pathlib import Path

import polars as pl

from tx_data.builds._wgs_tumour_files import ID_COLS, build_per_tumour, read_tsv, require, sample_hash
from tx_data.normalize import canonical_chr_expr
from tx_data.sources import resolve_source

TABLE = "wgs_conipher_mutations"
TREE_TABLE_PATTERN = "relative_path_pattern_tree_table"

SOURCE = {  # SCoutput column -> canonical (name, type)
    "SAMPLE": ("sample", pl.String), "CHR": ("chr", pl.String), "POS": ("pos", pl.Int64),
    "REF": ("ref", pl.String), "ALT": ("alt", pl.String), "mutation_id": ("mutation_id", pl.String),
    "CLUSTER": ("cluster", pl.Int64), "MUT_COPY": ("mut_copy", pl.Float64),
    "COPY_NUMBER_A": ("cn_a", pl.Float64), "COPY_NUMBER_B": ("cn_b", pl.Float64),
}
SCHEMA = {**ID_COLS, "sample": pl.String, "sample_name_hash": pl.String, "chr": pl.String, "pos": pl.Int64,
          "ref": pl.String, "alt": pl.String, "mutation_id": pl.String, "cluster": pl.Int64,
          "tree_cluster": pl.Int64, "mut_copy": pl.Float64, "cn_a": pl.Float64, "cn_b": pl.Float64,
          "ploidy": pl.Float64}


def read_one(label: str, path: Path) -> pl.DataFrame:
    muts = read_tsv(path)
    require(muts, SOURCE, path)
    ploidy = pl.col("PLOIDY").cast(pl.Float64) if "PLOIDY" in muts.columns else pl.lit(None, pl.Float64)
    df = muts.select(*[pl.col(s).cast(t).alias(n) for s, (n, t) in SOURCE.items()], ploidy.alias("ploidy"))
    tt_path = resolve_source(TABLE, pattern_key=TREE_TABLE_PATTERN, wgs_tumour_label=label)
    if not tt_path.is_file():
        raise FileNotFoundError(f"{tt_path}: treeTable.tsv missing next to {path.name}")
    tt = read_tsv(tt_path)
    require(tt, ("mutation_id", "SAMPLE", "treeCLUSTER"), tt_path)
    tt = tt.select(pl.col("mutation_id"), pl.col("SAMPLE").alias("sample"),
                   pl.col("treeCLUSTER").cast(pl.Int64).alias("tree_cluster")).unique(["mutation_id", "sample"])
    return (df.join(tt, on=["mutation_id", "sample"], how="left")
              .with_columns(canonical_chr_expr("chr"), sample_hash(label)))


def build() -> pl.DataFrame:
    return build_per_tumour(TABLE, SCHEMA, read_one, ("tumour_id", "sample", "chr", "pos"))


if __name__ == "__main__":
    build()
