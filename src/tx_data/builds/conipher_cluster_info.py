"""Build conipher_cluster_info (WES) → canonical Parquet.

CONIPHER `clusterInfo.txt` of each WES tumour in the tx842 release tree (`<patient>/mutation_trees/
tumour_<n>/`; Piotr, 2026-11-24), reduced to the columns the WGD pipeline's CCF QC reads: the lowest
`meanCCF` over truncal tree clusters must be above 90%. One row per cluster (× sample when the file has
SAMPLE; R-style names converted to the dash form). WES counterpart of `wgs_conipher_cluster_info`.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl

from tx_data.builds._wes_release_files import ID_COLS, build_per_tumour, dash_sample
from tx_data.builds._wgs_tumour_files import r_bool, read_tsv, require

TABLE = "conipher_cluster_info"
SCHEMA = {**ID_COLS, "cluster_id": pl.Int64, "sample": pl.String, "truncal": pl.Boolean, "tree_clust": pl.Boolean,
          "mean_ccf": pl.Float64}


def read_one(patient: str, n: str, path: Path) -> pl.DataFrame:
    df = read_tsv(path)
    require(df, ("clusterID", "truncal", "treeClust", "meanCCF"), path)
    sample = pl.col("SAMPLE") if "SAMPLE" in df.columns else pl.lit(None, pl.String)
    return (df.with_columns(r_bool("truncal"), r_bool("treeClust"))
              .select(pl.col("clusterID").cast(pl.Int64).alias("cluster_id"), sample.alias("sample"),
                      pl.col("truncal"), pl.col("treeClust").alias("tree_clust"),
                      pl.col("meanCCF").cast(pl.Float64, strict=False).alias("mean_ccf"))
              .with_columns(dash_sample()))


def build() -> pl.DataFrame:
    return build_per_tumour(TABLE, SCHEMA, read_one, ("tumour_id", "cluster_id", "sample"))


if __name__ == "__main__":
    build()
