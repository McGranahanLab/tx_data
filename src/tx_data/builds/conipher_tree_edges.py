"""Build conipher_tree_edges (WES) → canonical Parquet.

The default (first) tree of each WES tumour's CONIPHER `allTrees.txt` in the tx842 release tree
(`<patient>/mutation_trees/tumour_<n>/`; Piotr, 2026-11-24), one row per edge parent → child cluster.
WES counterpart of `wgs_conipher_tree_edges`, parsed the same way (`# tree k` blocks; a single-cluster
tree is CONIPHER's self-loop `1 1`, kept as such for WGD step 00). `n_trees` counts the alternative trees.
"""

from __future__ import annotations

import polars as pl

from tx_data.builds._wes_release_files import ID_COLS, build_per_tumour
from tx_data.builds.wgs_conipher_tree_edges import read_one as read_all_trees

TABLE = "conipher_tree_edges"
SCHEMA = {**ID_COLS, "parent_cluster": pl.Int64, "child_cluster": pl.Int64, "n_trees": pl.Int64}


def build() -> pl.DataFrame:
    return build_per_tumour(TABLE, SCHEMA, lambda p, n, path: read_all_trees(f"{p}_{n}", path),
                            ("tumour_id", "parent_cluster", "child_cluster"))


if __name__ == "__main__":
    build()
