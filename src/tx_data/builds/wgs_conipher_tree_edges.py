"""Build wgs_conipher_tree_edges → canonical Parquet.

The default (first) tree of each WGS tumour's CONIPHER `allTrees.txt` (2026_08_mean_max_absent
run), one row per edge parent → child cluster. A single-cluster tree is written by CONIPHER as
the self-loop `1 1` and kept as such (WGD step 00 reads it as single-clone). `n_trees` counts the
alternative trees in the file.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl

from tx_data.builds._wgs_tumour_files import ID_COLS, build_per_tumour

TABLE = "wgs_conipher_tree_edges"
SCHEMA = {**ID_COLS, "parent_cluster": pl.Int64, "child_cluster": pl.Int64, "n_trees": pl.Int64}


def read_one(label: str, path: Path) -> pl.DataFrame:
    lines = path.read_text().splitlines()
    starts = [i for i, ln in enumerate(lines) if ln.startswith("# tree ")]
    if starts:
        end = starts[1] if len(starts) > 1 else len(lines)
        lines = lines[starts[0] + 1:end]
    edges = [ln.split("\t")[:2] for ln in lines if ln.strip() and not ln.startswith("#")]
    if not edges or any(len(e) < 2 for e in edges):
        raise ValueError(f"{path}: no readable edges in the first tree")
    return pl.DataFrame({"parent_cluster": [int(p) for p, _ in edges], "child_cluster": [int(c) for _, c in edges],
                         "n_trees": [max(len(starts), 1)] * len(edges)})


def build() -> pl.DataFrame:
    return build_per_tumour(TABLE, SCHEMA, read_one, ("tumour_id", "parent_cluster", "child_cluster"))


if __name__ == "__main__":
    build()
