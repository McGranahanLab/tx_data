"""Build wgs_alpaca_trees (WGS) → canonical Parquet: one row per (tumour × clone), child -> parent.

The tree ALPACA used, from `wgs_resources/alpaca/input/<label>/tree_paths.json` (the same files as
`wgs_tree_clones`, here in the three-column form of `alpaca_trees`). Not CONIPHER's output (Piotr,
2026-11-19). The MRCA's parent is `diploid`; tumour_id from the inventory.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl

from tx_data.builds._alpaca_tree import SCHEMA, edges
from tx_data.builds._wgs_tumour_files import build_per_tumour

TABLE = "wgs_alpaca_trees"


def read_one(label: str, path: Path) -> pl.DataFrame:
    e = edges(path)
    return pl.DataFrame({"child": [c for c, _ in e], "parent": [p for _, p in e]})


def build() -> pl.DataFrame:
    return build_per_tumour(TABLE, SCHEMA, read_one, ("tumour_id", "child"))

if __name__ == "__main__":
    build()
