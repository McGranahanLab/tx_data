"""Read ALPACA's `tree_paths.json` (the tree ALPACA was given) into child -> parent edges.

`tree_paths.json` is a JSON list of trunk -> leaf paths of clone IDs (`[["clone1","clone2"],
["clone1","clone3"]]`; single clone: `[["clone1"]]`). CONIPHER writes several alternative trees; the
lab selects the optimal one before ALPACA, so ALPACA's input holds the selected tree and CONIPHER's
own output is not used for trees (Piotr, 2026-11-19).
"""

from __future__ import annotations

import json
from pathlib import Path

import polars as pl

SCHEMA = {"tumour_id": pl.String, "child": pl.String, "parent": pl.String}


def edges(path: Path) -> list[tuple[str, str]]:
    """(child, parent) for every clone; the MRCA's parent is `diploid` (as in wgd_ratio_scores)."""
    paths = json.loads(Path(path).read_text())
    if not isinstance(paths, list) or not paths or not all(isinstance(p, list) and p for p in paths):
        raise ValueError(f"{path}: expected a non-empty list of non-empty clone paths")
    if len({str(p[0]) for p in paths}) != 1:
        raise ValueError(f"{path}: paths start at different clones")
    parent: dict[str, str] = {}
    for p in paths:
        for i, c in enumerate(map(str, p)):
            par = str(p[i - 1]) if i else "diploid"
            if parent.setdefault(c, par) != par:
                raise ValueError(f"{path}: clone {c} has two parents ({parent[c]}, {par})")
    return list(parent.items())
