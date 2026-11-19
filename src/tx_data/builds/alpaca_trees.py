"""Build alpaca_trees (WES) → canonical Parquet: one row per (tumour × clone), child -> parent.

The tree ALPACA used, from `ALPACA/input/<tumour_id>/tree_paths.json` (next to cp_table.csv). Not
CONIPHER's output: CONIPHER gives several alternative trees and the optimal one is selected before
ALPACA (Piotr, 2026-11-19). The MRCA's parent is `diploid`.
"""

from __future__ import annotations

import polars as pl

from tx_data.builds._alpaca_tree import SCHEMA, edges
from tx_data.builds._base import canonical_output_path, log
from tx_data.sources import resolve_source

TABLE = "alpaca_trees"


def build() -> pl.DataFrame:
    probe = resolve_source(TABLE, tumour_id="__probe__")
    root = probe.parent.parent
    rows, missing = [], []
    for d in sorted(p for p in root.iterdir() if p.is_dir()) if root.is_dir() else []:
        f = d / probe.name
        if not f.is_file():
            missing.append(d.name)
            continue
        rows += [{"tumour_id": d.name, "child": c, "parent": p} for c, p in edges(f)]
    log(f"{TABLE}: {len({r['tumour_id'] for r in rows})} tumours from {root}")
    if missing:
        log(f"warning: {TABLE}: no {probe.name} for {len(missing)} tumour folder(s): {missing}")
    df = pl.DataFrame(rows, schema=SCHEMA).sort("tumour_id", "child")
    out = canonical_output_path(TABLE)
    df.write_parquet(out)
    log(f"wrote {out} ({df.height} rows × {df.width} cols)")
    return df


if __name__ == "__main__":
    build()
