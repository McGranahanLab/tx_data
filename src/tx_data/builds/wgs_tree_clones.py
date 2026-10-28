"""Build wgs_tree_clones → canonical Parquet.

Clones on each WGS tumour's ALPACA input phylogeny, read from
`wgs_resources/alpaca/input/<label>/tree_paths.json`: a JSON list of paths,
each a list of clone IDs from the trunk (MRCA) to a leaf
(`[["clone1", "clone2"], ["clone1", "clone3"]]`; single clone: `[["clone1"]]`).

One row per (tumour × clone on the tree). A tumour is single-clone when it
has exactly one clone here (Piotr, 2026-10-28: count clones on the tree, not
in cp_table); `n_clones_tree` carries that count on every row.

Tumours are discovered from the built `wgs_battenberg_inventory`.
"""

from __future__ import annotations

import json

import polars as pl

from tx_data import wgs_samples
from tx_data.builds._base import canonical_output_path, log
from tx_data.sources import resolve_source

TABLE = "wgs_tree_clones"

SCHEMA = {
    "tumour_id": pl.String,
    "patient_id": pl.String,
    "wgs_tumour_label": pl.String,
    "clone": pl.String,
    "is_trunk": pl.Boolean,
    "is_leaf": pl.Boolean,
    "depth": pl.Int64,
    "parent_clone": pl.String,
    "n_clones_tree": pl.Int64,
}


def _read_one(label: str, path) -> list[dict]:
    paths = json.loads(path.read_text())
    if (not isinstance(paths, list) or not paths
            or not all(isinstance(p, list) and p for p in paths)):
        raise ValueError(f"{path}: expected a non-empty list of non-empty clone paths")
    trunks = {str(p[0]) for p in paths}
    if len(trunks) != 1:
        raise ValueError(f"{path}: paths start at different clones {sorted(trunks)}")
    clones: dict[str, dict] = {}
    for p in paths:
        p = [str(c) for c in p]
        for depth, c in enumerate(p):
            parent = p[depth - 1] if depth else None
            prev = clones.get(c)
            if prev and (prev["depth"], prev["parent_clone"]) != (depth, parent):
                raise ValueError(f"{path}: clone {c} has inconsistent parents across paths")
            clones.setdefault(c, {"clone": c, "depth": depth, "parent_clone": parent,
                                  "is_trunk": depth == 0, "is_leaf": False})
        clones[p[-1]]["is_leaf"] = True
    for r in clones.values():
        r["wgs_tumour_label"] = label
        r["n_clones_tree"] = len(clones)
    return list(clones.values())


def build() -> pl.DataFrame:
    wgs_samples.clear_cache()
    tumours = wgs_samples.tumours()

    rows: list[dict] = []
    absent: list[str] = []
    n_read = 0
    for label in tumours.get_column("wgs_tumour_label").to_list():
        path = resolve_source(TABLE, wgs_tumour_label=label)
        if not path.is_file():
            absent.append(label)
            continue
        rows.extend(_read_one(label, path))
        n_read += 1
    log(f"read tree_paths.json for {n_read}/{tumours.height} WGS tumours")
    if absent:
        log(f"warning: no tree_paths.json for {len(absent)} tumour(s): {absent[:10]}")

    if not rows:
        df = pl.DataFrame(schema=SCHEMA)
    else:
        df = (
            pl.DataFrame(rows, schema={k: v for k, v in SCHEMA.items()
                                       if k not in ("tumour_id", "patient_id")})
            .join(tumours, on="wgs_tumour_label", how="left")
            .select(list(SCHEMA))
        )
        single = df.filter(pl.col("n_clones_tree") == 1)["tumour_id"].n_unique()
        log(f"{single}/{n_read} tumour(s) single-clone on the tree")

    df = df.sort("tumour_id", "depth", "clone")
    out = canonical_output_path(TABLE)
    df.write_parquet(out)
    log(f"wrote {out} ({df.height} rows × {df.width} cols)")
    return df


if __name__ == "__main__":
    build()
