"""WGS sample / tumour lookups backed by the built `wgs_battenberg_inventory`.

The inventory is the authoritative list of WGS tumour samples. Other WGS
sources name the same samples differently:

  inventory `TUMOUR_SAMPLE`   LTX0001_SU_T1-R1--38da873ff7
  tree_input `SAMPLE`         LTX0001_SU_T1.R1      (stem, '-' → '.')
  ALPACA cp_table columns     same as WES: full name (assumed, TBC)

`match_sample()` resolves any of those forms to the inventory
`sample_name_hash`, scoped to one tumour so short forms can't collide
across tumours. Name tokens are otherwise not interpreted.
"""

from __future__ import annotations

from functools import lru_cache

import polars as pl

from tx_data.builds._base import canonical_output_path

INVENTORY_TABLE = "wgs_battenberg_inventory"


def sample_keys(sample_name_hash: str) -> set[str]:
    """Every form a WGS sample name is known to appear in."""
    stem = sample_name_hash.split("--", 1)[0]
    return {sample_name_hash, stem, stem.replace("-", ".")}


@lru_cache(maxsize=1)
def inventory() -> pl.DataFrame:
    """The built inventory parquet (build `wgs_battenberg_inventory` first)."""
    p = canonical_output_path(INVENTORY_TABLE)
    if not p.is_file():
        raise FileNotFoundError(
            f"{p} not found. Build `{INVENTORY_TABLE}` before other WGS tables."
        )
    return pl.read_parquet(p)


@lru_cache(maxsize=1)
def _index() -> dict[tuple[str, str], str]:
    idx: dict[tuple[str, str], str] = {}
    for label, name in inventory().select("wgs_tumour_label", "sample_name_hash").iter_rows():
        for key in sample_keys(name):
            prev = idx.setdefault((label, key), name)
            if prev != name:
                raise ValueError(
                    f"Ambiguous WGS sample key {key!r} in tumour {label!r}: "
                    f"{prev!r} vs {name!r}"
                )
    return idx


def match_sample(wgs_tumour_label: str, source_name: str) -> str | None:
    """Inventory `sample_name_hash` for a source sample name, or None."""
    return _index().get((wgs_tumour_label, source_name))


def tumours() -> pl.DataFrame:
    """One row per WGS tumour: wgs_tumour_label, patient_id, tumour_id."""
    return (
        inventory()
        .select("wgs_tumour_label", "patient_id", "tumour_id")
        .unique()
        .sort("tumour_id")
    )


def clear_cache() -> None:
    """Drop cached inventory (after a rebuild in the same process)."""
    inventory.cache_clear()
    _index.cache_clear()
