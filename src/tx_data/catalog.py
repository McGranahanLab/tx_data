"""Load per-table catalog YAMLs from ``catalog/``.

Files starting with ``_`` are treated as registries/metadata (e.g.
``_sources.yml``) and skipped by ``load_all_catalogs()``.
"""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import yaml

from tx_data.paths import REPO_ROOT

CATALOG_DIR = REPO_ROOT / "catalog"


def _catalog_path(name: str) -> Path:
    return CATALOG_DIR / f"{name}.yml"


@lru_cache(maxsize=32)
def load_catalog(name: str) -> dict:
    """Load one canonical-table catalog file by name (e.g. ``muttable``)."""
    p = _catalog_path(name)
    if not p.is_file():
        raise FileNotFoundError(f"No catalog entry at {p}")
    doc = yaml.safe_load(p.read_text())
    if not isinstance(doc, dict):
        raise ValueError(f"Catalog {name} is not a mapping")
    return doc


def load_all_catalogs() -> dict[str, dict]:
    """Load every table catalog (skipping files that start with `_`)."""
    out: dict[str, dict] = {}
    for p in sorted(CATALOG_DIR.glob("*.yml")):
        if p.name.startswith("_"):
            continue
        out[p.stem] = yaml.safe_load(p.read_text())
    return out


# Back-compat alias for pre-existing callers/tests.
load_all = load_all_catalogs


def catalog_column_names(name: str) -> list[str]:
    """Return the ordered list of source-schema column names for a table."""
    cat = load_catalog(name)
    cols = cat.get("columns", [])
    return [c["name"] for c in cols if not c.get("derived", False)]


def catalog_dtype_hints(name: str) -> dict[str, str]:
    """Return ``{column: sql_type}`` mapping from the catalog columns list."""
    cat = load_catalog(name)
    return {c["name"]: c.get("type", "VARCHAR") for c in cat.get("columns", [])}
