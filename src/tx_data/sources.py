"""Load the ``catalog/_sources.yml`` registry of HPC source locations."""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import yaml

from tx_data.paths import (
    REPO_ROOT,
    resolve_source_path,
    resolve_source_pattern,
)

SOURCES_YML = REPO_ROOT / "catalog" / "_sources.yml"


@lru_cache(maxsize=1)
def load_sources() -> dict:
    """Load the sources registry as a dict."""
    if not SOURCES_YML.is_file():
        raise FileNotFoundError(f"Missing sources registry: {SOURCES_YML}")
    doc = yaml.safe_load(SOURCES_YML.read_text())
    if not isinstance(doc, dict):
        raise ValueError("_sources.yml is not a mapping")
    return doc


def resolve_source(name: str, **fmt: str) -> Path:
    """Resolve a table's source path (or pattern), applying `{placeholders}`.

    Uses `alt_source_roots` from the registry when the table declares
    `source_root:`.
    """
    doc = load_sources()
    tables = doc.get("sources", {})
    if name not in tables:
        raise KeyError(f"Table {name!r} not in _sources.yml `sources:`")

    entry = tables[name]
    alt_roots = doc.get("alt_source_roots", {})

    if "relative_path" in entry:
        return resolve_source_path(
            entry["relative_path"],
            source_root=entry.get("source_root"),
            alt_source_roots=alt_roots,
        )
    if "relative_path_pattern" in entry:
        return resolve_source_pattern(
            entry["relative_path_pattern"],
            source_root=entry.get("source_root"),
            alt_source_roots=alt_roots,
            **fmt,
        )
    raise ValueError(
        f"Table {name!r} has neither `relative_path` nor `relative_path_pattern`"
    )


def source_entry(name: str) -> dict:
    """Return the raw registry entry for a table (for format/comment-prefix etc.)."""
    doc = load_sources()
    tables = doc.get("sources", {})
    if name not in tables:
        raise KeyError(f"Table {name!r} not in _sources.yml `sources:`")
    return tables[name]
