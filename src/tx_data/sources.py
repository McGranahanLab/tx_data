"""Load the ``catalog/_sources.yml`` registry of HPC source locations."""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import yaml

from tx_data.paths import (
    REPO_ROOT,
    nemo_path,
    resolve_source_path,
    resolve_source_pattern,
)

SOURCES_YML = REPO_ROOT / "catalog" / "_sources.yml"

# `cohort: reference` marks cohort-independent reference data (hg19 tables).
REFERENCE_COHORT = "reference"


@lru_cache(maxsize=1)
def load_sources() -> dict:
    """Load the sources registry as a dict."""
    if not SOURCES_YML.is_file():
        raise FileNotFoundError(f"Missing sources registry: {SOURCES_YML}")
    doc = yaml.safe_load(SOURCES_YML.read_text())
    if not isinstance(doc, dict):
        raise ValueError("_sources.yml is not a mapping")
    return doc


def resolve_source(name: str, pattern_key: str = "relative_path_pattern", **fmt: str) -> Path:
    """Resolve a table's source path (or pattern), applying `{placeholders}`.

    `pattern_key` selects an alternative pattern declared on the same entry
    (e.g. `relative_path_pattern_single_sample`).

    Uses `alt_source_roots` from the registry when the table declares
    `source_root:`.
    """
    doc = load_sources()
    tables = doc.get("sources", {})
    if name not in tables:
        raise KeyError(f"Table {name!r} not in _sources.yml `sources:`")

    entry = tables[name]
    alt_roots = doc.get("alt_source_roots", {})
    cohort = source_cohort(name)
    cohort_root = (
        None if cohort == REFERENCE_COHORT
        else doc["cohorts"][cohort]["source_root_default"]
    )

    if "relative_path" in entry:
        return resolve_source_path(
            entry["relative_path"],
            source_root=entry.get("source_root"),
            cohort_root=cohort_root,
            alt_source_roots=alt_roots,
        )
    if pattern_key in entry:
        return resolve_source_pattern(
            entry[pattern_key],
            source_root=entry.get("source_root"),
            cohort_root=cohort_root,
            alt_source_roots=alt_roots,
            **fmt,
        )
    raise ValueError(
        f"Table {name!r} has neither `relative_path` nor `{pattern_key}`"
    )


def source_entry(name: str) -> dict:
    """Return the raw registry entry for a table (for format/comment-prefix etc.)."""
    doc = load_sources()
    tables = doc.get("sources", {})
    if name not in tables:
        raise KeyError(f"Table {name!r} not in _sources.yml `sources:`")
    return tables[name]


def source_cohort(name: str) -> str:
    """Return the cohort key (`wes`, `wgs`, ...) or `reference` for a table."""
    doc = load_sources()
    entry = source_entry(name)
    cohort = entry.get("cohort")
    known = set(doc.get("cohorts", {})) | {REFERENCE_COHORT}
    if cohort not in known:
        raise ValueError(
            f"Table {name!r} has cohort {cohort!r}; expected one of {sorted(known)}"
        )
    return cohort


def cohort_root(cohort: str) -> Path:
    """Resolve a cohort's main source tree under the current `nemo_root()`."""
    cohorts = load_sources().get("cohorts", {})
    if cohort not in cohorts:
        raise KeyError(f"Cohort {cohort!r} not in _sources.yml `cohorts:`")
    return nemo_path(cohorts[cohort]["source_root_default"])
