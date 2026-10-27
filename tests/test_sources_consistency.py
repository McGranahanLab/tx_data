"""Cross-check `_sources.yml` (live registry) against per-table `catalog/*.yml`.

The build code only reads `_sources.yml` to resolve a table's on-disk
location. The per-table `catalog/<name>.yml` is documentation. Nothing
prevents the two from drifting apart silently when a source file moves on
Nemo — this test fails fast when they do.

For every entry in `_sources.yml::sources` with a concrete `relative_path`
(not a per-tumour pattern), we build the fully-resolved HPC path and
require it to appear verbatim in the corresponding catalog's
`provenance.hpc_source`.
"""

from __future__ import annotations

import yaml

from tx_data.paths import REPO_ROOT

CATALOG_DIR = REPO_ROOT / "catalog"
SOURCES_YML = CATALOG_DIR / "_sources.yml"

SOURCE_TO_CATALOG = {
    "hg19_refgene__refgene": "hg19_refgene",
    "hg19_refgene__cytoband": "hg19_refgene",
    "kallisto_counts": "kallisto",
    "kallisto_tpm": "kallisto",
}


def _resolve(entry: dict, cohorts: dict, alt_roots: dict[str, str]) -> str | None:
    rel = entry.get("relative_path")
    if rel is None:
        return None
    root_alias = entry.get("source_root")
    if root_alias is None:
        root = cohorts[entry["cohort"]]["source_root_default"]
    elif root_alias == "in_repo":
        # Static reference data committed directly in tx_data; catalog
        # documents the upstream origin (URLs) rather than a mirrored path,
        # so skip equality-vs-catalog here.
        return None
    elif root_alias in alt_roots:
        root = alt_roots[root_alias]
    else:
        raise AssertionError(
            f"source_root alias {root_alias!r} missing from alt_source_roots"
        )
    return f"{root.rstrip('/')}/{rel}"


def test_sources_vs_catalog_paths_agree():
    sources_doc = yaml.safe_load(SOURCES_YML.read_text())
    cohorts = sources_doc["cohorts"]
    alt_roots = sources_doc.get("alt_source_roots", {}) or {}
    sources = sources_doc.get("sources", {}) or {}

    mismatches: list[str] = []

    for src_name, entry in sources.items():
        expected = _resolve(entry, cohorts, alt_roots)
        if expected is None:
            # Per-tumour pattern or no relative_path — skip.
            continue

        catalog_name = SOURCE_TO_CATALOG.get(src_name, src_name)
        catalog_path = CATALOG_DIR / f"{catalog_name}.yml"
        if not catalog_path.is_file():
            mismatches.append(
                f"{src_name}: expected catalog file {catalog_path} does not exist"
            )
            continue

        catalog_doc = yaml.safe_load(catalog_path.read_text())
        hpc_source = (catalog_doc.get("provenance") or {}).get("hpc_source") or ""
        if not hpc_source:
            # Catalog deliberately leaves hpc_source blank — fine.
            continue
        if "http://" in hpc_source or "https://" in hpc_source:
            # Catalog documents upstream URLs (e.g. UCSC downloads) rather than
            # a Nemo cache path — the authoritative source lives upstream, so
            # a mismatch against _sources.yml's cache location is expected.
            continue
        if expected not in hpc_source:
            mismatches.append(
                f"{src_name}: _sources.yml resolves to\n    {expected}\n"
                f"  but catalog/{catalog_name}.yml::provenance.hpc_source is\n"
                f"    {hpc_source.strip()}"
            )

    assert not mismatches, "Catalog ↔ _sources.yml drift:\n\n" + "\n\n".join(mismatches)


def test_every_source_declares_a_known_cohort():
    sources_doc = yaml.safe_load(SOURCES_YML.read_text())
    known = set(sources_doc["cohorts"]) | {"reference"}
    bad = {
        name: entry.get("cohort")
        for name, entry in sources_doc["sources"].items()
        if entry.get("cohort") not in known
    }
    assert not bad, f"sources with missing/unknown cohort: {bad}"


def test_catalog_cohort_matches_sources():
    sources_doc = yaml.safe_load(SOURCES_YML.read_text())
    mismatches = []
    for src_name, entry in sources_doc["sources"].items():
        catalog_name = SOURCE_TO_CATALOG.get(src_name, src_name)
        catalog_doc = yaml.safe_load((CATALOG_DIR / f"{catalog_name}.yml").read_text())
        if catalog_doc.get("cohort") != entry.get("cohort"):
            mismatches.append(
                f"{src_name}: _sources.yml cohort={entry.get('cohort')!r}, "
                f"catalog/{catalog_name}.yml cohort={catalog_doc.get('cohort')!r}"
            )
    assert not mismatches, "\n".join(mismatches)


def test_every_catalog_declares_cohort():
    sources_doc = yaml.safe_load(SOURCES_YML.read_text())
    known = set(sources_doc["cohorts"]) | {"reference"}
    missing = [
        p.name for p in sorted(CATALOG_DIR.glob("*.yml"))
        if not p.name.startswith("_")
        and yaml.safe_load(p.read_text()).get("cohort") not in known
    ]
    assert not missing, f"catalog files without a valid cohort: {missing}"


def test_wes_root_constant_matches_registry():
    from tx_data.paths import WES_RELEASE_ROOT

    sources_doc = yaml.safe_load(SOURCES_YML.read_text())
    assert str(WES_RELEASE_ROOT) == sources_doc["cohorts"]["wes"]["source_root_default"]


def test_every_nemo_source_exists_in_mock(monkeypatch):
    """nemo_mock/ mirrors /nemo: each non-pattern source must exist in it."""
    from tx_data.paths import DEFAULT_MOCK_ROOT
    from tx_data.sources import resolve_source

    monkeypatch.delenv("TX_NEMO_ROOT", raising=False)
    sources_doc = yaml.safe_load(SOURCES_YML.read_text())
    missing = []
    for name, entry in sources_doc["sources"].items():
        if "relative_path" not in entry:
            continue
        p = resolve_source(name)
        if entry.get("source_root") != "in_repo":
            assert DEFAULT_MOCK_ROOT in p.parents, (name, p)
        if not p.exists():
            missing.append(f"{name}: {p}")
    assert not missing, "missing from nemo_mock:\n" + "\n".join(missing)
