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


def _resolve(entry: dict, cohort_root: str, alt_roots: dict[str, str]) -> str | None:
    rel = entry.get("relative_path")
    if rel is None:
        return None
    root_alias = entry.get("source_root")
    if root_alias is None:
        root = cohort_root
    elif root_alias in alt_roots:
        root = alt_roots[root_alias]
    else:
        raise AssertionError(
            f"source_root alias {root_alias!r} missing from alt_source_roots"
        )
    return f"{root.rstrip('/')}/{rel}"


def test_sources_vs_catalog_paths_agree():
    sources_doc = yaml.safe_load(SOURCES_YML.read_text())
    cohort_root = sources_doc["cohort"]["source_root_default"]
    alt_roots = sources_doc.get("alt_source_roots", {}) or {}
    sources = sources_doc.get("sources", {}) or {}

    mismatches: list[str] = []

    for src_name, entry in sources.items():
        expected = _resolve(entry, cohort_root, alt_roots)
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
