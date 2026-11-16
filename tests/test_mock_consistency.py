"""Cross-source agreement of the committed WES mock (nemo_mock/).

Two invariants, checked on the raw mock files (what the builds and the
downstream analyses read):

* Tree clones. For every tumour with a tree, the clones are the same in
  ALPACA (``cp_table.csv`` + ``all_tumours_combined.csv``), the CONIPHER RDS
  (``ccf_cluster_table`` / ``clonality_table_corrected`` rows) and
  ``clone_info`` (``treeClones == TRUE``). Non-tree clusters are the same in
  ``clone_info`` (``treeClones == FALSE``), the RDS ``tree_removed_clusters``
  and the PGDD scenario truth (``in_tree = FALSE``); PGDD's on-tree clusters
  are the tree clones.
* Sample hashes. Each DNA region has one ``sample_name_hash`` across the
  muttable, purity_ploidy, cp_table, the CONIPHER RDS (R-style ``.``/``..``
  separators) and the PGDD outputs, and every one is in the sample-hash map.

The cohort manifest is the generators' shared source; it is checked too.
Regenerate with ``pixi run mocks`` + the two R generators if this fails.
"""

from __future__ import annotations

import gzip
import json
import shutil
import subprocess
from pathlib import Path

import polars as pl
import pytest

from tx_data.paths import DEFAULT_MOCK_ROOT, REPO_ROOT

WES = DEFAULT_MOCK_ROOT / "project/proj-tracerx-lung/tctProjects/lungTx/Tx842/release/tracerx_842"
MANIFEST = DEFAULT_MOCK_ROOT / "_cohort_manifest.json"
PGDD_DIR = WES / "WGD/release_tx842/pgdd"
PGDD_TRUTH = REPO_ROOT / "synthetic/pgdd/pgdd_mock_scenarios.json"


def _canon(tumour: str) -> str:
    """LTX0001_1 / LTX0001_tumour1 / LTX0001-Tumour1 -> LTX0001-Tumour1."""
    if "-Tumour" in tumour:
        return tumour
    pid, _, n = tumour.replace("_tumour", "_").rpartition("_")
    return f"{pid}-Tumour{n}"


def _from_r_colname(col: str) -> str:
    """CONIPHER region column (LTX0001_SU_T1.R1..abc) -> muttable form (LTX0001_SU_T1-R1--abc)."""
    return col.replace(".", "-")


# ---- per-source readers: {tumour: set(...)} --------------------------------

@pytest.fixture(scope="module")
def manifest() -> dict:
    return {t["tumour_id_canonical"]: t for t in json.loads(MANIFEST.read_text())["tumours"]}


@pytest.fixture(scope="module")
def tree_tumours(manifest) -> list[str]:
    return sorted(t for t, m in manifest.items() if m["has_tree"])


@pytest.fixture(scope="module")
def rds() -> dict[str, dict]:
    rscript = shutil.which("Rscript")
    if rscript is None:
        pytest.skip("Rscript not on PATH (run via `pixi run test`)")
    expr = (
        "x <- readRDS(commandArgs(TRUE)[1]); for (b in x) cat(paste(c("
        "b$parameters$sampleID,"
        "paste(rownames(b$nested_pyclone$ccf_cluster_table), collapse=','),"
        "paste(rownames(b$clonality_out$clonality_table_corrected), collapse=','),"
        "paste(b$tree_removed_clusters, collapse=','),"
        "paste(colnames(b$nested_pyclone$ccf_cluster_table), collapse=','),"
        "paste(colnames(b$clonality_out$clonality_table_corrected), collapse=',')"
        "), collapse='\\t'), '\\n')"
    )
    out = subprocess.run([rscript, "-e", expr, str(WES / "tx842_trees.RDS")],
                         check=True, capture_output=True, text=True).stdout
    res = {}
    for line in out.strip().splitlines():
        tid, ccf_rows, cln_rows, removed, ccf_cols, cln_cols = (f.strip() for f in line.split("\t"))
        split = lambda s: set(s.split(",")) if s else set()  # noqa: E731
        res[_canon(tid)] = dict(
            ccf_clones={f"clone{c}" for c in split(ccf_rows)},
            clonality_clones={f"clone{c}" for c in split(cln_rows)},
            removed=split(removed),
            ccf_samples={_from_r_colname(c) for c in split(ccf_cols)},
            clonality_samples={_from_r_colname(c) for c in split(cln_cols)},
        )
    return res


@pytest.fixture(scope="module")
def cp_tables() -> dict[str, pl.DataFrame]:
    return {p.parent.name: pl.read_csv(p) for p in sorted((WES / "ALPACA/input").glob("*/cp_table.csv"))}


@pytest.fixture(scope="module")
def clone_info() -> pl.DataFrame:
    df = pl.read_csv(WES / "primary_met/tx842_primary_met_cloneInfo.df.txt", separator="\t",
                     infer_schema=False)
    return df.with_columns(
        pl.col("tumour_id").map_elements(_canon, return_dtype=pl.String),
        pl.col("treeClones").str.to_uppercase() == "TRUE",
    )


@pytest.fixture(scope="module")
def muttable_samples() -> dict[str, set[str]]:
    with gzip.open(WES / "tx842_mutation_table.tsv.gz", "rt") as f:
        df = pl.read_csv(f.read().encode(), separator="\t", columns=["patient_tumour", "sample_name_hash"])
    return {_canon(t): set(g["sample_name_hash"]) for (t,), g in df.unique().group_by("patient_tumour")}


# ---- clone sets -------------------------------------------------------------

def test_tree_clones_agree_across_sources(manifest, tree_tumours, rds, cp_tables, clone_info):
    alpaca_out = pl.read_csv(WES / "ALPACA/output/cohort_results/all_tumours_combined.csv")
    ci_tree = clone_info.filter(pl.col("treeClones"))
    assert set(cp_tables) == set(rds) == set(tree_tumours) == set(alpaca_out["tumour_id"]) \
        == set(clone_info["tumour_id"])
    bad = []
    for t in tree_tumours:
        sets = {
            "manifest": {c["name"] for c in manifest[t]["clones"]},
            "cp_table": set(cp_tables[t]["clone"]),
            "alpaca_output": set(alpaca_out.filter(pl.col("tumour_id") == t)["clone"]),
            "rds_ccf_cluster_table": rds[t]["ccf_clones"],
            "rds_clonality_table": rds[t]["clonality_clones"],
            "clone_info_tree": {f"clone{c}" for c in ci_tree.filter(pl.col("tumour_id") == t)["mutation_cluster"]},
        }
        if len({frozenset(s) for s in sets.values()}) != 1:
            bad.append(f"{t}: {sets}")
    assert not bad, "tree clones disagree:\n" + "\n".join(bad)


def test_non_tree_clusters_agree_across_sources(manifest, tree_tumours, rds, clone_info):
    truth = {t["tumour_id"]: t for t in json.loads(PGDD_TRUTH.read_text())["tumours"]}
    ci_off = clone_info.filter(~pl.col("treeClones"))
    bad = []
    for t in tree_tumours:
        sets = {
            "manifest": {str(k) for k in manifest[t]["non_tree_clusters"]},
            "rds_tree_removed_clusters": rds[t]["removed"],
            "clone_info_non_tree": set(ci_off.filter(pl.col("tumour_id") == t)["mutation_cluster"]),
            "pgdd_in_tree_false": {c["cluster_id"] for c in truth[t]["clusters"] if not c["in_tree"]},
        }
        if len({frozenset(s) for s in sets.values()}) != 1:
            bad.append(f"{t}: {sets}")
        pgdd_on = {f"clone{c['cluster_id']}" for c in truth[t]["clusters"] if c["in_tree"]}
        if pgdd_on != rds[t]["ccf_clones"]:
            bad.append(f"{t}: PGDD on-tree clusters {pgdd_on} != tree clones {rds[t]['ccf_clones']}")
    assert not bad, "non-tree clusters disagree:\n" + "\n".join(bad)
    # The deliberate off-tree PGDD clusters (generate_mock_pgdd.R scenarios).
    assert {t: m["non_tree_clusters"] for t, m in manifest.items() if m["non_tree_clusters"]} == {
        "LTX0005-Tumour1": [6], "LTX0006-Tumour1": [7], "LTX0007-Tumour1": [12],
    }
    # Tumours without a tree: PGDD has calls, nothing else has clones.
    no_tree = {t for t, m in manifest.items() if not m["has_tree"]}
    assert no_tree == {"LTX0009-Tumour1", "LTX0010-Tumour1"}
    for t in no_tree:
        assert not any(c["in_tree"] for c in truth[t]["clusters"])


def test_muttable_clusters_are_tree_clones(manifest, tree_tumours):
    with gzip.open(WES / "tx842_mutation_table.tsv.gz", "rt") as f:
        df = pl.read_csv(f.read().encode(), separator="\t", columns=["patient_tumour", "mutation_cluster"])
    for t in tree_tumours:
        got = {f"clone{int(c)}" for c in df.filter(
            pl.col("patient_tumour").map_elements(_canon, return_dtype=pl.String) == t)["mutation_cluster"]}
        assert got == {c["name"] for c in manifest[t]["clones"]}, t


# ---- sample hashes ----------------------------------------------------------

def test_one_hash_per_region_across_sources(manifest, rds, cp_tables, muttable_samples):
    purity = pl.read_csv(WES / "tx842_purity_ploidy.csv")
    pgdd_regions = pl.read_csv(PGDD_DIR / "pgdd_GDs_per_region.tsv", separator="\t")
    pgdd_input = pl.read_csv(REPO_ROOT / "synthetic/pgdd/pgdd_input_mock.tsv", separator="\t",
                             columns=["tumour_id", "sample_id"])
    bad = []
    for t, m in manifest.items():
        expected = {s["sample_name_hash"] for s in m["samples"]}
        assert len(expected) == m["n_samples"]
        sets = {
            "muttable": muttable_samples[t],
            "pgdd_GDs_per_region": set(pgdd_regions.filter(pl.col("tumour_id") == t)["sample_id"]),
            "pgdd_input": set(pgdd_input.filter(pl.col("tumour_id") == t)["sample_id"]),
        }
        if m["has_tree"]:
            prefix = f"{m['patient_id']}_SU_T{m['tumour_ordinal']}-"
            sets |= {
                "purity_ploidy": {s for s in purity["sample_name_hash"] if s.startswith(prefix)},
                "cp_table": {c for c in cp_tables[t].columns if c != "clone"},
                "rds_ccf_cluster_table": rds[t]["ccf_samples"],
                "rds_clonality_table": rds[t]["clonality_samples"],
            }
        for name, got in sets.items():
            if got != expected:
                bad.append(f"{t} {name}: {sorted(got ^ expected)}")
    assert not bad, "sample hashes disagree:\n" + "\n".join(bad)


def test_pgdd_event_samples_use_region_hashes(manifest):
    known = {s["sample_name_hash"] for m in manifest.values() for s in m["samples"]}
    events = pl.read_csv(PGDD_DIR / "pgdd_GD_events.tsv", separator="\t")
    used = {s for row in events["samples"].drop_nulls() for s in row.split(",")}
    assert used and used <= known, sorted(used - known)


def test_region_hashes_in_sample_map_once(manifest):
    smap = pl.read_csv(WES / "mapping_files/tx842_samplehash_tumour_lesion_map.tsv", separator="\t")
    assert smap["sample_hash"].is_unique().all()
    by_hash = {r["sample_hash"]: r for r in smap.to_dicts()}
    for t, m in manifest.items():
        for s in m["samples"]:
            h = s["sample_name_hash"].rsplit("--", 1)[1]
            assert h in by_hash, (t, s)
            assert _canon(by_hash[h]["patient_tumour"]) == t
    # One DNA hash per region: no other `_SU_T<n>-R<r>` hash for the same region
    # outside the RNA (kallisto) columns.
    kallisto = set(pl.read_parquet(
        DEFAULT_MOCK_ROOT / "project/proj-tracerx-lung/tctProjects/lungTx/Tx842/release/"
        "tx842_mets_rnaseq/_aggregate/kallisto_mat/kallisto_tpm.parquet").columns)
    # Normals (`_SU_N<n>`) are not in the map by design (catalog/sample_hash_map.yml).
    rna_hashes = {c.rsplit("--", 1)[1] for c in kallisto if "--" in c and "_SU_N" not in c}
    dna = {s["sample_name_hash"].rsplit("--", 1)[1] for m in manifest.values() for s in m["samples"]}
    assert set(smap["sample_hash"]) == dna | rna_hashes
