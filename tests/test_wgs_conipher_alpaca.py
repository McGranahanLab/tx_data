"""wgs_conipher_* and wgs_alpaca* on the WGS mock (inputs of the WGS WGD pipeline, 2026-11-17)."""

from __future__ import annotations

import json
import runpy

import polars as pl
import pytest

from tx_data import catalog
from tx_data.paths import DEFAULT_MOCK_ROOT, REPO_ROOT

MANIFEST = DEFAULT_MOCK_ROOT / "_cohort_manifest_wgs.json"
TABLES = ("wgs_conipher_mutations", "wgs_conipher_tree_edges", "wgs_conipher_cluster_info", "wgs_alpaca",
          "wgs_alpaca_wgd_ratios")


@pytest.fixture(scope="module")
def built():
    mp = pytest.MonkeyPatch()
    mp.delenv("TX_NEMO_ROOT", raising=False)
    runpy.run_path(str(REPO_ROOT / "scripts" / "generate_mocks_wgs.py"), run_name="__main__")
    from tx_data.builds import BUILDERS

    BUILDERS["wgs_battenberg_inventory"]()
    out = {t: BUILDERS[t]() for t in TABLES}
    out["manifest"] = json.loads(MANIFEST.read_text())
    yield out
    mp.undo()


def test_every_tumour_in_every_table(built):
    ids = {t["tumour_id_canonical"] for t in built["manifest"]["tumours"]}
    for t in TABLES:
        assert set(built[t]["tumour_id"]) == ids, t


def test_catalog_columns_match(built):
    for t in TABLES:
        assert [c["name"] for c in catalog.load_catalog(t)["columns"]] == built[t].columns, t


def test_tree_edges_match_alpaca_tree(built):
    e = built["wgs_conipher_tree_edges"]
    for t in built["manifest"]["tumours"]:
        g = e.filter(pl.col("tumour_id") == t["tumour_id_canonical"])
        clones = {c for p in t["tree_paths"] for c in p}
        if len(clones) == 1:
            assert g.select("parent_cluster", "child_cluster").rows() == [(1, 1)]   # CONIPHER's single-clone form
        else:
            assert {f"clone{c}" for c in g["child_cluster"]} | {"clone1"} == clones


def test_mutations_samples_and_tree_cluster(built):
    m = built["wgs_conipher_mutations"]
    assert m["sample_name_hash"].null_count() == 0 and m["chr"].str.starts_with("chr").sum() == 0
    off = m.filter(pl.col("tree_cluster").is_null())
    assert set(off["tumour_id"]) == {"LTX0103-Tumour1"} and set(off["cluster"]) == {9}
    trunk = m.filter((pl.col("tumour_id") == "LTX0004-Tumour1") & (pl.col("cluster") == 1))
    assert trunk["mut_copy"].max() == 2.0                                      # WGD-like: truncal on 2 copies
    seg = pl.read_parquet(REPO_ROOT / "data" / "wgs_consensus_segments.parquet") if (
        REPO_ROOT / "data" / "wgs_consensus_segments.parquet").exists() else None
    if seg is not None:  # step 00 needs the same SAMPLE names in mutations and CN segments
        for (t,), g in m.group_by("tumour_id"):
            assert set(g["sample"]) == set(seg.filter(pl.col("tumour_id") == t)["sample"]), t


def test_cluster_info_ccf_qc_inputs(built):
    ci = built["wgs_conipher_cluster_info"]
    trunk = ci.filter(pl.col("truncal") & pl.col("tree_clust"))
    low = trunk.group_by("tumour_id").agg(pl.col("mean_ccf").min())
    assert low.filter(pl.col("mean_ccf") <= 90)["tumour_id"].to_list() == ["LTX0104-Tumour2"]
    assert ci.filter(pl.col("tree_clust").is_null())["cluster_id"].unique().to_list() == [9]


def test_alpaca_and_ratios(built):
    a, r = built["wgs_alpaca"], built["wgs_alpaca_wgd_ratios"]
    assert a["segment"].str.contains(r"^\d+_\d+_\d+$").all()
    assert set(a.select("tumour_id", "clone").unique().rows()) == set(r.select("tumour_id", "clone").rows())
    assert r.filter(pl.col("is_mrca"))["tumour_id"].n_unique() == r["tumour_id"].n_unique()
    wgd = r.filter((pl.col("clone") == "clone1") & (pl.max_horizontal("ratio_A", "ratio_B") > 1.5))
    assert {"LTX0001-Tumour1", "LTX0004-Tumour1"} <= set(wgd["tumour_id"])       # the WGD-like mock tumours


def test_manual_qc_mock_covers_wgs(built):
    qc = pl.read_csv(DEFAULT_MOCK_ROOT / "project" / "proj-tracerx-lung" / "tracerx" / "_RELEASE" / "release_tx842"
                     / "_aggregate" / "20250729_tx842_manual_qc_complete_post_refphase.tsv", separator="\t")
    ids = {f"{p}-Tumour{n}" for p, n in qc.select("patient", "tumour_id").iter_rows()}
    assert {t["tumour_id_canonical"] for t in built["manifest"]["tumours"]} <= ids
    assert qc.filter(~pl.col("FINALIZED"))["patient"].to_list() == ["LTX0105"]
