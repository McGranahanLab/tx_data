"""sample_clone_counts (WES), wgs_sample_clone_counts (WGS) and their cohort_metadata summary (2026-11-25)."""

from __future__ import annotations

import json

import polars as pl
import pytest

from tx_data import catalog
from tx_data.paths import DEFAULT_MOCK_ROOT


@pytest.fixture(scope="module")
def built():
    mp = pytest.MonkeyPatch()
    mp.delenv("TX_NEMO_ROOT", raising=False)
    from tx_data.builds import BUILDERS

    for t in ("clone_proportions", "conipher_mutations", "conipher_tree_edges", "wgs_battenberg_inventory",
              "wgs_clone_proportions", "wgs_conipher_mutations", "wgs_conipher_tree_edges"):
        BUILDERS[t]()
    out = {t: BUILDERS[t]() for t in ("sample_clone_counts", "wgs_sample_clone_counts", "cohort_metadata")}
    yield out
    mp.undo()


def test_catalog_columns(built):
    for t in ("sample_clone_counts", "wgs_sample_clone_counts"):
        assert [c["name"] for c in catalog.load_catalog(t)["columns"]] == built[t].columns, t


def test_wes_counts_match_manifest(built):
    """cp_table tumours: samples and clones as in the mock manifest; tumours without ALPACA from CONIPHER."""
    df = built["sample_clone_counts"]
    for t in json.loads((DEFAULT_MOCK_ROOT / "_cohort_manifest.json").read_text())["tumours"]:
        r = df.filter(pl.col("tumour_id") == t["tumour_id_canonical"]).row(0, named=True)
        assert r["n_samples"] == t["n_samples"], t["tumour_id_canonical"]
        if r["in_cp_table"]:
            assert r["source"] == "cp_table" and r["n_clones"] == len(t["clones"])
        else:   # no ALPACA input and no CONIPHER tree in the mock
            assert (r["source"], r["n_clones"], r["category"]) == ("conipher", None, "clones_unknown")


def test_wgs_single_sample_and_clone(built):
    df = built["wgs_sample_clone_counts"]
    m = json.loads((DEFAULT_MOCK_ROOT / "_cohort_manifest_wgs.json").read_text())["tumours"]
    for t in m:
        r = df.filter(pl.col("tumour_id") == t["tumour_id_canonical"]).row(0, named=True)
        clones = {c for p in t["tree_paths"] for c in p}
        assert r["single_sample"] == (len(t["samples"]) == 1), t["tumour_id_canonical"]
        assert r["single_clone"] == (len(clones) == 1), t["tumour_id_canonical"]
    assert set(df["category"]) <= {"single_sample_single_clone", "single_sample", "single_clone",
                                   "multi_sample_multi_clone"}


def test_cohort_metadata_summary(built):
    md, wgs = built["cohort_metadata"], built["wgs_sample_clone_counts"]
    row = md.filter((pl.col("cohort") == "wgs") & (pl.col("field") == "single_sample_single_clone_tumours")).row(
        0, named=True)
    ids = wgs.filter(pl.col("category") == "single_sample_single_clone")["tumour_id"].sort().to_list()
    assert row["value"] == str(len(ids)) and all(i in row["note"] for i in ids)
    wes = md.filter(pl.col("cohort") == "wes")
    assert set(wes["field"]) >= {"single_sample_tumours", "single_clone_tumours",
                                 "single_sample_single_clone_tumours", "clones_unknown_tumours",
                                 "single_sample_or_clone"}
