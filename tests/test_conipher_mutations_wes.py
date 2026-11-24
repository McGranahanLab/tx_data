"""conipher_mutations (WES) on the committed mock (release_tx842 tree; WGD_pipeline config.wes.toml, 2026-11-24)."""

from __future__ import annotations

import json

import polars as pl
import pytest

from tx_data import catalog
from tx_data.paths import DEFAULT_MOCK_ROOT

MANIFEST = json.loads((DEFAULT_MOCK_ROOT / "_cohort_manifest.json").read_text())["tumours"]


@pytest.fixture(scope="module")
def built():
    mp = pytest.MonkeyPatch()
    mp.delenv("TX_NEMO_ROOT", raising=False)
    from tx_data.builds import BUILDERS

    cp = BUILDERS["clone_proportions"]()
    yield BUILDERS["conipher_mutations"](), cp
    mp.undo()


def test_catalog_columns_match(built):
    m, _ = built
    assert [c["name"] for c in catalog.load_catalog("conipher_mutations")["columns"]] == m.columns


def test_every_manifest_tumour_present(built):
    m, _ = built
    assert set(m["tumour_id"]) == {t["tumour_id_canonical"] for t in MANIFEST}
    assert m["chr"].str.starts_with("chr").sum() == 0


def test_tree_cluster_off_tree_and_no_tree(built):
    m, _ = built
    for t in MANIFEST:
        g = m.filter(pl.col("tumour_id") == t["tumour_id_canonical"])
        off = set(g.filter(pl.col("tree_cluster").is_null())["cluster"])
        if t["has_tree"]:
            assert off == set(t["non_tree_clusters"]), t["tumour_id_canonical"]
            assert (g.drop_nulls("tree_cluster")["tree_cluster"] == g.drop_nulls("tree_cluster")["cluster"]).all()
        else:
            assert g["tree_cluster"].null_count() == g.height     # no treeTable.tsv: whole tumour null


def test_sample_converted_to_dash_form(built):
    """R-style SAMPLE (`.`/`..`) becomes the cp_table form (`-`/`--`) in both SCoutput and treeTable."""
    m, cp = built
    assert m["sample"].str.contains(r"\.").sum() == 0
    with_cp = set(cp["tumour_id"])
    pairs = set(m.filter(pl.col("tumour_id").is_in(with_cp)).select("tumour_id", "sample").unique().iter_rows())
    assert pairs == set(cp.filter(pl.col("tumour_id").is_in(set(m["tumour_id"])))
                          .select("tumour_id", "sample_id").unique().iter_rows())
    assert m.filter(pl.col("tumour_id") == "LTX0005-Tumour1", pl.col("cluster") != 6)["tree_cluster"].null_count() == 0


def test_clonal_gd_truncal_mut_copy(built):
    m, _ = built
    for t in MANIFEST:
        trunk = m.filter((pl.col("tumour_id") == t["tumour_id_canonical"]) & (pl.col("cluster") == 1))
        want = 2.0 if t["clones"] and t["clones"][0]["gds_at_clone"] > 0 else 1.0
        assert set(trunk["mut_copy"]) == {want}, t["tumour_id_canonical"]


def test_case_id_mismatch_fails(tmp_path, monkeypatch):
    from tx_data.builds import conipher_mutations as cm

    f = tmp_path / "LTX0001_1.SCoutput.CLEAN.tsv"
    f.write_text("\t".join(("CASE_ID", *cm.SOURCE)) + "\n" + "\t".join(("LTX0002_1", "S", "1", "10", "A", "C",
                                                                         "m", "1", "1", "1", "1")) + "\n")
    with pytest.raises(ValueError, match="CASE_ID"):
        cm.read_one("LTX0001", "1", f)
