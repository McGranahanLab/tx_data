"""alpaca_trees / wgs_alpaca_trees, cohort_metadata, chr_table_hg19 / hg38 (2026-11-19)."""

from __future__ import annotations

import runpy

import polars as pl
import pytest
import yaml

from tx_data import catalog
from tx_data.paths import REPO_ROOT

TABLES = ("alpaca_trees", "wgs_alpaca_trees", "cohort_metadata", "chr_table_hg19", "chr_table_hg38")


@pytest.fixture(scope="module")
def built():
    mp = pytest.MonkeyPatch()
    mp.delenv("TX_NEMO_ROOT", raising=False)
    runpy.run_path(str(REPO_ROOT / "scripts" / "generate_mocks_wgs.py"), run_name="__main__")
    from tx_data.builds import BUILDERS

    BUILDERS["wgs_battenberg_inventory"]()
    out = {"alpaca_trees": BUILDERS["alpaca_trees"](), "wgs_alpaca_trees": BUILDERS["wgs_alpaca_trees"](),
           "cohort_metadata": BUILDERS["cohort_metadata"]()}
    chr_ = BUILDERS["chr_table"]()
    out |= {"chr_table_hg19": chr_["hg19"], "chr_table_hg38": chr_["hg38"]}
    out["wgs_ratios"] = BUILDERS["wgs_alpaca_wgd_ratios"]()
    out["wes_ratios"] = BUILDERS["alpaca_wgd_ratios"]()
    yield out
    mp.undo()


def test_catalog_columns_match(built):
    for t in TABLES:
        assert [c["name"] for c in catalog.load_catalog(t)["columns"]] == built[t].columns, t


def test_trees_agree_with_alpaca_ratio_parents(built):
    """The input tree and ALPACA's own output (wgd_ratio_scores clone -> parent) are the same tree."""
    for tree, ratios in (("alpaca_trees", "wes_ratios"), ("wgs_alpaca_trees", "wgs_ratios")):
        t = set(built[tree].select("tumour_id", "child", "parent").rows())
        r = set(built[ratios].select("tumour_id", "clone", "parent").rows())
        assert r <= t and {x[0] for x in t} >= {x[0] for x in r}, tree
        assert built[tree].group_by("tumour_id").agg((pl.col("parent") == "diploid").sum().alias("n"))["n"].to_list() \
            == [1] * built[tree]["tumour_id"].n_unique()                      # one MRCA per tumour


def test_chr_tables(built):
    for b, chr1, chrY in (("hg19", 249250621, 59373566), ("hg38", 248956422, 57227415)):
        d = built[f"chr_table_{b}"]
        assert d.columns == ["chr", "len"] and d.height == 24 and d["chr"].to_list()[0] == "1"
        assert dict(d.iter_rows())["1"] == chr1 and dict(d.iter_rows())["Y"] == chrY


def test_cohort_metadata_matches_registry(built):
    reg = yaml.safe_load((REPO_ROOT / "catalog" / "_sources.yml").read_text())["cohorts"]
    m = built["cohort_metadata"].filter(pl.col("field") == "genome_build")
    assert dict(m.select("cohort", "value").iter_rows()) == {c: reg[c]["genome_build"] for c in ("wes", "wgs")}
    assert dict(m.select("cohort", "value").iter_rows()) == {"wes": "hg19", "wgs": "hg38"}
