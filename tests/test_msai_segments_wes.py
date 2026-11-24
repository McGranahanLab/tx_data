"""msai_segments (WES, from the ALPACA input) on the committed mock (WGD_pipeline config.wes.toml, 2026-11-24)."""

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
    yield BUILDERS["msai_segments"](), BUILDERS["conipher_mutations"](), cp
    mp.undo()


def test_catalog_columns_and_tumours(built):
    seg, _, cp = built
    assert [c["name"] for c in catalog.load_catalog("msai_segments")["columns"]] == seg.columns
    assert set(seg["tumour_id"]) == set(cp["tumour_id"])
    assert seg["chr"].str.starts_with("chr").sum() == 0 and (seg["end"] >= seg["start"]).all()
    assert (seg["cn_major"] >= seg["cn_minor"]).all()


def test_samples_match_conipher_mutations(built):
    """WGD step 00 stops unless SCoutput, treeTable and the segments have the same samples."""
    seg, mut, _ = built
    for t in set(seg["tumour_id"]):
        a = set(seg.filter(pl.col("tumour_id") == t)["sample"])
        assert a == set(mut.filter(pl.col("tumour_id") == t)["sample"]), t


def test_clonal_gd_tumours_are_doubled(built):
    """Every sample of a clonal-GD tumour has > 50% of its genome at major CN >= 2 (step 00's rule)."""
    seg, _, _ = built
    frac = (seg.with_columns((pl.col("end") - pl.col("start") + 1).alias("len"))
               .group_by("tumour_id", "sample")
               .agg((((pl.col("cn_major").round(0) >= 2) * pl.col("len")).sum() / pl.col("len").sum()).alias("f")))
    for t in MANIFEST:
        g = frac.filter(pl.col("tumour_id") == t["tumour_id_canonical"])
        if g.height and t["clones"][0]["gds_at_clone"] > 0:
            assert g["f"].min() > 0.5, t["tumour_id_canonical"]


def test_tumour_id_mismatch_and_columns(tmp_path):
    from tx_data.builds import msai_segments as ms

    f = tmp_path / "ALPACA_input_table.csv"
    f.write_text("tumour_id,segment,sample,cpnA,cpnB\nLTX0002-Tumour1,1_10_20,S.R1..h,2.1,0.9\n")
    with pytest.raises(ValueError, match="does not match the folder"):
        ms.read_one("LTX0001-Tumour1", f)
    row = ms.read_one("LTX0002-Tumour1", f).row(0, named=True)
    assert (row["chr"], row["start"], row["end"], row["sample"], row["cn_major"], row["cn_minor"]) == \
        ("1", 10, 20, "S-R1--h", 2.1, 0.9)
    f.write_text("segment,sample,cpn_A,cpn_B\n1_10_20,S,1,1\n")
    with pytest.raises(ValueError, match="missing required columns"):
        ms.read_one("LTX0002-Tumour1", f)


def test_alpaca_symlinks_are_logged(tmp_path, capsys):
    from tx_data.builds import msai_segments as ms

    (tmp_path / "runs" / "v2" / "input").mkdir(parents=True)
    (tmp_path / "runs" / "v1" / "output").mkdir(parents=True)
    (tmp_path / "ALPACA").mkdir()
    (tmp_path / "ALPACA" / "input").symlink_to(tmp_path / "runs" / "v2" / "input")
    (tmp_path / "ALPACA" / "output").symlink_to(tmp_path / "runs" / "v1" / "output")
    ms.log_alpaca_links(tmp_path / "ALPACA" / "input")
    out = capsys.readouterr().out
    assert "runs/v2/input" in out and "runs/v1/output" in out and out.count("is a symlink") == 2
