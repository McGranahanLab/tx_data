"""wgs_consensus_segments + wgs_clone_proportions on the WGS mock.

Self-contained: regenerates the WGS mock and builds the inventory first,
since both tables discover tumours and samples from it.
"""

from __future__ import annotations

import json
import runpy

import polars as pl
import pytest

from tx_data.paths import DEFAULT_MOCK_ROOT, REPO_ROOT

MANIFEST = DEFAULT_MOCK_ROOT / "_cohort_manifest_wgs.json"


@pytest.fixture(scope="module")
def built():
    mp = pytest.MonkeyPatch()
    mp.delenv("TX_NEMO_ROOT", raising=False)
    mp.delenv("TX_ID_CONVENTION", raising=False)
    runpy.run_path(str(REPO_ROOT / "scripts" / "generate_mocks_wgs.py"), run_name="__main__")
    from tx_data.builds import BUILDERS

    BUILDERS["wgs_battenberg_inventory"]()
    out = {
        "seg": BUILDERS["wgs_consensus_segments"](),
        "cp": BUILDERS["wgs_clone_proportions"](),
        "manifest": json.loads(MANIFEST.read_text()),
    }
    yield out
    mp.undo()


def test_segments_cover_every_tumour_and_sample(built):
    seg, m = built["seg"], built["manifest"]
    assert seg["sample_name_hash"].null_count() == 0
    for t in m["tumours"]:
        sub = seg.filter(pl.col("tumour_id") == t["tumour_id_canonical"])
        assert sorted(sub["sample"].unique()) == sorted(t["segments_sample_names"])
        assert sub.height == t["n_segments"] * len(t["samples"])
        assert set(sub["sample_name_hash"]) == {s["sample_name_hash"] for s in t["samples"]}


def test_segments_schema_and_values(built):
    seg = built["seg"]
    assert seg.schema["cn_major"] == pl.Float64 and seg.schema["start"] == pl.Int64
    assert not seg["chr"].str.starts_with("chr").any()
    assert (seg["end"] >= seg["start"]).all()
    assert ((seg["cn_total"] - seg["cn_major"] - seg["cn_minor"]).abs() < 1e-12).all()
    # Values are fractional, not rounded.
    assert ((seg["cn_major"] - seg["cn_major"].round()).abs() > 1e-6).mean() > 0.9


def test_legacy_a_b_header_read(built):
    seg, m = built["seg"], built["manifest"]
    legacy = [t for t in m["tumours"] if t["segments_header"] == "legacy_A_B"]
    assert legacy
    for t in legacy:
        assert seg.filter(pl.col("wgs_tumour_label") == t["wgs_tumour_label"])["cn_major"].null_count() == 0


def test_multi_tumour_labels_standardised(built):
    seg = built["seg"]
    sub = seg.filter(pl.col("patient_id") == "LTX0104").select(
        "wgs_tumour_label", "tumour_id").unique().sort("tumour_id")
    assert sub.rows() == [("LTX0104_tumour1", "LTX0104-Tumour1"),
                          ("LTX0104_tumour2", "LTX0104-Tumour2")]


def test_clone_proportions_match_manifest(built):
    cp, m = built["cp"], built["manifest"]
    assert cp["sample_name_hash"].null_count() == 0
    for t in m["tumours"]:
        sub = cp.filter(pl.col("tumour_id") == t["tumour_id_canonical"])
        assert sub["clone"].n_unique() == t["n_clones"]
        for sample, props in t["clone_proportions"].items():
            got = dict(sub.filter(pl.col("sample_id") == sample).select("clone", "proportion").rows())
            assert got == pytest.approx(props)


def test_single_clone_tumours_identifiable(built):
    cp, m = built["cp"], built["manifest"]
    n = cp.group_by("tumour_id").agg(pl.col("clone").n_unique().alias("n"))
    single = set(n.filter(pl.col("n") == 1)["tumour_id"])
    assert single == {t["tumour_id_canonical"] for t in m["tumours"] if t["n_clones"] == 1}
