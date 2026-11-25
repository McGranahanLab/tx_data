"""wgs_wgd_calls (WGD_pipeline WGS output) on the mock, and the --source-dir parameter."""

from __future__ import annotations

import polars as pl
import pytest

from tx_data import catalog
from tx_data.builds import wgs_wgd_calls


@pytest.fixture(scope="module")
def built():
    mp = pytest.MonkeyPatch()
    mp.delenv("TX_NEMO_ROOT", raising=False)
    mp.delenv(wgs_wgd_calls.ENV_DIR, raising=False)
    yield wgs_wgd_calls.build()
    mp.undo()


def test_columns_match_catalog(built):
    assert [c["name"] for c in catalog.load_catalog("wgs_wgd_calls")["columns"]] == built.columns


def test_derived_null_without_calls(built):
    nc = built.filter(pl.col("status").is_in(["qc_failed", "no_alpaca_tree"]))
    assert nc.height and nc["cumulative_gds_to_clone"].null_count() == nc.height
    calls = built.filter(~pl.col("status").is_in(["qc_failed", "no_alpaca_tree"]))
    assert calls["cumulative_gds_to_clone"].null_count() == 0
    assert (calls["has_wgd_ancestor"] == (calls["cumulative_gds_to_clone"] > 0)).all()


def test_quality_pass_kept_and_all_calls_present(built):
    assert built.filter(pl.col("quality_pass") == False).height >= 1          # noqa: E712  a GD below the cut-off
    assert built.filter(pl.col("quality_pass") == False)["number_of_gds_at_clone"].min() >= 1  # noqa: E712


def test_source_dir_parameter(tmp_path, built):
    src = wgs_wgd_calls._source(None)
    (tmp_path / "wgd_calls.tsv").write_bytes(src.read_bytes())
    assert wgs_wgd_calls._source(tmp_path) == tmp_path / "wgd_calls.tsv"
    bad = pl.read_csv(src, separator="\t").drop("quality_pass")
    bad.write_csv(tmp_path / "wgd_calls.tsv", separator="\t")
    with pytest.raises(ValueError, match="quality_pass"):
        wgs_wgd_calls.build(tmp_path)
    wgs_wgd_calls.build()  # restore data/ from the mock


def test_missing_default_says_what_to_set(tmp_path, monkeypatch):
    """Nemo 2026-11-25: on prod the mock-layout default doesn't exist; the error names TX_WGS_WGD_CALLS_DIR."""
    from tx_data.builds import wgs_wgd_calls as w

    monkeypatch.delenv(w.ENV_DIR, raising=False)
    monkeypatch.setenv("TX_NEMO_ROOT", str(tmp_path))
    with pytest.raises(FileNotFoundError, match="TX_WGS_WGD_CALLS_DIR=<darwins_analysis>/WGD_pipeline"):
        w._source(None)
