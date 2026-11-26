"""wgd_calls_pipeline (WGD_pipeline WES output) on the mock, and its source-dir parameter (2026-11-26)."""

from __future__ import annotations

import polars as pl
import pytest

from tx_data import catalog
from tx_data.builds import wgd_calls_pipeline as w


@pytest.fixture(scope="module")
def built():
    mp = pytest.MonkeyPatch()
    mp.delenv("TX_NEMO_ROOT", raising=False)
    mp.delenv(w.ENV_DIR, raising=False)
    yield w.build()
    mp.undo()


def test_columns_match_catalog_and_wgs(built):
    assert [c["name"] for c in catalog.load_catalog("wgd_calls_pipeline")["columns"]] == built.columns
    assert built.columns == [c["name"] for c in catalog.load_catalog("wgs_wgd_calls")["columns"]]


def test_calls_and_derived(built):
    assert built["tumour_id"].n_unique() == 9 and set(built["status"]) <= {"resolved", "unresolved"}
    assert built["cumulative_gds_to_clone"].null_count() == 0
    assert built.filter(pl.col("quality_pass") == False)["number_of_gds_at_clone"].min() >= 1  # noqa: E712


def test_source_dir_and_env(tmp_path, monkeypatch, built):
    src = w._source(None, w.TABLE, w.ENV_DIR, "wes")
    (tmp_path / "wgd_calls.tsv").write_bytes(src.read_bytes())
    monkeypatch.setenv(w.ENV_DIR, str(tmp_path))
    assert w.build().height == built.height
    monkeypatch.delenv(w.ENV_DIR)
    monkeypatch.setenv("TX_NEMO_ROOT", str(tmp_path / "nowhere"))
    with pytest.raises(FileNotFoundError, match="TX_WES_WGD_CALLS_DIR=<darwins_analysis>/WGD_pipeline/.*output/wes"):
        w.build()
