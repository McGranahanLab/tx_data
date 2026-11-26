"""wgs_alpaca_corrected / wgs_alpaca_wgd_ratios_corrected: WGS ALPACA correction in its own tables (2026-11-26)."""

from __future__ import annotations

import polars as pl
import pytest

from tx_data import catalog
from tx_data.builds import alpaca_corrected as ac
from tx_data.builds._base import canonical_output_path


@pytest.fixture(scope="module")
def built():
    mp = pytest.MonkeyPatch()
    mp.delenv("TX_NEMO_ROOT", raising=False)
    for e in (ac.ENV_DIR, ac.WGS_ENV_DIR):
        mp.delenv(e, raising=False)
    wes = ac.build()
    wgs = ac.build_wgs()
    yield wes, wgs
    mp.undo()


def test_catalogs_match(built):
    _, wgs = built
    assert [c["name"] for c in catalog.load_catalog("wgs_alpaca_corrected")["columns"]] == wgs.columns
    r = pl.read_parquet(canonical_output_path("wgs_alpaca_wgd_ratios_corrected"))
    assert [c["name"] for c in catalog.load_catalog("wgs_alpaca_wgd_ratios_corrected")["columns"]] == r.columns


def test_wes_and_wgs_kept_apart(built):
    wes, wgs = built
    assert pl.read_parquet(canonical_output_path("alpaca_corrected")).height == wes.height
    assert wgs.height
    wgs_ids = set(pl.read_parquet(canonical_output_path("wgs_alpaca"))["tumour_id"])
    assert set(wgs["tumour_id"]) <= wgs_ids


def test_wgs_source_dir(tmp_path, monkeypatch, built):
    src, rsrc = ac._sources(None, ac.WGS_TABLE, ac.WGS_RATIOS, ac.WGS_ENV_DIR)
    (tmp_path / "alpaca_corrected.parquet").write_bytes(src.read_bytes())
    (tmp_path / "wgd_ratios_corrected.tsv").write_bytes(rsrc.read_bytes())
    monkeypatch.setenv(ac.WGS_ENV_DIR, str(tmp_path))
    assert ac.build_wgs().height == built[1].height
