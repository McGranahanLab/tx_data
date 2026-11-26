"""Per-tumour failures don't stop a build: the tumour is left out and the reason recorded (Piotr, 2026-11-26)."""

from __future__ import annotations

import polars as pl
import pytest

from tx_data.builds._base import ISSUES_DIR


@pytest.fixture
def mock_root(monkeypatch):
    monkeypatch.delenv("TX_NEMO_ROOT", raising=False)


def _issues(table: str) -> pl.DataFrame:
    return pl.read_csv(ISSUES_DIR / f"{table}.tsv", separator="\t", infer_schema_length=0)


@pytest.mark.parametrize("module,table,bad", [
    ("conipher_mutations", "conipher_mutations", "LTX0002-Tumour1"),
    ("msai_segments", "msai_segments", "LTX0002-Tumour1"),
    ("conipher_cluster_info", "conipher_cluster_info", "LTX0002-Tumour1"),
])
def test_unreadable_tumour_is_skipped_and_recorded(mock_root, monkeypatch, module, table, bad):
    import importlib

    from tx_data.builds import BUILDERS

    BUILDERS["clone_proportions"]()
    mod = importlib.import_module(f"tx_data.builds.{module}")
    real = mod.read_one

    def flaky(*args):
        tid = args[0] if module == "msai_segments" else f"{args[0]}-Tumour{args[1]}"
        if tid == bad:
            raise ValueError("synthetic broken file")
        return real(*args)

    monkeypatch.setattr(mod, "read_one", flaky)
    df = BUILDERS[table]()
    assert bad not in set(df["tumour_id"]) and df["tumour_id"].n_unique() > 1
    iss = _issues(table).filter(pl.col("tumour_id") == bad)
    assert iss.row(0, named=True)["severity"] == "skipped" and "synthetic broken file" in iss["reason"][0]


def test_missing_file_is_recorded(mock_root, tmp_path, monkeypatch):
    """A tumour folder without the table's file is skipped with the path in the record."""
    from tx_data.builds import msai_segments as ms
    from tx_data.builds._base import Issues

    root = tmp_path / "ALPACA" / "input"
    (root / "LTX9999-Tumour1").mkdir(parents=True)
    monkeypatch.setattr(ms, "resolve_source", lambda *a, **k: root / k["tumour_id"] / "ALPACA_input_table.csv")
    issues = Issues("msai_segments_test")
    assert ms.discover(issues) == []
    assert issues.rows[0]["tumour_id"] == "LTX9999-Tumour1" and issues.rows[0]["severity"] == "skipped"
