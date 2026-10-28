"""wgs_battenberg_inventory: build on the WGS mock + unit checks on parsing.

Self-contained: regenerates the WGS mock and builds this one table, so it
doesn't depend on test_pipeline.py having run first.
"""

from __future__ import annotations

import json

import polars as pl
import pytest

from tx_data.builds import wgs_battenberg_inventory as inv
from tx_data.paths import DEFAULT_MOCK_ROOT, nemo_path

MANIFEST = DEFAULT_MOCK_ROOT / "_cohort_manifest_wgs.json"


@pytest.fixture(scope="module")
def table(monkeypatch_module) -> pl.DataFrame:
    import runpy

    from tx_data.paths import REPO_ROOT

    runpy.run_path(str(REPO_ROOT / "scripts" / "generate_mocks_wgs.py"), run_name="__main__")
    return inv.build()


@pytest.fixture(scope="module")
def monkeypatch_module():
    mp = pytest.MonkeyPatch()
    mp.delenv("TX_NEMO_ROOT", raising=False)
    mp.delenv("TX_ID_CONVENTION", raising=False)
    yield mp
    mp.undo()


def test_columns_and_key(table):
    assert table.columns == inv.OUTPUT_COLUMNS
    assert table["sample_name_hash"].is_unique().all()


def test_matches_manifest(table):
    m = json.loads(MANIFEST.read_text())
    # tumour_id always comes from the tumour label; non-standard sample
    # names only lose `region`.
    expected = {
        s["sample_name_hash"]: (
            t["patient_id"],
            t["tumour_id_canonical"],
            t["wgs_tumour_label"],
            s["region"] if s["standard_name"] else None,
            s["purity"],
        )
        for t in m["tumours"] for s in t["samples"]
    }
    got = {
        r["sample_name_hash"]: (
            r["patient_id"], r["tumour_id"], r["wgs_tumour_label"], r["region"], r["purity"]
        )
        for r in table.iter_rows(named=True)
    }
    assert got == expected


def test_paths_normalised_and_resolvable(table):
    for col in ("results_dir", "msai_path", "mosdepth_summary_path", "subclones_path"):
        vals = table[col].drop_nulls()
        assert vals.str.starts_with("/nemo/").all(), col
        assert not vals.str.contains("//").any(), col
        assert not vals.str.ends_with("/").any(), col
    # Each path maps into nemo_mock/ in a sandbox build.
    p = nemo_path(table["msai_path"][0])
    assert DEFAULT_MOCK_ROOT in p.parents


def test_identity_wes_consistent(table):
    assert table["patient_id"].null_count() == 0
    assert table["tumour_id"].null_count() == 0
    assert table["tumour_id"].str.contains(r"^LTX\d+-Tumour\d+$").all()
    assert table["patient_id"].str.contains(r"^LTX\d+$").all()
    # Multi-tumour patient: suffixed labels; single-tumour: bare patient id.
    multi = table.filter(pl.col("patient_id") == "LTX0104").sort("tumour_id")
    assert multi["wgs_tumour_label"].to_list() == ["LTX0104_tumour1", "LTX0104_tumour2"]
    assert multi["tumour_id"].to_list() == ["LTX0104-Tumour1", "LTX0104-Tumour2"]
    single = table.filter(pl.col("patient_id") == "LTX0001")
    assert single["wgs_tumour_label"].to_list() == ["LTX0001"]
    assert single["tumour_id"].to_list() == ["LTX0001-Tumour1"]
    # Non-standard sample name keeps label-derived identity, loses region.
    ln = table.filter(pl.col("patient_id") == "LTX0106")
    assert ln["tumour_id"].to_list() == ["LTX0106-Tumour1"] and ln["region"].is_null().all()


def test_battenberg_dirs_use_tumour_label(table):
    for label, d in table.select("wgs_tumour_label", "results_dir").iter_rows():
        assert d.endswith(f"/output/{label}/results"), (label, d)


def test_no_rho_psi_columns(table):
    # purity / ploidy are the values to use; label tokens are not surfaced.
    assert not {"solution_psi", "solution_rho"} & set(table.columns)
    assert table["subclones_file"].str.contains(table["optimal_solution"], literal=True).all()


def test_wgs_hashes_disjoint_from_wes_hashes(table):
    from tx_data.sources import resolve_source

    wes = pl.read_csv(resolve_source("sample_hash_map"), separator="\t")
    wgs_suffixes = set(table["sample_name_hash"].str.split("--").list.last())
    assert not wgs_suffixes & set(wes["sample_hash"].cast(pl.String))


def test_wgs_tumour_label_parsing():
    import tx_data.builds  # noqa: F401  (tx_data.ids has a circular import if loaded first)
    from tx_data.ids import parse_wgs_tumour_label, wgs_tumour_id

    assert wgs_tumour_id("LTX0001") == "LTX0001-Tumour1"
    assert wgs_tumour_id("LTX0002_tumour1") == "LTX0002-Tumour1"
    assert wgs_tumour_id("LTX0002_tumour2") == "LTX0002-Tumour2"
    assert wgs_tumour_id("LTX0002_tumour_2") == "LTX0002-Tumour2"
    assert parse_wgs_tumour_label("CRUK0001") == ("CRUK0001", 1)
    for bad in ("LTX0001-Tumour1", "LTX0001_SU_T1", "patient1", ""):
        with pytest.raises(ValueError):
            parse_wgs_tumour_label(bad)


def test_parse_sample_rejects_unexpected_names():
    assert inv._parse_sample("LTX0001_SU_T1-R1--38da873ff7") == ("LTX0001", 1, "R1")
    for bad in ("LTX0001_SU_LN1--abc", "LTX0001-T1-R1--abc", "LTX0001_SU_T1-R1"):
        with pytest.raises(ValueError):
            inv._parse_sample(bad)

