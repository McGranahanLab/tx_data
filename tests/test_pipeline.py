"""End-to-end smoke test: mock → build → assert invariants on canonical Parquets."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import polars as pl
import pytest

REPO = Path(__file__).resolve().parents[1]
SRC = REPO / "src"
SCRIPTS = REPO / "scripts"
DATA = REPO / "data"


def _run(script: Path) -> None:
    r = subprocess.run(
        [sys.executable, str(script)],
        cwd=REPO,
        # Keep the caller's PATH: the conipher_trees build needs Rscript from the pixi env.
        env={"PYTHONPATH": str(SRC), "PATH": os.environ.get("PATH", "/usr/bin:/bin")},
        check=False,
        capture_output=True,
    )
    # check=False so a failing build reports its own stderr (e.g. Rscript not
    # on PATH when pytest runs outside `pixi run`) instead of a bare
    # CalledProcessError on whichever test runs first.
    assert r.returncode == 0, f"{script.name} exited {r.returncode}:\n{r.stderr.decode()[-3000:]}"


@pytest.fixture(scope="module", autouse=True)
def _build_all():
    _run(SCRIPTS / "generate_mocks.py")
    _run(SCRIPTS / "generate_mocks_wgs.py")
    _run(SCRIPTS / "build_all.py")


def test_muttable_shape():
    df = pl.read_parquet(DATA / "muttable.parquet")
    # Canonical column names shared with every other per-tumour table.
    assert "tumour_id" in df.columns
    assert "clone" in df.columns
    assert "alt" in df.columns
    assert "patient_tumour" not in df.columns
    assert "mutation_cluster" not in df.columns
    assert "var" not in df.columns
    # Canonical tumour_id form applied.
    assert df["tumour_id"].str.contains("-Tumour").all()
    # Clone form: string 'cloneN' or null.
    assert df["clone"].drop_nulls().str.starts_with("clone").all()
    # `chr` stripped of any `chr` prefix (values like "17", "X", "MT").
    assert not df["chr"].str.starts_with("chr").any()


def test_wgd_calls_derived_columns():
    df = pl.read_parquet(DATA / "wgd_calls.parquet")
    for col in ("is_root", "is_wgd_event_clone",
                "cumulative_gds_to_clone", "has_wgd_ancestor"):
        assert col in df.columns
    # Cumulative WGD at MRCA equals number_of_gds_at_clone at MRCA.
    roots = df.filter(pl.col("is_root"))
    assert (roots["cumulative_gds_to_clone"] == roots["number_of_gds_at_clone"]).all()
    # Cumulative WGD is non-decreasing along parent→child (spot-check on one tumour).
    for tumour_id, sub in df.group_by("tumour_id"):
        cmap = dict(zip(sub["clone"], sub["cumulative_gds_to_clone"]))
        pmap = dict(zip(sub["clone"], sub["parent"]))
        for c, p in pmap.items():
            if p == "diploid":
                continue
            assert cmap[c] >= cmap[p], f"{tumour_id[0]}: {c}({cmap[c]}) < {p}({cmap[p]})"


def test_wgd_calls_class_values():
    df = pl.read_parquet(DATA / "wgd_calls.parquet")
    resolved = df.filter(pl.col("status") == "resolved")
    seen_classes = set(resolved["class"].drop_nulls().unique().to_list())
    # Mock cohort covers all four class values.
    assert {"no_wgd", "clonal_wgd", "mut_supported", "ploidy_only"}.issubset(seen_classes)


def test_alphamissense_dedup():
    df = pl.read_parquet(DATA / "alphamissense.parquet")
    # One row per (tumour_id, chr, pos, alt).
    n = df.select(["tumour_id", "chr", "pos", "alt"]).n_unique()
    assert n == df.height


def test_clinical_per_tumour_unpivot():
    per_p = pl.read_parquet(DATA / "clinical.parquet")
    per_t = pl.read_parquet(DATA / "clinical_per_tumour.parquet")
    # Multi-tumour patient LTX0008 should have 2 rows in per_tumour.
    assert per_p.filter(pl.col("Patient_ID") == "LTX0008").height == 1
    assert per_t.filter(pl.col("Patient_ID") == "LTX0008").height == 2
    # tumour_id in canonical form.
    assert per_t["tumour_id"].str.contains("-Tumour").all()


def test_id_normalisation_consistency():
    """The tumour_id form in muttable must match wgd_calls after normalisation."""
    mut = pl.read_parquet(DATA / "muttable.parquet")
    wgd = pl.read_parquet(DATA / "wgd_calls.parquet")
    mut_tumours = set(mut["tumour_id"].unique().to_list())
    wgd_tumours = set(wgd["tumour_id"].unique().to_list())
    # muttable side already normalised by builder.
    common = mut_tumours & wgd_tumours
    assert common, f"no overlap; mut={list(mut_tumours)[:3]} wgd={list(wgd_tumours)[:3]}"


def test_chr_columns_have_no_prefix():
    """Every table with a chr column should carry values without `chr` prefix."""
    for path, cols in [
        (DATA / "muttable.parquet", ["chr"]),
        (DATA / "alphamissense.parquet", ["chr"]),
        (DATA / "alpaca.parquet", ["segment_chr"]),
        (DATA / "driver_list.parquet", ["chr_hg19", "chr_hg38"]),
    ]:
        df = pl.read_parquet(path)
        for col in cols:
            starts = df[col].drop_nulls().str.starts_with("chr")
            assert not starts.any(), f"{path.name}:{col} has `chr` prefix on some rows"


def test_clone_proportions_shape():
    """Long-form with the expected columns and sensible invariants."""
    df = pl.read_parquet(DATA / "clone_proportions.parquet")
    assert set(df.columns) == {"tumour_id", "clone", "sample_id", "proportion"}
    # Proportions are in [0, 1].
    assert df["proportion"].min() >= 0.0
    assert df["proportion"].max() <= 1.0 + 1e-9
    # Per (tumour, sample) they sum to ~1.
    sums = df.group_by(["tumour_id", "sample_id"]).agg(pl.col("proportion").sum().alias("s"))
    assert ((sums["s"] - 1.0).abs() < 1e-6).all()
    # tumour_id in canonical form; sample_id retains full hashed suffix.
    assert df["tumour_id"].str.contains("-Tumour").all()
    assert df["sample_id"].str.contains("--").all()


def test_clone_proportions_joins_wgd_calls():
    """Every (tumour_id, clone) in clone_proportions should exist in wgd_calls."""
    cp = pl.read_parquet(DATA / "clone_proportions.parquet")
    wgd = pl.read_parquet(DATA / "wgd_calls.parquet")
    cp_keys = cp.select(["tumour_id", "clone"]).unique()
    wgd_keys = wgd.select(["tumour_id", "clone"]).unique()
    orphans = cp_keys.join(wgd_keys, on=["tumour_id", "clone"], how="anti")
    assert orphans.height == 0, f"clones in cp_table not in wgd_calls: {orphans}"


def test_alpaca_wgd_ratios_schema_and_pk():
    """Columns/types follow the catalog; one row per (tumour_id, clone); one MRCA per tumour."""
    from tx_data.catalog import load_catalog

    df = pl.read_parquet(DATA / "alpaca_wgd_ratios.parquet")
    cat = load_catalog("alpaca_wgd_ratios")
    sql_to_pl = {"VARCHAR": pl.String, "DOUBLE": pl.Float64, "BIGINT": pl.Int64, "BOOLEAN": pl.Boolean}
    assert df.schema == {c["name"]: sql_to_pl[c["type"]] for c in cat["columns"]}
    assert df.height > 0
    assert df.select(cat["primary_key"]).is_unique().all()
    assert df["tumour_id"].str.contains("-Tumour").all()
    assert df["clone"].str.starts_with("clone").all()
    mrca = df.filter(pl.col("is_mrca"))
    assert (mrca["parent"] == "diploid").all()
    assert sorted(mrca["tumour_id"]) == sorted(df["tumour_id"].unique())
    # Every (tumour, clone) and every non-diploid parent is a clone in ALPACA.
    alpaca = pl.read_parquet(DATA / "alpaca.parquet").select("tumour_id", "clone").unique()
    assert df.join(alpaca, on=["tumour_id", "clone"], how="anti").height == 0
    parents = df.filter(~pl.col("is_mrca")).select("tumour_id", pl.col("parent").alias("clone"))
    assert parents.join(alpaca, on=["tumour_id", "clone"], how="anti").height == 0


def test_alpaca_wgd_ratios_reads_legacy_four_column_file(tmp_path):
    """Older ALPACA wrote only tumour_id, clone, parent, ratio_A: optional columns → null."""
    from tx_data.builds.alpaca_wgd_ratios import OPTIONAL, SCHEMA, read_wgd_ratio_file

    p = tmp_path / "wgd_ratio_scores.csv"
    p.write_text("tumour_id,clone,parent,ratio_A\n"
                 "LTX0001-Tumour1,clone1,diploid,2.0\nLTX0001-Tumour1,clone2,clone1,\n")
    df, missing = read_wgd_ratio_file(p, "LTX0001-Tumour1")
    assert missing == list(OPTIONAL)
    assert df.schema == SCHEMA
    assert df["ratio_A"].to_list() == [2.0, None]
    assert df.select(list(OPTIONAL)).null_count().row(0) == (2,) * len(OPTIONAL)
    with pytest.raises(ValueError, match="directory"):
        read_wgd_ratio_file(p, "LTX0002-Tumour1")
