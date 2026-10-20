"""Synthetic CONIPHER cases: determinism + simulation invariants."""

from __future__ import annotations

import csv
import hashlib
import importlib.util
import random
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "scripts" / "generate_conipher_synthetic.py"


@pytest.fixture(scope="module")
def gen():
    spec = importlib.util.spec_from_file_location("generate_conipher_synthetic", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module  # dataclasses resolve their module via sys.modules
    spec.loader.exec_module(module)
    return module


def _digest(case_dir: Path) -> dict[str, str]:
    return {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(case_dir.iterdir())}


def _read(path: Path) -> list[dict]:
    with path.open() as f:
        return list(csv.DictReader(f, delimiter="\t"))


def test_regeneration_is_byte_identical(gen, tmp_path):
    for case in gen.CASES:
        first = _digest(gen.write_case(case, tmp_path / "a"))
        second = _digest(gen.write_case(case, tmp_path / "b"))
        assert first == second, case.case_id


def test_committed_files_match_generator(gen, tmp_path):
    for case in gen.CASES:
        fresh = _digest(gen.write_case(case, tmp_path))
        committed = _digest(gen.OUT_ROOT / case.case_id)
        assert fresh == committed, f"{case.case_id}: re-run scripts/generate_conipher_synthetic.py"


def test_input_table_invariants(gen, tmp_path):
    for case in gen.CASES:
        case_dir = gen.write_case(case, tmp_path)
        rows = _read(case_dir / "input.tsv")
        truth = _read(case_dir / "truth_mutations.tsv")
        assert list(rows[0]) == gen.INPUT_COLUMNS
        assert len(rows) == len(truth) * case.n_regions
        # Every chromosome 1..22 is present (CONIPHER plotting needs no gaps).
        assert {r["CHR"] for r in rows} == set(gen.ALL_CHROMS)
        for r in rows:
            assert int(r["REF_COUNT"]) + int(r["VAR_COUNT"]) == int(r["DEPTH"])
        # A clone absent from a region contributes no variant reads there.
        absent = {
            (c.name, f"{case.case_id}_R{i + 1}")
            for c in case.clones
            for i, ccf in enumerate(c.ccf)
            if ccf == 0
        }
        clone_of = {(t["CHR"], t["POS"]): t["CLONE"] for t in truth}
        for r in rows:
            if (clone_of[(r["CHR"], r["POS"])], r["SAMPLE"]) in absent:
                assert r["VAR_COUNT"] == "0"


def test_validate_rejects_sum_rule_violation(gen):
    bad = gen.Case(
        "BAD",
        "children exceed parent",
        purity=(0.5,),
        clones=(
            gen.Clone("clone1", "diploid", 5, (1.0,)),
            gen.Clone("clone2", "clone1", 5, (0.7,)),
            gen.Clone("clone3", "clone1", 5, (0.6,)),
        ),
    )
    with pytest.raises(AssertionError, match="exceed"):
        gen.simulate(bad, random.Random(0))
