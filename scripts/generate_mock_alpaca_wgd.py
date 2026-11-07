"""Generate the ALPACA per-clone WGD ratio mock (`wgd_ratio_scores.csv`).

Real tool: ``alpaca wgd`` (McGranahanLab/ALPACA-model, tested at 9d3c2f0),
``alpaca/scripts.py::run_calculate_wgd`` → ``alpaca/analysis.py::calculate_wgd_ratios``.
It runs on ONE tumour at a time (single-tumour ALPACA output + that tumour's
``tree_paths.json``) and writes ``<output_directory>/wgd_ratio_scores.csv``.

Unlike a hand-written mock, this script imports and runs the real
``calculate_wgd_ratios`` (and ALPACA's own ``save_dataframe``) on the existing
WES mock:

  * CN: ``ALPACA/output/cohort_results/all_tumours_combined.csv``
    (from ``generate_mocks.py``), split per tumour, as ``alpaca wgd`` expects.
  * Tree: the mock tree in ``nemo_mock/_cohort_manifest.json`` (clone →
    parent, the same tree as wgd_calls / the CONIPHER RDS / cp_table),
    converted to ALPACA's ``tree_paths.json`` format (trunk → leaf paths).

Both intermediate inputs are written to a temp dir and discarded; only the
per-tumour outputs land in nemo_mock, at the assumed release path
``ALPACA/output/tumour_results/<tumour_id>/wgd_ratio_scores.csv``
(catalog/alpaca_wgd_ratios.yml; real path unconfirmed).

Requires pandas + scipy (ALPACA's deps) and an ALPACA-model checkout:

    ALPACA_SRC=/path/to/ALPACA-model pixi run -e alpaca mocks-alpaca-wgd

Re-run after ``pixi run mocks`` whenever the ALPACA CN mock or the tree changes.
Deterministic (no randomness).
"""

from __future__ import annotations

import csv
import json
import os
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
MANIFEST = REPO_ROOT / "nemo_mock" / "_cohort_manifest.json"
WES_ROOT = (
    REPO_ROOT / "nemo_mock" / "project" / "proj-tracerx-lung" / "tctProjects"
    / "lungTx" / "Tx842" / "release" / "tracerx_842"
)
ALPACA_COMBINED = WES_ROOT / "ALPACA" / "output" / "cohort_results" / "all_tumours_combined.csv"
OUT_ROOT = WES_ROOT / "ALPACA" / "output" / "tumour_results"
OUT_NAME = "wgd_ratio_scores.csv"

ALPACA_SRC = Path(os.environ.get("ALPACA_SRC", "/tmp/alpaca_model"))


def _import_alpaca():
    if not (ALPACA_SRC / "alpaca" / "analysis.py").is_file():
        sys.exit(f"ALPACA-model checkout not found at {ALPACA_SRC} (set ALPACA_SRC)")
    sys.path.insert(0, str(ALPACA_SRC))
    from alpaca.analysis import calculate_wgd_ratios
    from alpaca.utils import save_dataframe

    return calculate_wgd_ratios, save_dataframe


def tree_paths(clones: list[dict]) -> list[list[str]]:
    """Clone → parent edges → ALPACA tree_paths.json (one trunk → leaf path per leaf)."""
    parent = {c["name"]: c["parent"] for c in clones}
    roots = [c for c, p in parent.items() if p == "diploid"]
    if len(roots) != 1:
        raise ValueError(f"expected one MRCA, got {roots}")
    parents = set(parent.values())
    paths = []
    for leaf in (c["name"] for c in clones if c["name"] not in parents):
        path = [leaf]
        while parent[path[0]] != "diploid":
            path.insert(0, parent[path[0]])
        paths.append(path)
    return paths


def main() -> None:
    calculate_wgd_ratios, save_dataframe = _import_alpaca()
    manifest = json.loads(MANIFEST.read_text())

    with ALPACA_COMBINED.open() as f:
        reader = csv.reader(f)
        header = next(reader)
        rows_by_tumour: dict[str, list[list[str]]] = {}
        for row in reader:
            rows_by_tumour.setdefault(row[header.index("tumour_id")], []).append(row)

    n = 0
    with tempfile.TemporaryDirectory() as tmp:
        for t in manifest["tumours"]:
            tid = t["tumour_id_canonical"]
            if not t["has_tree"]:
                continue  # ALPACA is not run on tumours without a tree
            if tid not in rows_by_tumour:
                raise ValueError(f"{tid} has a tree but no ALPACA CN rows")
            cn_path = Path(tmp) / f"{tid}_alpaca_output.csv"
            with cn_path.open("w", newline="") as f:
                w = csv.writer(f)
                w.writerow(header)
                w.writerows(rows_by_tumour[tid])
            tree_path = Path(tmp) / f"{tid}_tree_paths.json"
            tree_path.write_text(json.dumps(tree_paths(t["clones"])))

            ratios = calculate_wgd_ratios(results_path=str(cn_path), tree_path=str(tree_path))
            save_dataframe(ratios, str(OUT_ROOT / tid), OUT_NAME)
            n += 1
            print(f"wrote {OUT_ROOT / tid / OUT_NAME} ({len(ratios)} clones)")
    print(f"done: {n} tumours")


if __name__ == "__main__":
    main()
