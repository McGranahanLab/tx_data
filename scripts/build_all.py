"""Build every canonical Parquet from source, in dependency order.

Writes a `data/_build_manifest.json` at the end recording the nemo_root
and WES data_root used, the TX_NEMO_ROOT env var seen, whether this was a
mock build, the timestamp, per-table cohort (wes / wgs / reference) and
per-table row counts. Downstream consumers (e.g.
mut_essential_wgd/run_analysis.py) read this to prove which cohort a
given set of parquets came from.
"""
from __future__ import annotations

import argparse
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import polars as pl

from tx_data.paths import build_root_for_mode

TABLES = (
    "muttable",
    "clinical",
    "wgd_calls",
    "alphamissense",
    "driver_list",
    "alpaca",
    "kallisto",
    "clone_proportions",
    # id_crosswalk must come before clone_info / seeding_table — the latter
    # two read it (via tx_data.ids) to translate CRUK ↔ LTX patient
    # identifiers when the source and target conventions differ.
    "id_crosswalk",
    "clone_info",
    "seeding_table",
    "hg19_refgene",
    # conipher_trees emits TWO parquets (conipher_ccf, conipher_clonality);
    # _write_manifest discovers both via its special-case fallbacks below.
    "conipher_trees",
    "purity_ploidy",
)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "mode",
        nargs="?",
        choices=("sandbox", "prod"),
        help="Override the source root used for the build.",
    )
    return parser.parse_args()


def _apply_build_mode(mode: str | None) -> None:
    root = build_root_for_mode(mode)
    if root is None:
        return
    os.environ["TX_NEMO_ROOT"] = str(root)


def _write_manifest(
    data_dir: Path,
    data_root_fn,
    data_root_is_mock_fn,
    log_fn,
) -> None:
    row_counts: dict[str, int] = {}
    for name in TABLES:
        p = data_dir / f"{name}.parquet"
        if p.is_file():
            row_counts[name] = pl.read_parquet(p).height
        # clinical builder additionally emits clinical_per_tumour; capture it.
        if name == "clinical":
            p2 = data_dir / "clinical_per_tumour.parquet"
            if p2.is_file():
                row_counts["clinical_per_tumour"] = pl.read_parquet(p2).height
        # kallisto builder additionally emits kallisto_paired_normals.
        if name == "kallisto":
            p2 = data_dir / "kallisto_paired_normals.parquet"
            if p2.is_file():
                row_counts["kallisto_paired_normals"] = pl.read_parquet(p2).height
        # conipher_trees builder emits two parquets, neither named after the table.
        if name == "conipher_trees":
            for sub in ("conipher_ccf", "conipher_clonality"):
                p2 = data_dir / f"{sub}.parquet"
                if p2.is_file():
                    row_counts[sub] = pl.read_parquet(p2).height

    from tx_data.catalog import load_catalog
    from tx_data.paths import nemo_root

    manifest = {
        "built_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "nemo_root": str(nemo_root()),
        "data_root": str(data_root_fn()),
        "data_root_is_mock": data_root_is_mock_fn(),
        "env_TX_NEMO_ROOT": os.environ.get("TX_NEMO_ROOT"),
        "table_cohorts": {name: load_catalog(name)["cohort"] for name in TABLES},
        "table_row_counts": row_counts,
    }
    manifest_path = data_dir / "_build_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    log_fn(f"wrote {manifest_path}")
    log_fn(
        f"nemo_root={manifest['nemo_root']} "
        f"is_mock={manifest['data_root_is_mock']} "
        f"tables={row_counts}"
    )


if __name__ == "__main__":
    args = _parse_args()
    _apply_build_mode(args.mode)

    from tx_data.builds import BUILDERS
    from tx_data.builds._base import DATA_DIR, log
    from tx_data.paths import data_root, data_root_is_mock

    for name in TABLES:
        BUILDERS[name]()
    _write_manifest(DATA_DIR, data_root, data_root_is_mock, log)
