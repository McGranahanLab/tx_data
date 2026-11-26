"""Build every canonical Parquet from source, in dependency order.

Writes a `data/_build_manifest.json` at the end recording the nemo_root
and WES data_root used, the TX_NEMO_ROOT env var seen, whether this was a
mock build, the timestamp, per-table cohort (wes / wgs / reference) and
per-table row counts. Downstream consumers (e.g.
mut_essential_wgd/run_analysis.py) read this to prove which cohort a
given set of parquets came from.

The build is resilient: a builder that raises is logged, recorded in the
manifest under ``failed_tables``, and the pipeline continues with the
remaining builders. A human-readable summary is also written to
``data/_build_report.txt`` so missing inputs are easy to spot.
"""
from __future__ import annotations

import argparse
import json
import os
import traceback
from datetime import datetime, timezone
from pathlib import Path

import polars as pl

from tx_data.paths import build_root_for_mode

TABLES = (
    "muttable",
    "clinical",
    "wgd_calls",              # WGD_pipeline WES output (run the pipeline first); TX_WES_WGD_CALLS_DIR (else the mock)
    "alphamissense",
    "driver_list",
    "alpaca",
    "alpaca_wgd_ratios",
    # alpaca_corrected emits TWO parquets (alpaca_corrected, alpaca_wgd_ratios_corrected)
    "alpaca_corrected",
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
    # WES CONIPHER outputs from the tx842 release tree + WES segments from the ALPACA input (WGD_pipeline,
    # config.wes.toml);
    # after clone_proportions, against which conipher_mutations checks its samples
    "conipher_mutations",
    "conipher_tree_edges",
    "conipher_cluster_info",
    "msai_segments",
    "purity_ploidy",
    # ---- WGS cohort ----
    # Inventory first: the other WGS builds discover tumours and samples from it.
    "wgs_battenberg_inventory",
    "wgs_consensus_segments",
    "wgs_clone_proportions",
    "wgs_tree_clones",
    # CONIPHER (2026_08_mean_max_absent run) and ALPACA output for the WGS WGD pipeline
    "wgs_conipher_mutations",
    "wgs_conipher_tree_edges",
    "wgs_conipher_cluster_info",
    "wgs_alpaca",
    "wgs_alpaca_wgd_ratios",
    # trees (ALPACA input), cohort facts, chromosome lengths (2026-11-19)
    "alpaca_trees",
    "wgs_alpaca_trees",
    # single-sample / single-clone tumours per cohort (Piotr, 2026-11-25); after the cp_table and CONIPHER
    # tables they count from, before cohort_metadata, which summarises them
    "sample_clone_counts",
    "wgs_sample_clone_counts",
    "cohort_metadata",
    # chr_table emits TWO parquets (chr_table_hg19, chr_table_hg38)
    "chr_table",
    # WGD_pipeline WGS output (run the pipeline first); location from TX_WGS_WGD_CALLS_DIR (else the mock)
    "wgs_wgd_calls",
    # WGS ALPACA correction (WGD_pipeline WGS run, correct step); TX_WGS_ALPACA_CORRECTED_DIR (else the mock)
    "wgs_alpaca_corrected",
)

# Extra parquet names a builder may emit in addition to `<table>.parquet`.
# Used both for row-count discovery and for reporting "expected but missing"
# outputs when a build silently produces a subset.
EXTRA_OUTPUTS: dict[str, tuple[str, ...]] = {
    "clinical": ("clinical_per_tumour",),
    "kallisto": ("kallisto_paired_normals",),
    "conipher_trees": ("conipher_ccf", "conipher_clonality"),
    "chr_table": ("chr_table_hg19", "chr_table_hg38"),
    "wgs_alpaca_corrected": ("wgs_alpaca_wgd_ratios_corrected",),
}
# Builders with no parquet named after them: their primary output.
PRIMARY_OUTPUT = {"conipher_trees": "conipher_ccf", "chr_table": "chr_table_hg19"}


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


def _primary_output(name: str) -> str:
    """Primary parquet basename — same as table name except for conipher_trees."""
    # conipher_trees / chr_table have no file named after the builder; use their primary output.
    return PRIMARY_OUTPUT.get(name, name)


def _expected_outputs(name: str) -> tuple[str, ...]:
    """All parquet basenames a builder is expected to produce."""
    if name in PRIMARY_OUTPUT:
        return EXTRA_OUTPUTS[name]
    return (name,) + EXTRA_OUTPUTS.get(name, ())


def _collect_row_counts(
    data_dir: Path, failures: dict[str, dict[str, str]]
) -> dict[str, int]:
    """Row counts for every parquet produced *this run*.

    Skips any table listed in ``failures`` so that a stale parquet from a
    previous successful run is not reported as if freshly rebuilt.
    """
    row_counts: dict[str, int] = {}
    for name in TABLES:
        if name in failures:
            continue
        for sub in _expected_outputs(name):
            p = data_dir / f"{sub}.parquet"
            if p.is_file():
                row_counts[sub] = pl.read_parquet(p).height
    return row_counts


def _missing_outputs(data_dir: Path, name: str) -> list[str]:
    return [
        sub for sub in _expected_outputs(name)
        if not (data_dir / f"{sub}.parquet").is_file()
    ]


def _combined_issues(data_dir: Path, failures: dict[str, dict[str, str]]) -> pl.DataFrame:
    """Per-tumour issues (skipped / partial / warning) of the tables built this run -> data/build_issues.tsv."""
    from tx_data.builds._base import ISSUE_COLS, ISSUES_DIR

    frames = [pl.read_csv(ISSUES_DIR / f"{n}.tsv", separator="\t", infer_schema_length=0)
              for n in TABLES if n not in failures and (ISSUES_DIR / f"{n}.tsv").is_file()]
    df = (pl.concat(frames, how="diagonal_relaxed") if frames
          else pl.DataFrame(schema={c: pl.String for c in ISSUE_COLS})).fill_null("")
    df.write_csv(data_dir / "build_issues.tsv", separator="\t")
    return df


def _write_report(
    data_dir: Path,
    manifest: dict,
    failures: dict[str, dict[str, str]],
    log_fn,
) -> None:
    lines = []
    lines.append("tx_data build report")
    lines.append(f"built_at_utc: {manifest['built_at_utc']}")
    lines.append(f"nemo_root:    {manifest['nemo_root']}")
    lines.append(f"data_root:    {manifest['data_root']}")
    lines.append(f"is_mock:      {manifest['data_root_is_mock']}")
    lines.append("")

    lines.append("Per-table status:")
    row_counts = manifest["table_row_counts"]
    for name in TABLES:
        missing = _missing_outputs(data_dir, name)
        if name in failures:
            f = failures[name]
            lines.append(
                f"  FAIL    {name:<28} {f['error_type']}: {f['message']}"
            )
        elif missing:
            lines.append(
                f"  MISSING {name:<28} no output parquet(s): {', '.join(missing)}"
            )
        else:
            rows = row_counts.get(_primary_output(name), 0)
            extras = [
                f"{sub}={row_counts.get(sub, 0)}"
                for sub in EXTRA_OUTPUTS.get(name, ())
                if sub in row_counts
            ]
            extra_str = f"  [{', '.join(extras)}]" if extras else ""
            lines.append(f"  OK      {name:<28} rows={rows}{extra_str}")

    if failures:
        lines.append("")
        lines.append("Failures (with tracebacks):")
        for name, f in failures.items():
            lines.append(f"  - {name}: {f['error_type']}: {f['message']}")
            for tb_line in f["traceback"].rstrip().splitlines():
                lines.append(f"      {tb_line}")

    issues = _combined_issues(data_dir, failures)
    lines.append("")
    if issues.height:
        lines.append(f"Tumour issues ({issues.height}; full list in {data_dir / 'build_issues.tsv'}):")
        for (table, sev), g in issues.group_by("table", "severity", maintain_order=True):
            lines.append(f"  {table:<28} {sev:<8} {g.height:>4}: " + ", ".join(g["tumour_id"].head(10).to_list())
                         + (" ..." if g.height > 10 else ""))
        lines.append("  By reason:")
        for (reason,), g in issues.group_by("reason", maintain_order=True):
            lines.append(f"    {g.height:>4}  {reason[:140]}")
    else:
        lines.append("Tumour issues: none")

    report_path = data_dir / "_build_report.txt"
    report_path.write_text("\n".join(lines) + "\n")
    log_fn(f"wrote {report_path}")


def _write_manifest(
    data_dir: Path,
    data_root_fn,
    data_root_is_mock_fn,
    failures: dict[str, dict[str, str]],
    log_fn,
) -> dict:
    row_counts = _collect_row_counts(data_dir, failures)

    from tx_data.catalog import load_catalog
    from tx_data.paths import nemo_root

    manifest = {
        "built_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "nemo_root": str(nemo_root()),
        "data_root": str(data_root_fn()),
        "data_root_is_mock": data_root_is_mock_fn(),
        "env_TX_NEMO_ROOT": os.environ.get("TX_NEMO_ROOT"),
        "table_cohorts": {name: load_catalog(name if name != "chr_table" else PRIMARY_OUTPUT[name])["cohort"]
                          for name in TABLES},
        "table_row_counts": row_counts,
        "failed_tables": {
            name: {"error_type": f["error_type"], "message": f["message"]}
            for name, f in failures.items()
        },
    }
    manifest_path = data_dir / "_build_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    log_fn(f"wrote {manifest_path}")
    log_fn(
        f"nemo_root={manifest['nemo_root']} "
        f"is_mock={manifest['data_root_is_mock']} "
        f"tables={row_counts}"
    )
    return manifest


if __name__ == "__main__":
    args = _parse_args()
    _apply_build_mode(args.mode)

    from tx_data.builds import BUILDERS
    from tx_data.builds._base import DATA_DIR, log
    from tx_data.paths import data_root, data_root_is_mock

    failures: dict[str, dict[str, str]] = {}
    for name in TABLES:
        try:
            BUILDERS[name]()
        except Exception as e:
            tb = traceback.format_exc()
            failures[name] = {
                "error_type": type(e).__name__,
                "message": str(e),
                "traceback": tb,
            }
            log(f"FAILED {name}: {type(e).__name__}: {e}")

    manifest = _write_manifest(
        DATA_DIR, data_root, data_root_is_mock, failures, log,
    )
    _write_report(DATA_DIR, manifest, failures, log)

    if failures:
        log(
            f"build finished with {len(failures)} failure(s): "
            f"{', '.join(sorted(failures))}. "
            f"See {DATA_DIR / '_build_report.txt'} for details."
        )
    else:
        log("build finished with all tables succeeded.")
