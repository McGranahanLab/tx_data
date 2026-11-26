"""Build wgd_calls (WES) → canonical Parquet.

Since 2026-11-26 (Piotr) the WES calls come only from darwins_analysis `WGD_pipeline` run with
`config.wes.toml`: step 05 `wgd_calls.tsv`, read by the same code as `wgs_wgd_calls` (same columns). The
previous release table (`WGD/release_tx842/WGD_calls.tsv`) is no longer a tx_data source. Where the
pipeline output lives is a build parameter: `build(source_dir=...)`,
`scripts/build_wgd_calls.py --source-dir <darwins_analysis>/WGD_pipeline/output/wes/wgd/05_table`, or
TX_WES_WGD_CALLS_DIR (used by build_all); without one, the mock layout in nemo_mock/.

Derived columns (`add_derived`, also used by wgs_wgd_calls):
  - `is_root`: `parent == 'diploid'`
  - `is_wgd_event_clone`: `number_of_gds_at_clone > 0`
  - `cumulative_gds_to_clone`: sum along parent chain to MRCA (inclusive)
  - `has_wgd_ancestor`: `cumulative_gds_to_clone > 0`

Only `status == 'resolved'` rows are considered analytically usable — but we
keep unresolved rows in the canonical Parquet (with all derived columns still
computed) so downstream can filter as it sees fit.
"""

from __future__ import annotations

import polars as pl


TABLE = "wgd_calls"


def _cumulative_gds_per_tumour(pdf: pl.DataFrame) -> pl.DataFrame:
    """Compute cumulative-GDs-to-clone within a single tumour via tree walk."""
    # Map clone → (parent, ngds).
    rows = pdf.select(["clone", "parent", "number_of_gds_at_clone"]).to_dicts()
    parent_of = {r["clone"]: r["parent"] for r in rows}
    ngds_of = {r["clone"]: int(r["number_of_gds_at_clone"] or 0) for r in rows}

    cache: dict[str, int] = {}

    def cum(clone: str) -> int:
        if clone in cache:
            return cache[clone]
        if clone == "diploid":
            cache[clone] = 0
            return 0
        parent = parent_of.get(clone)
        if parent is None:
            # Malformed tree: clone with parent not in table.
            cache[clone] = ngds_of.get(clone, 0)
            return cache[clone]
        cache[clone] = ngds_of.get(clone, 0) + cum(parent)
        return cache[clone]

    cum_map = {c: cum(c) for c in parent_of.keys()}
    return pdf.with_columns(
        pl.col("clone")
        .replace_strict(cum_map, return_dtype=pl.Int64)
        .alias("cumulative_gds_to_clone")
    )


def add_derived(df: pl.DataFrame) -> pl.DataFrame:
    """is_root, is_wgd_event_clone, cumulative_gds_to_clone, has_wgd_ancestor (shared with wgs_wgd_calls)."""
    # Base derived columns (per-row, no tree walk).
    df = df.with_columns(
        [
            (pl.col("parent") == "diploid").alias("is_root"),
            (pl.col("number_of_gds_at_clone") > 0).alias("is_wgd_event_clone"),
        ]
    )

    # Tree-walked derived column: cumulative_gds_to_clone (per tumour).
    parts = []
    for tumour_id, sub in df.group_by("tumour_id"):
        parts.append(_cumulative_gds_per_tumour(sub))
    df = pl.concat(parts, how="vertical")

    df = df.with_columns(
        (pl.col("cumulative_gds_to_clone") > 0).alias("has_wgd_ancestor")
    )
    return df


ENV_DIR = "TX_WES_WGD_CALLS_DIR"


def build(source_dir=None) -> pl.DataFrame:
    """WES calls = step 05 wgd_calls.tsv of WGD_pipeline (config.wes.toml); shared reader with wgs_wgd_calls."""
    from tx_data.builds.wgs_wgd_calls import _source, build_table   # lazy: wgs_wgd_calls imports add_derived

    return build_table(TABLE, _source(source_dir, TABLE, ENV_DIR, "wes"))


if __name__ == "__main__":
    build()
