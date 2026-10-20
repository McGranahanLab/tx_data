"""Build wgd_calls → canonical Parquet.

Adds derived columns:
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

from tx_data.builds._base import canonical_output_path, log
from tx_data.sources import resolve_source

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


def build() -> pl.DataFrame:
    src = resolve_source(TABLE)
    log(f"read {src}")

    df = pl.read_csv(
        src, separator="\t", null_values=["NA", "", "nan"], infer_schema_length=10000
    )

    # Ensure integer types where the catalog says so.
    df = df.with_columns(
        [
            pl.col("number_of_gds_at_clone").cast(pl.Int64),
            pl.col("total_gds_per_tumour").cast(pl.Int64),
        ]
    )

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

    out = canonical_output_path(TABLE)
    df.write_parquet(out)
    log(f"wrote {out} ({df.height} rows × {df.width} cols)")
    return df


if __name__ == "__main__":
    build()
