"""Build sample_clone_counts (WES) and wgs_sample_clone_counts (WGS) → canonical Parquet.

Per tumour: how many samples and how many clones it has, and whether it is single-sample, single-clone or
both (Piotr, 2026-11-25). Such tumours can behave unusually: an analysis designed for multi-sample or
multi-clone tumours may fail on them. The two cohorts differ (different samples were sequenced), so each
has its own table; cohort_metadata summarises both.

Source: the ALPACA input clone composition `cp_table.csv`, i.e. the built `clone_proportions` /
`wgs_clone_proportions` (samples = its sample columns, clones = its rows). Tumours with CONIPHER output
but no cp_table (e.g. not run through ALPACA) are counted from CONIPHER instead: samples from the
`conipher_mutations` table, clones from the default tree in `conipher_tree_edges` (a `1 -> 1` self-loop is
one clone). `source` says which was used; `in_cp_table` is false for those tumours.
"""

from __future__ import annotations

import polars as pl

from tx_data.builds._base import canonical_output_path, log

SCHEMA = {"tumour_id": pl.String, "n_samples": pl.Int64, "n_clones": pl.Int64, "single_sample": pl.Boolean,
          "single_clone": pl.Boolean, "category": pl.String, "in_cp_table": pl.Boolean, "source": pl.String}
INPUTS = {  # table -> (clone_proportions, conipher_mutations, conipher_tree_edges)
    "sample_clone_counts": ("clone_proportions", "conipher_mutations", "conipher_tree_edges"),
    "wgs_sample_clone_counts": ("wgs_clone_proportions", "wgs_conipher_mutations", "wgs_conipher_tree_edges"),
}
CATEGORIES = {(True, True): "single_sample_single_clone", (True, False): "single_sample",
              (False, True): "single_clone", (False, False): "multi_sample_multi_clone"}


def _read(name: str, columns: list[str]) -> pl.DataFrame | None:
    p = canonical_output_path(name)
    if not p.is_file():
        log(f"warning: {name} not built; skipped")
        return None
    return pl.read_parquet(p, columns=columns)


def counts(table: str) -> pl.DataFrame:
    cp_name, mut_name, edge_name = INPUTS[table]
    parts = []
    cp = _read(cp_name, ["tumour_id", "clone", "sample_id"])
    if cp is not None:
        parts.append(cp.group_by("tumour_id").agg(pl.col("sample_id").n_unique().alias("n_samples"),
                                                  pl.col("clone").n_unique().alias("n_clones"))
                       .with_columns(pl.lit("cp_table").alias("source")))
    in_cp = set(parts[0]["tumour_id"]) if parts else set()
    mut, edges = _read(mut_name, ["tumour_id", "sample"]), _read(edge_name, ["tumour_id", "parent_cluster",
                                                                             "child_cluster"])
    if mut is not None:
        samples = mut.filter(~pl.col("tumour_id").is_in(list(in_cp))).group_by("tumour_id").agg(
            pl.col("sample").n_unique().alias("n_samples"))
        if edges is not None:
            clones = (pl.concat([edges.select("tumour_id", pl.col("parent_cluster").alias("c")),
                                 edges.select("tumour_id", pl.col("child_cluster").alias("c"))])
                        .group_by("tumour_id").agg(pl.col("c").n_unique().alias("n_clones")))
            samples = samples.join(clones, on="tumour_id", how="left")
        else:
            samples = samples.with_columns(pl.lit(None, pl.Int64).alias("n_clones"))
        parts.append(samples.with_columns(pl.lit("conipher").alias("source")))
    if not parts:
        return pl.DataFrame(schema=SCHEMA)
    df = pl.concat(parts, how="diagonal_relaxed").with_columns(
        (pl.col("n_samples") == 1).alias("single_sample"), (pl.col("n_clones") == 1).alias("single_clone"),
        pl.col("tumour_id").is_in(list(in_cp)).alias("in_cp_table"))
    df = df.with_columns(pl.struct("single_sample", "single_clone").map_elements(
        lambda r: CATEGORIES.get((r["single_sample"], r["single_clone"])), return_dtype=pl.String)
        .fill_null("clones_unknown").alias("category"))           # CONIPHER tumour without a tree
    return df.select([pl.col(c).cast(t) for c, t in SCHEMA.items()]).sort("tumour_id")


def _build(table: str) -> pl.DataFrame:
    df = counts(table)
    out = canonical_output_path(table)
    df.write_parquet(out)
    summary = ", ".join(f"{c}={n}" for c, n in df.group_by("category").len().sort("category").iter_rows())
    log(f"{table}: {df.height} tumours ({summary}); {df.filter(~pl.col('in_cp_table')).height} without cp_table")
    log(f"wrote {out} ({df.height} rows × {df.width} cols)")
    return df


def build() -> pl.DataFrame:
    return _build("sample_clone_counts")


def build_wgs() -> pl.DataFrame:
    return _build("wgs_sample_clone_counts")


if __name__ == "__main__":
    build()
    build_wgs()
