"""Entry point: build sample_clone_counts (WES) and wgs_sample_clone_counts (WGS), then cohort_metadata.

Needs clone_proportions / wgs_clone_proportions (cp_table) built, and for tumours without a cp_table the
CONIPHER tables (conipher_mutations, conipher_tree_edges, wgs_conipher_*). Prints the single-sample /
single-clone tumours of each cohort.
"""
import polars as pl

from tx_data.builds import cohort_metadata
from tx_data.builds.sample_clone_counts import build, build_wgs

if __name__ == "__main__":
    for name, df in (("WES", build()), ("WGS", build_wgs())):
        odd = df.filter(pl.col("category") != "multi_sample_multi_clone")
        print(f"\n{name}: {df.height} tumours, {odd.height} single-sample and/or single-clone")
        with pl.Config(tbl_rows=-1):
            print(odd.select("tumour_id", "n_samples", "n_clones", "category", "in_cp_table", "source"))
    cohort_metadata.build()
