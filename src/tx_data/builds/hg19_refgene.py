"""Build hg19_refgene → canonical Parquet.

Ingests UCSC ``refGene.txt.gz`` (81 k transcript rows) + ``cytoBand.txt.gz``
(862 cytoband intervals), aggregates to one row per (gene × chromosome)
with min-max interval, and annotates every gene with the semicolon-joined
list of cytobands its interval overlaps.
"""

from __future__ import annotations

import polars as pl

from tx_data.builds._base import canonical_output_path, log
from tx_data.sources import resolve_source

TABLE = "hg19_refgene"

# UCSC refGene column order (no header row in the source file).
_REFGENE_COLS = (
    "bin",
    "name",         # transcript ID (NR_/NM_)
    "chrom",        # 'chrN'
    "strand",
    "txStart",
    "txEnd",
    "cdsStart",
    "cdsEnd",
    "exonCount",
    "exonStarts",
    "exonEnds",
    "score",
    "name2",        # gene symbol
    "cdsStartStat",
    "cdsEndStat",
    "exonFrames",
)

# UCSC cytoBand column order (no header row).
_CYTOBAND_COLS = ("chrom", "chromStart", "chromEnd", "band_name", "gie_stain")

# Chromosomes we keep. Everything else (chr*_random, chrUn_*, alt haplotype
# chromosomes) is dropped — paper analyses only run on the primary
# assembly.
_KEEP_CHR = {f"chr{c}" for c in list(range(1, 23)) + ["X", "Y", "M"]}


def _strip_chr(col: str = "chromosome") -> pl.Expr:
    return pl.col(col).str.replace(r"^chr", "").str.replace(r"^M$", "MT")


def build() -> pl.DataFrame:
    refgene_src = resolve_source(TABLE + "__refgene")
    cytoband_src = resolve_source(TABLE + "__cytoband")

    log(f"read {refgene_src}")
    refgene = pl.read_csv(
        refgene_src,
        separator="\t",
        has_header=False,
        new_columns=list(_REFGENE_COLS),
        null_values=["NA", ""],
        infer_schema_length=10000,
    )
    refgene = refgene.filter(pl.col("chrom").is_in(_KEEP_CHR))

    log(f"read {cytoband_src}")
    cytoband = pl.read_csv(
        cytoband_src,
        separator="\t",
        has_header=False,
        new_columns=list(_CYTOBAND_COLS),
    ).filter(pl.col("chrom").is_in(_KEEP_CHR))

    # Aggregate transcripts to one row per (gene_name, chrom).
    genes = (
        refgene.rename({"name2": "gene_name", "chrom": "chromosome"})
        .group_by(["gene_name", "chromosome"])
        .agg(
            pl.col("txStart").min().alias("start"),
            pl.col("txEnd").max().alias("end"),
            pl.col("strand").first().alias("strand"),
        )
    )

    # Cytoband overlap: for each gene, find every cytoBand row whose
    # interval intersects the gene's [start, end) on the same chromosome,
    # then semicolon-join their band names in sorted order.
    cytoband = cytoband.select(
        pl.col("chrom").alias("cyto_chrom"),
        pl.col("chromStart").alias("cyto_start"),
        pl.col("chromEnd").alias("cyto_end"),
        pl.col("band_name"),
    )
    # Cross-join per chromosome (small — ~26k genes × ~20 bands per chrom).
    joined = genes.join(
        cytoband, left_on="chromosome", right_on="cyto_chrom", how="left"
    ).filter(
        (pl.col("cyto_start") < pl.col("end"))
        & (pl.col("cyto_end") > pl.col("start"))
    )
    cyto_by_gene = joined.group_by(["gene_name", "chromosome"]).agg(
        pl.col("band_name")
        .unique()
        .sort()
        .str.concat(delimiter=";")
        .alias("cytoband_name")
    )

    genes = (
        genes.join(cyto_by_gene, on=["gene_name", "chromosome"], how="left")
        .with_columns(pl.col("cytoband_name").fill_null(""))
    )

    # Strip 'chr' prefix on chromosome to match tx_data canonical form
    # (matches muttable / alpaca / driver_list).
    genes = genes.with_columns(_strip_chr("chromosome"))

    genes = genes.select(
        ["gene_name", "chromosome", "start", "end", "strand", "cytoband_name"]
    ).sort(["chromosome", "start", "gene_name"])

    out = canonical_output_path(TABLE)
    genes.write_parquet(out)
    log(f"wrote {out} ({genes.height} rows × {genes.width} cols)")
    return genes


if __name__ == "__main__":
    build()
