"""Ad-hoc: add an hg38 `mutation_id_hg38` column to muttable.

Muttable coordinates are hg19 (see catalog/muttable.yml). This script lifts
the `pos` of every mutation over to hg38 with the UCSC hg19→hg38 chain and
rebuilds the mutation_id in the hg38 build:

    {patient}_{tumour_ordinal}:{chr_hg38}:{pos_hg38}:{ref}:{alt}

Rows whose position cannot be lifted (deleted / ambiguous / non-unique)
get a null `mutation_id_hg38`.

NOTE: This is a coordinate liftover only. For indels the local reference
sequence in hg38 may differ; the coordinate is correct but the (ref, alt)
alleles are not re-verified against the hg38 sequence. Fine for joining
against hg38-coordinate tables; not sufficient to redo variant calling.

Requires `pyliftover` (installed into the tx_data venv). On first run,
pyliftover downloads the ~200 MB chain to ~/.pyliftover/. Pass
`--chain PATH/hg19ToHg38.over.chain.gz` to skip the download.

Usage:
    python scripts/liftover_muttable.py                    # data/muttable.parquet → data/muttable_hg38_lifted.parquet
    python scripts/liftover_muttable.py --chain /path/to/hg19ToHg38.over.chain.gz
"""

from __future__ import annotations

import argparse
from pathlib import Path

import polars as pl
from pyliftover import LiftOver

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_IN = REPO_ROOT / "data" / "muttable.parquet"
DEFAULT_OUT = REPO_ROOT / "data" / "muttable_hg38_lifted.parquet"


def _lift_one(lo: LiftOver, chrom: str | None, pos: int | None) -> tuple[str | None, int | None]:
    """Lift a single (chr, 1-based pos) to hg38. Returns (None, None) on failure."""
    if chrom is None or pos is None:
        return None, None
    ucsc_chrom = chrom if chrom.startswith("chr") else f"chr{chrom}"
    # pyliftover is 0-based; muttable pos is 1-based.
    hits = lo.convert_coordinate(ucsc_chrom, pos - 1)
    if not hits:
        return None, None
    new_c, new_p, _strand, _score = hits[0]
    clean_c = new_c[3:] if new_c.startswith("chr") else new_c
    return clean_c, new_p + 1


def build_lookup(df: pl.DataFrame, lo: LiftOver) -> pl.DataFrame:
    """Build a (chr, pos) → (chr_hg38, pos_hg38) lookup for the muttable's unique loci."""
    uniq = df.select("chr", "pos").unique()
    chrs = uniq["chr"].to_list()
    poss = uniq["pos"].to_list()

    new_chr: list[str | None] = []
    new_pos: list[int | None] = []
    for c, p in zip(chrs, poss, strict=True):
        nc, np_ = _lift_one(lo, c, p)
        new_chr.append(nc)
        new_pos.append(np_)

    return pl.DataFrame(
        {
            "chr": chrs,
            "pos": poss,
            "chr_hg38": new_chr,
            "pos_hg38": new_pos,
        }
    )


def add_mutation_id_hg38(df: pl.DataFrame, lookup: pl.DataFrame) -> pl.DataFrame:
    """Join the lookup, rebuild mutation_id in hg38, drop intermediate columns."""
    df = df.join(lookup, on=["chr", "pos"], how="left")

    # mutation_id prefix = "{patient}_{ordinal}" (everything before the first ':')
    prefix = pl.col("mutation_id").str.split(":").list.get(0)

    mut_id_hg38 = (
        pl.when(pl.col("chr_hg38").is_null() | pl.col("pos_hg38").is_null())
        .then(None)
        .otherwise(
            prefix
            + pl.lit(":")
            + pl.col("chr_hg38")
            + pl.lit(":")
            + pl.col("pos_hg38").cast(pl.String)
            + pl.lit(":")
            + pl.col("ref")
            + pl.lit(":")
            + pl.col("alt")
        )
        .alias("mutation_id_hg38")
    )

    return df.with_columns(mut_id_hg38).drop("chr_hg38", "pos_hg38")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", type=Path, default=DEFAULT_IN, help="Input muttable parquet (hg19).")
    ap.add_argument("--output", type=Path, default=DEFAULT_OUT, help="Output parquet.")
    ap.add_argument(
        "--chain",
        type=Path,
        default=None,
        help="Path to hg19ToHg38.over.chain.gz. Omit to let pyliftover download it.",
    )
    args = ap.parse_args()

    print(f"[liftover] reading {args.input}")
    df = pl.read_parquet(args.input)

    print("[liftover] loading chain")
    lo = LiftOver(str(args.chain)) if args.chain else LiftOver("hg19", "hg38")

    n_unique = df.select("chr", "pos").unique().height
    print(f"[liftover] lifting {n_unique} unique (chr, pos) loci")
    lookup = build_lookup(df, lo)

    n_lifted = lookup.filter(pl.col("pos_hg38").is_not_null()).height
    pct_loci = 100.0 * n_lifted / n_unique if n_unique else 0.0
    print(f"[liftover]   {n_lifted}/{n_unique} loci mapped ({pct_loci:.1f}%)")

    df_out = add_mutation_id_hg38(df, lookup)

    n_rows_mapped = df_out.filter(pl.col("mutation_id_hg38").is_not_null()).height
    pct_rows = 100.0 * n_rows_mapped / df_out.height if df_out.height else 0.0
    print(f"[liftover]   {n_rows_mapped}/{df_out.height} rows carry a mutation_id_hg38 ({pct_rows:.1f}%)")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    df_out.write_parquet(args.output)
    print(f"[liftover] wrote {args.output} ({df_out.height} rows × {df_out.width} cols)")


if __name__ == "__main__":
    main()
