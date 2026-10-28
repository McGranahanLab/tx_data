"""Build wgs_consensus_segments → canonical Parquet.

MSAI-adjusted fractional allele-specific copy number per WGS sample on
consensus segments, read from each tumour's Battenberg `tree_input/` folder:

  battenberg_wgs/output/<label>/tree_input/<label>_copynumber_MSAI_adjusted_consensus_segments_tree_input.tsv

Consensus segmentation: every sample of a tumour has the same breakpoints;
the build logs any tumour where that doesn't hold.

Tumours are discovered from the built `wgs_battenberg_inventory` (one file
per `wgs_tumour_label`), not by globbing, so a missing file is reported
instead of silently absent. Source columns COPY_NUMBER_MAJOR / _MINOR are
also accepted under their legacy names COPY_NUMBER_A / _B (always
major/minor, Piotr 2026-10-28).

Source `SAMPLE` is a short form of the sample name (`LTX0001_SU_T1.R1`);
it is matched to the inventory `sample_name_hash` within the tumour via
`tx_data.wgs_samples`. Values are kept exactly as on source: no rounding,
no clipping, no major/minor swapping.
"""

from __future__ import annotations

import polars as pl

from tx_data import wgs_samples
from tx_data.builds._base import canonical_output_path, log
from tx_data.sources import resolve_source

TABLE = "wgs_consensus_segments"

CN_COLUMN_VARIANTS = (
    ("COPY_NUMBER_MAJOR", "COPY_NUMBER_MINOR"),
    ("COPY_NUMBER_A", "COPY_NUMBER_B"),        # older files
)

SCHEMA = {
    "tumour_id": pl.String,
    "patient_id": pl.String,
    "wgs_tumour_label": pl.String,
    "sample_name_hash": pl.String,
    "sample": pl.String,
    "chr": pl.String,
    "start": pl.Int64,
    "end": pl.Int64,
    "cn_major": pl.Float64,
    "cn_minor": pl.Float64,
    "cn_total": pl.Float64,
}


def _read_one(label: str, path) -> pl.DataFrame:
    df = pl.read_csv(
        path,
        separator="\t",
        null_values=["NA", "", "nan"],
        infer_schema_length=10000,
        schema_overrides={"SAMPLE": pl.String, "CHR": pl.String},
    )
    for major_col, minor_col in CN_COLUMN_VARIANTS:
        if {major_col, minor_col} <= set(df.columns):
            break
    else:
        raise ValueError(
            f"{path}: no copy-number columns among {CN_COLUMN_VARIANTS}; got {df.columns!r}"
        )
    if major_col != CN_COLUMN_VARIANTS[0][0]:
        log(f"{label}: legacy header {major_col}/{minor_col} read as major/minor")
    missing = {"SAMPLE", "CHR", "STARTPOS", "ENDPOS"} - set(df.columns)
    if missing:
        raise ValueError(f"{path}: missing columns {sorted(missing)}; got {df.columns!r}")

    return df.select(
        pl.lit(label).alias("wgs_tumour_label"),
        pl.col("SAMPLE").alias("sample"),
        pl.col("CHR").str.replace(r"^chr", "").alias("chr"),
        pl.col("STARTPOS").cast(pl.Int64).alias("start"),
        pl.col("ENDPOS").cast(pl.Int64).alias("end"),
        pl.col(major_col).cast(pl.Float64).alias("cn_major"),
        pl.col(minor_col).cast(pl.Float64).alias("cn_minor"),
    )


def build() -> pl.DataFrame:
    wgs_samples.clear_cache()
    tumours = wgs_samples.tumours()

    frames: list[pl.DataFrame] = []
    absent: list[str] = []
    for label in tumours.get_column("wgs_tumour_label").to_list():
        path = resolve_source(TABLE, wgs_tumour_label=label)
        if not path.is_file():
            absent.append(label)
            continue
        frames.append(_read_one(label, path))
    log(f"read consensus segments for {len(frames)}/{tumours.height} WGS tumours")
    if absent:
        log(f"warning: no consensus segments file for {len(absent)} tumour(s): {absent}")

    if not frames:
        df = pl.DataFrame(schema=SCHEMA)
    else:
        df = pl.concat(frames, how="vertical")
        pairs = df.select("wgs_tumour_label", "sample").unique()
        matched = {
            (lb, s): wgs_samples.match_sample(lb, s) for lb, s in pairs.iter_rows()
        }
        unmatched = [f"{lb}:{s}" for (lb, s), v in matched.items() if v is None]
        if unmatched:
            log(
                f"warning: {len(unmatched)} SAMPLE value(s) not matched to the inventory; "
                f"sample_name_hash left NULL. First: {unmatched[:5]}"
            )
        lookup = pl.DataFrame(
            [(lb, s, v) for (lb, s), v in matched.items()],
            schema={"wgs_tumour_label": pl.String, "sample": pl.String,
                    "sample_name_hash": pl.String},
            orient="row",
        )
        df = (
            df.join(lookup, on=["wgs_tumour_label", "sample"], how="left")
            .join(tumours, on="wgs_tumour_label", how="left")
            .with_columns((pl.col("cn_major") + pl.col("cn_minor")).alias("cn_total"))
            .select(list(SCHEMA))
        )

        inverted = df.filter(pl.col("cn_major") < pl.col("cn_minor")).height
        if inverted:
            log(f"note: {inverted} row(s) with cn_major < cn_minor (kept as on source)")
        per_sample = (
            df.group_by("tumour_id", "sample")
            .agg(pl.concat_str("chr", "start", "end", separator=":").sort().alias("segs"))
            .group_by("tumour_id").agg(pl.col("segs").n_unique().alias("n_layouts"))
            .filter(pl.col("n_layouts") > 1)
        )
        if per_sample.height:
            log(f"warning: samples don't share one segmentation in {per_sample.height} "
                f"tumour(s): {per_sample['tumour_id'].to_list()[:10]}")
        dupes = (
            df.group_by("wgs_tumour_label", "sample", "chr", "start")
            .len().filter(pl.col("len") > 1)
        )
        if dupes.height:
            raise ValueError(
                f"{TABLE}: {dupes.height} duplicated (sample, segment) rows; "
                f"first: {dupes.head(5).to_dicts()}"
            )

    df = df.sort("tumour_id", "sample", "chr", "start")
    out = canonical_output_path(TABLE)
    df.write_parquet(out)
    log(f"wrote {out} ({df.height} rows × {df.width} cols)")
    return df


if __name__ == "__main__":
    build()
