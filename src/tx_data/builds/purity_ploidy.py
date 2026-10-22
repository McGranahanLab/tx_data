"""Build purity_ploidy → canonical Parquet.

Source is a flat CSV keyed by `sample_name_hash`. The build script looks
each hash up in `sample_hash_map` to derive canonical `tumour_id` + `lesion`
so the table joins directly against every other per-tumour table while
preserving the original hash for lossless joins against muttable + kallisto.

Samples not present in the mapping (typically normals / buffy coats) are
dropped — purity/ploidy only applies to tumour tissue.
"""

from __future__ import annotations

import polars as pl

from tx_data.builds._base import canonical_output_path, log
from tx_data.ids import apply_convention_scalar
from tx_data.sample_map import extract_hash_suffix, extract_sample_stem, lookup
from tx_data.sources import resolve_source

TABLE = "purity_ploidy"


def build() -> pl.DataFrame:
    src = resolve_source(TABLE)
    log(f"read {src}")

    df = pl.read_csv(src, null_values=["NA", "", "nan"], infer_schema_length=10000)

    expected = {"sample_name_hash", "purity", "ploidy"}
    missing = expected - set(df.columns)
    if missing:
        raise ValueError(
            f"purity_ploidy source missing expected columns {sorted(missing)}; "
            f"got {df.columns!r}"
        )

    unique_hashes = df.get_column("sample_name_hash").unique().to_list()
    tumour_ids: dict[str, str] = {}
    lesions: dict[str, str] = {}
    for sh in unique_hashes:
        try:
            extract_hash_suffix(sh)
        except ValueError as exc:
            raise ValueError(
                f"purity_ploidy: cannot parse sample_name_hash {sh!r}: {exc}"
            ) from exc
        rec = lookup(sh)
        if rec is None:
            continue
        tumour_ids[sh] = apply_convention_scalar(rec.patient_tumour)
        lesions[sh] = rec.lesion

    df = df.with_columns(
        [
            pl.col("sample_name_hash").replace_strict(
                tumour_ids, default=None, return_dtype=pl.String
            ).alias("tumour_id"),
            pl.col("sample_name_hash").replace_strict(
                lesions, default=None, return_dtype=pl.String
            ).alias("lesion"),
            pl.col("sample_name_hash")
            .map_elements(extract_sample_stem, return_dtype=pl.String)
            .alias("sample_stem"),
            pl.col("purity").cast(pl.Float64),
            pl.col("ploidy").cast(pl.Float64),
        ]
    )

    unmapped = df.filter(pl.col("tumour_id").is_null())
    if unmapped.height:
        log(
            f"skipping {unmapped.height} sample(s) not in sample_hash_map "
            f"(likely non-tumour); first: "
            f"{unmapped.get_column('sample_name_hash').head(3).to_list()}"
        )
        df = df.filter(pl.col("tumour_id").is_not_null())

    dupes = df.group_by("sample_name_hash").len().filter(pl.col("len") > 1)
    if dupes.height:
        raise ValueError(
            f"purity_ploidy: {dupes.height} duplicate sample_name_hash rows. "
            f"First offenders: {dupes.head(5).to_dicts()}"
        )

    out_df = df.select(
        ["tumour_id", "lesion", "sample_stem", "sample_name_hash", "purity", "ploidy"]
    ).sort(["tumour_id", "lesion", "sample_name_hash"])

    out = canonical_output_path(TABLE)
    out_df.write_parquet(out)
    log(f"wrote {out} ({out_df.height} rows × {out_df.width} cols)")
    return out_df


if __name__ == "__main__":
    build()
