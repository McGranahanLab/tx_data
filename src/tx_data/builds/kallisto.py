"""Build kallisto → canonical Parquet.

Two source parquets — counts + TPM — are wide (rows = genes, columns = a
`gene_id` column + N sample columns). The builder:

1. Reads both, verifies gene and sample sets match.
2. Melts each wide table to long (gene_id, sample_name_hash, value).
3. Joins them into a single (sample_name_hash × gene_id) table with counts+tpm.
4. Classifies each sample as normal vs tumour via the `_N<n>` substring
   (`_N01`, `_N1`, …) in the sample identifier — the reliable normal
   signal per User (2026-10-21).
5. Parses the sample **stem** (everything before `--`) to recover
   `patient_id` / `tumour_ordinal` / `region_ordinal`. RNA and DNA share
   the stem but carry different hashes, so the stem — not the hash — is
   the RNA↔DNA join key. The stem pattern is reliable for kallisto per
   User (2026-10-22).
6. For tumour samples, derives `patient_tumour` (via
   `apply_convention_scalar`) and `lesion` from the parsed stem; looks
   up `histology_group` by `patient_tumour` against `sample_hash_map`.
7. Writes canonical Parquet + a paired-normals view keyed by patient_id.
"""

from __future__ import annotations

import re

import polars as pl

from tx_data.builds._base import canonical_output_path, log
from tx_data.ids import apply_convention_scalar
from tx_data.sample_map import (
    extract_sample_stem,
    lookup_tumour_histology,
    parse_tumour_stem,
)
from tx_data.sources import resolve_source

TABLE = "kallisto"

_NORMAL_RE = re.compile(r"_N\d+")
_PATIENT_RE = re.compile(r"^(?P<patient_id>(?:LTX|CRUK)\d+)_")


def _melt_wide(path, value_name: str) -> pl.DataFrame:
    df = pl.read_parquet(path)
    if "gene_id" not in df.columns:
        raise ValueError(f"{path}: expected `gene_id` column, got {df.columns[:5]}…")
    sample_cols = [c for c in df.columns if c != "gene_id"]
    return df.unpivot(
        index="gene_id",
        on=sample_cols,
        variable_name="sample_name_hash",
        value_name=value_name,
    )


def _parse_sample(sample_name_hash: str) -> dict:
    """Return canonical identifier fields for one kallisto sample."""
    try:
        stem = extract_sample_stem(sample_name_hash)
    except ValueError as exc:
        raise ValueError(
            f"kallisto: cannot parse sample_name_hash {sample_name_hash!r}: {exc}"
        ) from exc

    is_normal = bool(_NORMAL_RE.search(stem))

    if is_normal:
        pm = _PATIENT_RE.match(stem)
        if pm is None:
            raise ValueError(
                f"kallisto: normal sample missing (LTX|CRUK)<digits> prefix: "
                f"{sample_name_hash!r}"
            )
        return dict(
            sample_name_hash=sample_name_hash,
            sample_stem=stem,
            is_normal=True,
            patient_id=pm.group("patient_id"),
            patient_tumour=None,
            lesion=None,
            histology_group=None,
        )

    parts = parse_tumour_stem(stem)
    if parts is None:
        raise ValueError(
            f"kallisto: tumour sample stem does not match "
            f"'(LTX|CRUK)#_<TYPE>_T#-R#': {sample_name_hash!r}"
        )
    muttable_tumour_id = f"{parts.patient_id}_tumour{parts.tumour_ordinal}"
    return dict(
        sample_name_hash=sample_name_hash,
        sample_stem=stem,
        is_normal=False,
        patient_id=parts.patient_id,
        patient_tumour=apply_convention_scalar(muttable_tumour_id),
        lesion=f"lesion{parts.region_ordinal}",
        histology_group=lookup_tumour_histology(muttable_tumour_id),
    )


def _annotate_samples(df: pl.DataFrame) -> pl.DataFrame:
    unique_samples = df.get_column("sample_name_hash").unique().to_list()
    rows = [_parse_sample(s) for s in unique_samples]
    annotated = pl.DataFrame(rows)
    return df.join(annotated, on="sample_name_hash", how="left")


def build() -> pl.DataFrame:
    src_counts = resolve_source("kallisto_counts")
    src_tpm = resolve_source("kallisto_tpm")
    log(f"read {src_counts}")
    log(f"read {src_tpm}")

    counts_long = _melt_wide(src_counts, "counts")
    tpm_long = _melt_wide(src_tpm, "tpm")

    genes_c = set(counts_long["gene_id"].unique().to_list())
    genes_t = set(tpm_long["gene_id"].unique().to_list())
    if genes_c != genes_t:
        raise ValueError(
            f"kallisto counts vs tpm gene sets differ: "
            f"{len(genes_c - genes_t)} counts-only, {len(genes_t - genes_c)} tpm-only"
        )
    samples_c = set(counts_long["sample_name_hash"].unique().to_list())
    samples_t = set(tpm_long["sample_name_hash"].unique().to_list())
    if samples_c != samples_t:
        raise ValueError(
            f"kallisto counts vs tpm sample sets differ: "
            f"{len(samples_c - samples_t)} counts-only, {len(samples_t - samples_c)} tpm-only"
        )

    long = counts_long.join(tpm_long, on=["gene_id", "sample_name_hash"], how="full")
    long = _annotate_samples(long)

    long = long.select(
        [
            "sample_name_hash", "sample_stem", "patient_id", "patient_tumour",
            "lesion", "histology_group", "is_normal",
            "gene_id", "counts", "tpm",
        ]
    )

    missing_histology = long.filter(
        (pl.col("is_normal") == False)  # noqa: E712
        & pl.col("histology_group").is_null()
    )
    if missing_histology.height:
        unknown_tumours = (
            missing_histology.get_column("patient_tumour").unique().to_list()
        )
        log(
            f"kallisto: {len(unknown_tumours)} tumour(s) have no matching "
            f"patient_tumour in sample_hash_map (likely RNA-only patients); "
            f"histology_group left null for {missing_histology.height} rows. "
            f"First: {unknown_tumours[:5]}"
        )

    out = canonical_output_path(TABLE)
    long.write_parquet(out)
    log(f"wrote {out} ({long.height} rows × {long.width} cols)")

    normals = long.filter(pl.col("is_normal") == True)  # noqa: E712
    if normals.height > 0:
        paired = (
            normals.group_by(["patient_id", "gene_id"])
            .agg(
                [
                    pl.col("tpm").mean().alias("tpm_normal"),
                    pl.col("sample_name_hash").first().alias("normal_sample_name_hash"),
                ]
            )
        )
        paired_out = canonical_output_path("kallisto_paired_normals")
        paired.write_parquet(paired_out)
        log(f"wrote {paired_out} ({paired.height} rows) — {paired['patient_id'].n_unique()} patients")
    else:
        log("no normal samples found; skipping paired-normals view")

    return long


if __name__ == "__main__":
    build()
