"""Build kallisto → canonical Parquet.

Two source parquets — counts + TPM — are wide (rows = genes, columns = a
`gene_id` column + N sample columns). The builder:

1. Reads both, verifies gene and sample sets match.
2. Melts each wide table to long (gene_id, sample_name_hash, value).
3. Joins them into a single (sample_name_hash × gene_id) table with counts+tpm.
4. Classifies each sample as normal vs tumour via the `_N<n>` substring
   (`_N01`, `_N1`, …) in the sample identifier — the reliable normal signal
   per User (2026-10-21).
5. For tumour samples, looks up canonical identifiers in `sample_hash_map`:
   `patient_id`, `patient_tumour`, `lesion`, `histology_group`.
6. For normal samples, extracts `patient_id` from the leading
   `(LTX|CRUK)<digits>` prefix (the one piece of the name we're
   guaranteed) and leaves tumour-level columns null.
7. Writes canonical Parquet + a paired-normals view keyed by patient_id.
"""

from __future__ import annotations

import re

import polars as pl

from tx_data.builds._base import canonical_output_path, log
from tx_data.ids import apply_convention_scalar
from tx_data.sample_map import extract_hash_suffix, lookup
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
    """Return canonical identifier fields for one sample."""
    try:
        extract_hash_suffix(sample_name_hash)
    except ValueError as exc:
        raise ValueError(
            f"kallisto: cannot parse sample_name_hash {sample_name_hash!r}: {exc}"
        ) from exc

    is_normal = bool(_NORMAL_RE.search(sample_name_hash))

    if is_normal:
        pm = _PATIENT_RE.match(sample_name_hash)
        if pm is None:
            raise ValueError(
                f"kallisto: normal sample missing (LTX|CRUK)<digits> prefix: "
                f"{sample_name_hash!r}"
            )
        return dict(
            sample_name_hash=sample_name_hash,
            is_normal=True,
            patient_id=pm.group("patient_id"),
            patient_tumour=None,
            lesion=None,
            histology_group=None,
        )

    rec = lookup(sample_name_hash)
    if rec is None:
        # Not normal (no `_N<n>`) and not in the map: flag it clearly rather
        # than silently drop — unexpected shape wants a human decision.
        raise KeyError(
            f"kallisto: tumour sample not in sample_hash_map (and not a normal): "
            f"{sample_name_hash!r}"
        )
    return dict(
        sample_name_hash=sample_name_hash,
        is_normal=False,
        patient_id=rec.patient_id,
        patient_tumour=apply_convention_scalar(rec.patient_tumour),
        lesion=rec.lesion,
        histology_group=rec.histology_group,
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
            "sample_name_hash", "patient_id", "patient_tumour", "lesion",
            "histology_group", "is_normal",
            "gene_id", "counts", "tpm",
        ]
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
