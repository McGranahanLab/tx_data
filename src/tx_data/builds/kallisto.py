"""Build kallisto → canonical Parquet.

Two source parquets — counts + TPM — are wide (rows = genes, columns = a
`gene_id` column + N sample columns). The builder:

1. Reads both, verifies gene and sample sets match.
2. Melts each wide table to long (gene_id, sample_name_hash, value).
3. Joins them into a single (sample_name_hash × gene_id) table with counts+tpm.
4. Classifies each sample as normal vs tumour via the `_N<n>` substring
   (`_N01`, `_N1`, …) in the sample identifier — the reliable normal
   signal per User (2026-10-21).
5. Looks every sample up in `sample_hash_map` by its hash suffix. The
   mapping is the source of truth and covers every tumour sample
   category the release recognises (primary multi-region, lymph-node
   metastases, distant metastases, …), per User (2026-10-23). This
   replaces the earlier stem-regex parsing — real samples carry CODE
   segments the regex can't predict (`_LN02`, …).
6. For samples absent from the mapping (RNA-only patients, mapping lag,
   …) the build keeps the row and leaves tumour-level columns null, with
   a summary log line. Normals pass through on the `_N<n>` substring
   even when absent from the mapping.
7. Emits `sample_stem` (part before `--`) as the explicit RNA↔DNA join
   key for downstream per-biopsy matching with DNA tables.
8. Writes canonical Parquet + a paired-normals view keyed by patient_id.
"""

from __future__ import annotations

import re

import polars as pl

from tx_data.builds._base import canonical_output_path, log
from tx_data.ids import apply_convention_scalar
from tx_data.sample_map import extract_sample_stem, lookup
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


def _patient_from_stem(stem: str) -> str | None:
    m = _PATIENT_RE.match(stem)
    return m.group("patient_id") if m else None


def _parse_sample(sample_name_hash: str) -> dict:
    """Return canonical identifier fields for one kallisto sample.

    Primary strategy: look the hash up in `sample_hash_map`. The mapping
    is the source of truth for which (patient_tumour, lesion) each hash
    corresponds to — no name-token parsing required.

    Fallbacks for samples absent from the mapping:
      * `_N<n>` substring → treated as a normal (patient_id from stem
        prefix, tumour-level columns null);
      * otherwise → patient_id from stem prefix if parseable, else null,
        and tumour-level columns null. The build logs how many such
        samples it saw.
    """
    try:
        stem = extract_sample_stem(sample_name_hash)
    except ValueError as exc:
        raise ValueError(
            f"kallisto: cannot parse sample_name_hash {sample_name_hash!r}: {exc}"
        ) from exc

    is_normal = bool(_NORMAL_RE.search(stem))
    patient_id_from_stem = _patient_from_stem(stem)

    rec = lookup(sample_name_hash)
    if rec is not None:
        return dict(
            sample_name_hash=sample_name_hash,
            sample_stem=stem,
            patient_id=rec.patient_id,
            patient_tumour=apply_convention_scalar(rec.patient_tumour),
            lesion=rec.lesion,
            histology_group=rec.histology_group,
            is_normal=is_normal,
            in_sample_map=True,
        )

    # Not in the mapping. Pass the sample through with whatever we can
    # recover from the stem — never crash on a shape we don't recognise.
    return dict(
        sample_name_hash=sample_name_hash,
        sample_stem=stem,
        patient_id=patient_id_from_stem,
        patient_tumour=None,
        lesion=None,
        histology_group=None,
        is_normal=is_normal,
        in_sample_map=False,
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
            "lesion", "histology_group", "is_normal", "in_sample_map",
            "gene_id", "counts", "tpm",
        ]
    )

    unmapped = (
        long.filter(pl.col("in_sample_map") == False)  # noqa: E712
        .select(["sample_name_hash", "sample_stem", "is_normal"])
        .unique()
    )
    if unmapped.height:
        n_normals = unmapped.filter(pl.col("is_normal") == True).height  # noqa: E712
        n_tumour_ish = unmapped.height - n_normals
        log(
            f"kallisto: {unmapped.height} sample(s) not in sample_hash_map "
            f"({n_normals} normal / {n_tumour_ish} tumour-ish); first stems: "
            f"{unmapped.get_column('sample_stem').head(5).to_list()}. "
            f"Tumour-level columns left null for these."
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
