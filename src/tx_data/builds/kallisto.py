"""Build kallisto → canonical Parquet.

Two source parquets — counts + TPM — are wide (rows = genes, columns = a
`gene_id` column + N sample columns). The builder:

1. Reads both, verifies gene and sample sets match.
2. Melts each wide table to long (gene_id, sample_name_hash, value).
3. Joins them into a single (sample_name_hash × gene_id) table with counts+tpm.
4. Parses each sample name into (patient_id, sample_type, region_code,
   is_normal, tumour_ordinal, region_ordinal) — matches muttable.sample_name_hash.
5. Writes canonical Parquet + a paired-normals view keyed by patient_id.
"""

from __future__ import annotations

import re

import polars as pl

from tx_data.builds._base import canonical_output_path, log
from tx_data.sources import resolve_source

TABLE = "kallisto"

# LTX{patient}_{SU|BS}_{T{n}-R{n} | N{n}}--{hash}
_SAMPLE_RE = re.compile(
    r"^(?P<patient_id>LTX\d+)"
    r"_(?P<sample_type>[A-Z]+)"
    r"_(?P<region_code>[^-]+(?:-[^-]+)?)"
    r"--(?P<hash>[0-9a-f]+)$"
)


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


def _parse_samples(df: pl.DataFrame) -> pl.DataFrame:
    """Add parsed sample-name columns to the long table."""
    unique_samples = df.select("sample_name_hash").unique()
    parsed_rows = []
    for row in unique_samples.iter_rows(named=True):
        s = row["sample_name_hash"]
        m = _SAMPLE_RE.match(s)
        if not m:
            parsed_rows.append(
                dict(
                    sample_name_hash=s,
                    patient_id=None, sample_type=None, region_code=None,
                    is_normal=None, tumour_ordinal=None, region_ordinal=None,
                )
            )
            continue
        region_code = m.group("region_code")
        is_normal = region_code.startswith("N")
        tumour_ordinal = None
        region_ordinal = None
        if not is_normal and "-" in region_code:
            t_part, r_part = region_code.split("-", 1)
            try:
                tumour_ordinal = int(t_part[1:])
                region_ordinal = int(r_part[1:])
            except ValueError:
                pass
        parsed_rows.append(
            dict(
                sample_name_hash=s,
                patient_id=m.group("patient_id"),
                sample_type=m.group("sample_type"),
                region_code=region_code,
                is_normal=is_normal,
                tumour_ordinal=tumour_ordinal,
                region_ordinal=region_ordinal,
            )
        )
    parsed = pl.DataFrame(parsed_rows)
    return df.join(parsed, on="sample_name_hash", how="left")


def build() -> pl.DataFrame:
    src_counts = resolve_source("kallisto_counts")
    src_tpm = resolve_source("kallisto_tpm")
    log(f"read {src_counts}")
    log(f"read {src_tpm}")

    counts_long = _melt_wide(src_counts, "counts")
    tpm_long = _melt_wide(src_tpm, "tpm")

    # Verify the two sources agree on genes and samples.
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
    long = _parse_samples(long)

    # Reorder for readability.
    long = long.select(
        [
            "sample_name_hash", "patient_id", "sample_type", "region_code",
            "is_normal", "tumour_ordinal", "region_ordinal",
            "gene_id", "counts", "tpm",
        ]
    )

    out = canonical_output_path(TABLE)
    long.write_parquet(out)
    log(f"wrote {out} ({long.height} rows × {long.width} cols)")

    # Paired-normals view: for patients with ≥1 normal sample, per-gene TPM.
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
