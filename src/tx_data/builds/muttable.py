"""Build muttable → canonical Parquet.

Applies:
  - value normalisation:
      * ``patient_tumour``: ``LTX0001_tumour1`` → ``LTX0001-Tumour1``
      * ``mutation_cluster``: float → ``clone{N}`` string
      * ``chr``: strip ``chr`` prefix (``chr17`` → ``17``), coerce to String
  - column renames to canonical tx_data schema (see catalog/muttable.yml):
      * ``patient_tumour`` → ``tumour_id``
      * ``mutation_cluster`` → ``clone``
      * ``var`` → ``alt``

Source: TSV gz. Output: single Parquet at `data/muttable.parquet`.
"""

from __future__ import annotations

import polars as pl

from tx_data.builds._base import canonical_output_path, log
from tx_data.normalize import (
    canonical_chr_expr,
    canonical_clone_expr,
    canonical_tumour_id_expr,
)
from tx_data.sources import resolve_source

TABLE = "muttable"


def build() -> pl.DataFrame:
    src = resolve_source(TABLE)
    log(f"read {src}")

    # Force `chr` to String at READ time. If we let polars infer, the first
    # 10 000 rows may be all autosomes → i64 inferred → parse fails when it
    # hits X / Y / MT further into the file.
    df = pl.read_csv(
        src,
        separator="\t",
        null_values=["NA", ""],
        infer_schema_length=10000,
        schema_overrides={"chr": pl.String},
    )

    df = df.with_columns(
        [
            # Canonicalise tumour_id from muttable's `patient_tumour` form.
            canonical_tumour_id_expr("patient_tumour"),
            # Canonicalise clone (was float mutation_cluster).
            canonical_clone_expr("mutation_cluster"),
            # Canonicalise chromosome (strip `chr` prefix, coerce to string).
            canonical_chr_expr("chr"),
        ]
    )

    # Source has BOTH a numeric `tumour_id` (the ordinal `1`, `2`, ...) and a
    # string `patient_tumour` (`LTX0001_tumour1`). To make the canonical
    # `tumour_id` the string form (matching every other per-tumour table), we
    # first free the name by promoting the numeric column to `tumour_ordinal`.
    renames: dict[str, str] = {"mutation_cluster": "clone", "var": "alt"}
    if "tumour_id" in df.columns:
        renames["tumour_id"] = "tumour_ordinal"
    renames["patient_tumour"] = "tumour_id"
    df = df.rename(renames)

    out = canonical_output_path(TABLE)
    df.write_parquet(out)
    log(f"wrote {out} ({df.height} rows × {df.width} cols)")
    return df


if __name__ == "__main__":
    build()
