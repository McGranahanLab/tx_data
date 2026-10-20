"""Build clinical → canonical Parquet, plus a per-tumour unpivoted view.

Emits two files:
  - `data/clinical.parquet` — one row per patient (source shape).
  - `data/clinical_per_tumour.parquet` — one row per (patient, tumour_ordinal),
    with lesion1_/lesion2_/histology1_/histology2_ column pairs unpivoted.
"""

from __future__ import annotations

import polars as pl

from tx_data.builds._base import canonical_output_path, log
from tx_data.sources import resolve_source

TABLE = "clinical"

_LESION_STEMS = [
    ("Lesion{n}_site_central.reviewed", "lesion_site"),
    # PRIMARY histology column for our analyses is the detailed one
    # (histology1_central.reviewed) per User 2026-09-01. The `group` one
    # is kept as a coarser rollup.
    ("histology{n}_central.reviewed", "histology"),
    ("histology{n}_group_central.reviewed", "histology_group"),
    ("pTStage_v8_lesion{n}_central.reviewed", "pT_stage"),
    ("pNStage_lesion{n}_central.reviewed", "pN_stage"),
    # PRIMARY staging column: pTNMStage_v8_lesion1_central.reviewed → pTNM_stage.
    ("pTNMStage_v8_lesion{n}_central.reviewed", "pTNM_stage"),
    ("SizePath_lesion{n}_central.reviewed", "size_path_mm"),
    ("PathPleuInv_lesion{n}_central.reviewed", "pleural_invasion"),
    ("margin_status_lesion{n}_central.reviewed", "margin_status"),
]


def _unpivot_lesions(df: pl.DataFrame) -> pl.DataFrame:
    """Two-lesion column sets → tall (patient × tumour_ordinal) table."""
    out_parts = []
    for n in (1, 2):
        cols_map = {}
        for pattern, tidy in _LESION_STEMS:
            src_col = pattern.format(n=n)
            if src_col in df.columns:
                cols_map[src_col] = tidy
        # Build the per-lesion slice.
        keep = ["Patient_ID"] + list(cols_map.keys())
        part = df.select([c for c in keep if c in df.columns])
        # Rename source → tidy, add tumour_ordinal.
        renames = {src: tidy for src, tidy in cols_map.items()}
        part = part.rename(renames).with_columns(pl.lit(n).alias("tumour_ordinal"))
        out_parts.append(part)

    tall = pl.concat(out_parts, how="diagonal_relaxed")

    # Drop rows where the whole tumour block is empty (patient has no lesionN).
    non_id_cols = [c for c in tall.columns if c not in ("Patient_ID", "tumour_ordinal")]
    if non_id_cols:
        keep_mask = pl.any_horizontal([pl.col(c).is_not_null() for c in non_id_cols])
        tall = tall.filter(keep_mask)

    # Add canonical tumour_id via native format expression.
    tall = tall.with_columns(
        pl.format("{}-Tumour{}", pl.col("Patient_ID"), pl.col("tumour_ordinal"))
        .alias("tumour_id")
    )
    return tall


def build() -> tuple[pl.DataFrame, pl.DataFrame]:
    src = resolve_source(TABLE)
    log(f"read {src}")
    df = pl.read_csv(src, separator="\t", null_values=["NA", ""], infer_schema_length=10000)

    per_patient_out = canonical_output_path(TABLE)
    df.write_parquet(per_patient_out)
    log(f"wrote {per_patient_out} ({df.height} rows × {df.width} cols)")

    per_tumour = _unpivot_lesions(df)
    per_tumour_out = canonical_output_path(f"{TABLE}_per_tumour")
    per_tumour.write_parquet(per_tumour_out)
    log(f"wrote {per_tumour_out} ({per_tumour.height} rows × {per_tumour.width} cols)")

    return df, per_tumour


if __name__ == "__main__":
    build()
