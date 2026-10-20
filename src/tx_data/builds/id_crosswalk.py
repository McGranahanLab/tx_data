"""Build id_crosswalk → canonical Parquet.

Patient-level LTX ↔ CRUK identifier mapping. See ``catalog/id_crosswalk.yml``.
"""

from __future__ import annotations

import re

import polars as pl

from tx_data.builds._base import canonical_output_path, log
from tx_data.sources import resolve_source

TABLE = "id_crosswalk"

_LTX_RE = re.compile(r"^LTX\d+$")
_CRUK_RE = re.compile(r"^CRUK\d+$")


def build() -> pl.DataFrame:
    src = resolve_source(TABLE)
    log(f"read {src}")

    df = pl.read_csv(src, null_values=["NA", "", "nan"])

    # Canonical column names.
    df = df.rename({"PatientID": "ltx_patient_id", "CRUK_ID": "cruk_patient_id"})

    # Trim any stray whitespace on IDs (defensive — R write.csv can leave
    # trailing spaces in some releases).
    df = df.with_columns(
        pl.col("ltx_patient_id").str.strip_chars(),
        pl.col("cruk_patient_id").str.strip_chars(),
    )

    # Validate shape. Fail loudly if anything is off — this table is the
    # backbone for every CRUK↔LTX translation downstream.
    ltx_bad = df.filter(~pl.col("ltx_patient_id").str.contains(r"^LTX\d+$"))
    if ltx_bad.height:
        raise ValueError(
            f"id_crosswalk: {ltx_bad.height} rows have non-canonical ltx_patient_id.\n"
            f"First offenders: {ltx_bad.head(5).to_dicts()}"
        )
    cruk_bad = df.filter(~pl.col("cruk_patient_id").str.contains(r"^CRUK\d+$"))
    if cruk_bad.height:
        raise ValueError(
            f"id_crosswalk: {cruk_bad.height} rows have non-canonical cruk_patient_id.\n"
            f"First offenders: {cruk_bad.head(5).to_dicts()}"
        )

    if df["ltx_patient_id"].n_unique() != df.height:
        raise ValueError("id_crosswalk: ltx_patient_id column has duplicates")
    if df["cruk_patient_id"].n_unique() != df.height:
        raise ValueError("id_crosswalk: cruk_patient_id column has duplicates")

    df = df.select(["ltx_patient_id", "cruk_patient_id"]).sort("ltx_patient_id")

    out = canonical_output_path(TABLE)
    df.write_parquet(out)
    log(f"wrote {out} ({df.height} rows × {df.width} cols)")
    return df


if __name__ == "__main__":
    build()
