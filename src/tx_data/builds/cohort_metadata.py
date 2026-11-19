"""Build cohort_metadata → canonical Parquet: one row per (cohort × field) from
references/cohorts/cohort_metadata.yml (Piotr, 2026-11-19), e.g. the genome build of each cohort.
Fields ending `_note` become the `note` of the field they extend."""

from __future__ import annotations

import polars as pl
import yaml

from tx_data.builds._base import canonical_output_path, log
from tx_data.sources import resolve_source

TABLE = "cohort_metadata"
SCHEMA = {"cohort": pl.String, "field": pl.String, "value": pl.String, "note": pl.String}


def build() -> pl.DataFrame:
    src = resolve_source(TABLE)
    doc = yaml.safe_load(src.read_text())
    rows = []
    for cohort, fields in doc.items():
        for k, v in fields.items():
            if k.endswith("_note"):
                continue
            rows.append({"cohort": cohort, "field": k, "value": str(v), "note": fields.get(f"{k}_note")})
    df = pl.DataFrame(rows, schema=SCHEMA)
    out = canonical_output_path(TABLE)
    df.write_parquet(out)
    log(f"wrote {out} ({df.height} rows, cohorts {sorted(doc)})")
    return df


if __name__ == "__main__":
    build()
