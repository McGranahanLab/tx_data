"""Build cohort_metadata → canonical Parquet: one row per (cohort × field) from
references/cohorts/cohort_metadata.yml (Piotr, 2026-11-19), e.g. the genome build of each cohort.
Fields ending `_note` become the `note` of the field they extend.

Computed fields (Piotr, 2026-11-25), when `sample_clone_counts` / `wgs_sample_clone_counts` are built:
per cohort, the number of single-sample, single-clone and single-sample-and-single-clone tumours
(`value`) with their tumour ids (`note`). Built after those tables (build_all order)."""

from __future__ import annotations

import polars as pl
import yaml

from tx_data.builds._base import canonical_output_path, log
from tx_data.builds.sample_clone_counts import CATEGORIES
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
    rows += _structure_rows()
    df = pl.DataFrame(rows, schema=SCHEMA)
    out = canonical_output_path(TABLE)
    df.write_parquet(out)
    log(f"wrote {out} ({df.height} rows, cohorts {sorted(doc)})")
    return df


COUNT_TABLES = {"wes": "sample_clone_counts", "wgs": "wgs_sample_clone_counts"}


def _structure_rows() -> list[dict]:
    """Single-sample / single-clone tumours per cohort, from the sample_clone_counts tables."""
    rows = []
    for cohort, table in COUNT_TABLES.items():
        p = canonical_output_path(table)
        if not p.is_file():
            log(f"warning: {table} not built; cohort_metadata without {cohort} single-sample/clone counts")
            continue
        df = pl.read_parquet(p)
        for cat in [c for c in CATEGORIES.values() if c != "multi_sample_multi_clone"] + ["clones_unknown"]:
            ids = df.filter(pl.col("category") == cat)["tumour_id"].sort().to_list()
            if cat == "clones_unknown" and not ids:
                continue
            rows.append({"cohort": cohort, "field": f"{cat}_tumours", "value": str(len(ids)),
                         "note": f"computed from {table} ({df.height} tumours): " + (", ".join(ids) or "none")})
    return rows


if __name__ == "__main__":
    build()
