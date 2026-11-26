"""Shared reading for the per-tumour WGS sources (CONIPHER run folders, ALPACA tumour_results).

Each WGS table below is one file (or two) per tumour, located by `wgs_tumour_label`
(`LTX0001`, `LTX0104_tumour1`); tumours come from the built `wgs_battenberg_inventory`.
The canonical `tumour_id` / `patient_id` always come from the inventory, never from the
file, so a label-form id inside a file can't leak through.
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import polars as pl

from tx_data import wgs_samples
from tx_data.builds._base import Issues, canonical_output_path, log
from tx_data.sources import resolve_source

ID_COLS = {"tumour_id": pl.String, "patient_id": pl.String, "wgs_tumour_label": pl.String}


def build_per_tumour(table: str, schema: dict, read_one: Callable[[str, Path], pl.DataFrame],
                     sort: tuple[str, ...]) -> pl.DataFrame:
    """Read every tumour's file with `read_one(label, path)` (frame without the id columns),
    attach tumour_id / patient_id / label, check the schema, write the canonical parquet."""
    wgs_samples.clear_cache()
    tumours = wgs_samples.tumours()
    issues = Issues(table)
    tid = dict(tumours.select("wgs_tumour_label", "tumour_id").iter_rows())
    frames = []
    for label in tumours["wgs_tumour_label"].to_list():
        path = resolve_source(table, wgs_tumour_label=label)
        if not path.is_file():
            issues.add(tid[label], "skipped", f"no {path.name} for WGS label {label}", path)
            continue
        try:
            frames.append(read_one(label, path).with_columns(pl.lit(label).alias("wgs_tumour_label")))
        except Exception as e:
            issues.add(tid[label], "skipped", f"unreadable: {type(e).__name__}: {e}", path)
    log(f"{table}: read {len(frames)}/{tumours.height} WGS tumours")
    if frames:
        df = (pl.concat(frames, how="diagonal_relaxed")
                .join(tumours, on="wgs_tumour_label", how="left")
                .select([pl.col(c).cast(t) for c, t in schema.items()]))
    else:
        df = pl.DataFrame(schema=schema)
    df = df.sort(*sort, nulls_last=True)
    out = canonical_output_path(table)
    df.write_parquet(out)
    log(f"wrote {out} ({df.height} rows × {df.width} cols)")
    issues.write()
    return df


def read_tsv(path: Path) -> pl.DataFrame:
    """Tab-separated, every column as text (cast by the caller); R's NA/NaN read as null."""
    return pl.read_csv(path, separator="\t", infer_schema_length=0, null_values=["NA", "NaN", ""])


def require(df: pl.DataFrame, cols, path: Path) -> None:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(f"{path}: missing required columns {missing}")


def r_int(df: pl.DataFrame, cols, path: Path) -> pl.DataFrame:
    """Integer columns written by R, which may use scientific notation (`2e+07` for 20000000, Nemo
    2026-11-24): parse as float, check every value is whole, cast to Int64. A fractional value fails."""
    out = df.with_columns(pl.col(c).cast(pl.Float64) for c in cols)
    for c in cols:
        bad = out.filter(pl.col(c).is_not_null() & (pl.col(c) != pl.col(c).round(0)))[c]
        if bad.len():
            raise ValueError(f"{path}: column {c} has non-integer values {bad.unique().to_list()[:5]}")
    return out.with_columns(pl.col(c).cast(pl.Int64) for c in cols)


def r_bool(col: str) -> pl.Expr:
    """R logicals (`TRUE`/`FALSE`, also `1`/`0`); anything else → null."""
    v = pl.col(col).str.to_uppercase()
    return (pl.when(v.is_in(["TRUE", "T", "1"])).then(True)
              .when(v.is_in(["FALSE", "F", "0"])).then(False)
              .otherwise(None).alias(col))


def sample_hash(label: str) -> pl.Expr:
    """Inventory `sample_name_hash` for a source SAMPLE (any known short form), within the tumour."""
    return pl.col("sample").map_elements(lambda s: wgs_samples.match_sample(label, s),
                                         return_dtype=pl.String).alias("sample_name_hash")
