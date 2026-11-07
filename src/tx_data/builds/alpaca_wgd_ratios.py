"""Build alpaca_wgd_ratios → canonical Parquet.

Reads every per-tumour ALPACA `wgd_ratio_scores.csv` (output of
`alpaca wgd`, one row per clone) and concatenates them into one long
table keyed by (tumour_id, clone). Adds `is_mrca` (parent == 'diploid').

Required source columns: tumour_id, clone, parent, ratio_A. The rest
(ratio_B, ratio_total, segments_used_*, weight_sum_*) are filled with null
when a file lacks them: older ALPACA versions wrote only the first four
(catalog/alpaca_wgd_ratios.yml open_questions).
"""

from __future__ import annotations

from pathlib import Path

import polars as pl

from tx_data.builds._base import canonical_output_path, log
from tx_data.sources import resolve_source

TABLE = "alpaca_wgd_ratios"

REQUIRED = {
    "tumour_id": pl.String,
    "clone": pl.String,
    "parent": pl.String,
    "ratio_A": pl.Float64,
}
OPTIONAL = {
    "ratio_B": pl.Float64,
    "ratio_total": pl.Float64,
    "segments_used_A": pl.Int64,
    "segments_used_B": pl.Int64,
    "segments_used_total": pl.Int64,
    "weight_sum_A": pl.Float64,
    "weight_sum_B": pl.Float64,
    "weight_sum_total": pl.Float64,
}
SCHEMA = {**REQUIRED, **OPTIONAL}


def _discover() -> list[tuple[str, Path]]:
    """Return ``[(tumour_id_from_dir, path), ...]`` for every tumour with a file."""
    probe = resolve_source(TABLE, tumour_id="__probe__")
    root = probe.parent.parent              # `.../ALPACA/output/tumour_results/`
    if not root.is_dir():
        log(f"{TABLE}: tumour_results root does not exist: {root}")
        return []
    return [
        (sub.name, sub / probe.name)
        for sub in sorted(root.iterdir())
        if sub.is_dir() and (sub / probe.name).is_file()
    ]


def read_wgd_ratio_file(path: Path, tumour_dir: str) -> tuple[pl.DataFrame, list[str]]:
    """Read one `wgd_ratio_scores.csv` → (canonical frame, missing optional columns)."""
    df = pl.read_csv(
        path,
        null_values=["NA", "", "nan", "NaN"],
        infer_schema_length=10000,
        schema_overrides={c: t for c, t in REQUIRED.items()},
    )
    missing_req = [c for c in REQUIRED if c not in df.columns]
    if missing_req:
        raise ValueError(f"{path}: missing required columns {missing_req}")
    missing_opt = [c for c in OPTIONAL if c not in df.columns]
    df = df.with_columns(pl.lit(None, dtype=OPTIONAL[c]).alias(c) for c in missing_opt)
    df = df.select(pl.col(c).cast(t, strict=True) for c, t in SCHEMA.items())
    # pandas writes NaN as an empty field (read as null); ratio NaN → null.
    df = df.with_columns(pl.col(pl.Float64).fill_nan(None))

    ids = set(df["tumour_id"].unique())
    if ids != {tumour_dir}:
        raise ValueError(f"{path}: tumour_id {sorted(ids)} != directory {tumour_dir!r}")
    # Clone IDs are already canonical (`clone1`) in ALPACA output, as in the
    # `alpaca` table; kept verbatim.
    return df, missing_opt


def build() -> pl.DataFrame:
    files = _discover()
    log(f"{TABLE}: found {len(files)} tumour wgd_ratio_scores.csv files")

    frames: list[pl.DataFrame] = []
    for tumour_dir, path in files:
        df, missing = read_wgd_ratio_file(path, tumour_dir)
        if missing:
            log(f"{TABLE}: {tumour_dir} lacks {missing} (older ALPACA?) — filled with null")
        frames.append(df)

    df = (
        pl.concat(frames, how="vertical")
        if frames
        else pl.DataFrame(schema=SCHEMA)
    )
    df = df.with_columns((pl.col("parent") == "diploid").alias("is_mrca"))

    dup = df.filter(pl.struct("tumour_id", "clone").is_duplicated())
    if dup.height:
        raise ValueError(f"{TABLE}: duplicate (tumour_id, clone):\n{dup}")

    out = canonical_output_path(TABLE)
    df.write_parquet(out)
    log(f"wrote {out} ({df.height} rows × {df.width} cols)")
    return df


if __name__ == "__main__":
    build()
