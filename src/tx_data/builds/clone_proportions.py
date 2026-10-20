"""Build clone_proportions → canonical Parquet.

Reads per-tumour `ALPACA/input/<tumour_id>/cp_table.csv` files. Source is
wide-format (N sample columns + trailing `clone` column); this builder
melts to long canonical form `(tumour_id, clone, sample_id, proportion)`
and concatenates every tumour into one parquet.

tumour_id is derived from the containing directory name — the source path
uses canonical hyphen-and-capital-T form (`LTX0001-Tumour1`), no coercion.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl

from tx_data.builds._base import canonical_output_path, log
from tx_data.sources import resolve_source

TABLE = "clone_proportions"


def _discover_cp_files() -> list[tuple[str, Path]]:
    """Return ``[(tumour_id, cp_table_path), ...]`` for every present tumour."""
    probe = resolve_source(TABLE, tumour_id="__probe__")
    input_root = probe.parent.parent          # `.../ALPACA/input/`
    filename = probe.name                     # `cp_table.csv`

    if not input_root.is_dir():
        log(f"clone_proportions input root does not exist: {input_root}")
        return []

    found: list[tuple[str, Path]] = []
    for sub in sorted(input_root.iterdir()):
        if not sub.is_dir():
            continue
        f = sub / filename
        if f.is_file():
            found.append((sub.name, f))
    return found


def _read_cp_table(tumour_id: str, path: Path) -> pl.DataFrame:
    """Read one wide cp_table.csv and melt to long form."""
    df = pl.read_csv(
        source=path,
        separator=",",
        null_values=["NA", ""],
        infer_schema_length=10000,
    )
    if "clone" not in df.columns:
        raise ValueError(f"{path}: expected trailing `clone` column, got {df.columns}")

    sample_cols = [c for c in df.columns if c != "clone"]
    long = df.unpivot(
        index=["clone"],
        on=sample_cols,
        variable_name="sample_id",
        value_name="proportion",
    )
    return long.with_columns(
        [
            pl.lit(tumour_id).alias("tumour_id"),
            pl.col("clone").cast(pl.String),
            pl.col("sample_id").cast(pl.String),
            pl.col("proportion").cast(pl.Float64),
        ]
    ).select(["tumour_id", "clone", "sample_id", "proportion"])


def build() -> pl.DataFrame:
    files = _discover_cp_files()
    log(f"clone_proportions: found {len(files)} tumour cp_table.csv files")

    frames: list[pl.DataFrame] = []
    for tumour_id, path in files:
        try:
            frames.append(_read_cp_table(tumour_id, path))
        except Exception as e:
            log(f"clone_proportions: skipping {tumour_id}: {e}")

    if not frames:
        df = pl.DataFrame(
            schema={
                "tumour_id": pl.String,
                "clone": pl.String,
                "sample_id": pl.String,
                "proportion": pl.Float64,
            }
        )
    else:
        df = pl.concat(frames, how="vertical_relaxed")

    out = canonical_output_path(TABLE)
    df.write_parquet(out)
    log(f"wrote {out} ({df.height} rows × {df.width} cols)")
    return df


if __name__ == "__main__":
    build()
