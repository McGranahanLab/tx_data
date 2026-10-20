"""Build seeding_table → canonical Parquet.

Per-tumour metastatic seeding characterisation. See
``catalog/seeding_table.yml``.

Source format: space-separated text (R ``write.table`` default), header row.
``tumour_id`` flavour auto-detected via ``tx_data.ids.apply_convention_scalar``.
"""

from __future__ import annotations

import polars as pl

from tx_data.builds._base import canonical_output_path, log
from tx_data.ids import apply_convention_expr
from tx_data.paths import id_convention
from tx_data.sources import resolve_source

TABLE = "seeding_table"

_RENAMES = {
    "seedingClonality": "seeding_clonality",
    "seedingPhyletic": "seeding_phyletic",
    "individualSeedingClonality": "individual_seeding_clonality",
    "individualSeedingPhyletic": "individual_seeding_phyletic",
    "nRegions": "n_regions",
}


def build() -> pl.DataFrame:
    src = resolve_source(TABLE)
    log(f"read {src} (space-separated)")

    df = pl.read_csv(src, separator=" ", null_values=["NA", "", "nan"])

    missing = [c for c in _RENAMES.keys() if c not in df.columns]
    if missing:
        raise ValueError(
            f"seeding_table source is missing expected columns: {missing}. "
            f"Got columns: {df.columns!r}"
        )
    df = df.rename(_RENAMES)

    if df.schema["n_regions"] != pl.Int64:
        df = df.with_columns(pl.col("n_regions").cast(pl.Int64))

    # tumour_id → configured convention.
    df = df.with_columns(apply_convention_expr("tumour_id"))

    # Uniqueness: one row per tumour.
    dupe = df.group_by("tumour_id").len().filter(pl.col("len") > 1)
    if dupe.height:
        raise ValueError(
            f"seeding_table: {dupe.height} tumour_id values have >1 row.\n"
            f"First offenders: {dupe.head(5).to_dicts()}"
        )

    # Length invariant: n_regions must match split lengths in the two
    # individual_* columns.
    mismatched = df.with_columns(
        pl.col("individual_seeding_clonality").str.split(";").list.len().alias("_n_clon"),
        pl.col("individual_seeding_phyletic").str.split(";").list.len().alias("_n_phyl"),
    ).filter(
        (pl.col("_n_clon") != pl.col("n_regions"))
        | (pl.col("_n_phyl") != pl.col("n_regions"))
    )
    if mismatched.height:
        raise ValueError(
            f"seeding_table: {mismatched.height} rows have length mismatch between "
            f"n_regions and individual_seeding_* split length.\n"
            f"First offenders: "
            f"{mismatched.select(['tumour_id', 'n_regions', '_n_clon', '_n_phyl']).head(5).to_dicts()}"
        )

    df = df.select(
        [
            "tumour_id",
            "seeding_clonality",
            "seeding_phyletic",
            "individual_seeding_clonality",
            "individual_seeding_phyletic",
            "n_regions",
        ]
    ).sort("tumour_id")

    out = canonical_output_path(TABLE)
    df.write_parquet(out)
    log(f"wrote {out} ({df.height} rows × {df.width} cols) [id_convention={id_convention()}]")
    return df


if __name__ == "__main__":
    build()
