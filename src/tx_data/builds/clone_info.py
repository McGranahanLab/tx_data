"""Build clone_info → canonical Parquet.

Per (tumour × clone) classification. See ``catalog/clone_info.yml``.

Source ``tumour_id`` flavour is auto-detected — we've seen CRUK ``_N``,
CRUK ``-ClusterN``, LTX ``-TumourN``, LTX ``_tumourN`` and (as of the
2026-09-23 release) LTX ``_N``. ``tx_data.ids.apply_convention_scalar``
handles all of them and emits the form dictated by
``config/paths.yml:id_convention``. Only ONE ``tumour_id`` column is
emitted — no companion CRUK/LTX field. Downstream code that needs the
other form uses ``tx_data.ids`` (which reads the ``id_crosswalk`` table).

All source rows are preserved. ``is_tree_clone == False`` clusters are
real mutations that CONIPHER did not fit onto the tree; individual analyses
decide whether to include them (they must filter locally).
"""

from __future__ import annotations

import polars as pl

from tx_data.builds._base import canonical_output_path, log
from tx_data.ids import apply_convention_expr
from tx_data.paths import id_convention
from tx_data.sources import resolve_source

TABLE = "clone_info"

# Rename map: source-side → canonical.
_RENAMES = {
    "mutation_cluster": "clone",
    "treeClones": "is_tree_clone",
    "seedingClones": "is_seeding_clone",
    "clonalClones": "is_clonal",
    "sharedClones": "is_shared",
    "primaryClones": "is_primary_clone",
    "metClones": "is_met_clone",
}


def _canonical_clone(value) -> str:
    """Coerce a clone identifier (int / string / clone-prefixed) to `cloneN`."""
    if value is None:
        raise ValueError("clone_info: null clone id encountered")
    s = str(value).strip()
    if s.startswith("clone"):
        return s
    try:
        n = int(float(s))
    except ValueError as exc:
        raise ValueError(f"Un-parseable clone id: {value!r}") from exc
    return f"clone{n}"


def build() -> pl.DataFrame:
    src = resolve_source(TABLE)
    log(f"read {src} (tab-separated)")

    df = pl.read_csv(src, separator="\t", null_values=["NA", "", "nan"])

    missing = [c for c in _RENAMES.keys() if c not in df.columns]
    if missing:
        raise ValueError(
            f"clone_info source is missing expected columns: {missing}. "
            f"Got columns: {df.columns!r}"
        )
    df = df.rename(_RENAMES)

    # Boolean coercion — R writes TRUE / FALSE as strings.
    bool_cols = [
        "is_tree_clone",
        "is_seeding_clone",
        "is_clonal",
        "is_shared",
        "is_primary_clone",
        "is_met_clone",
    ]
    for col in bool_cols:
        if df.schema[col] == pl.Boolean:
            continue
        df = df.with_columns(
            pl.when(
                pl.col(col).cast(pl.String).str.to_uppercase().is_in(["TRUE", "T", "1"])
            )
            .then(True)
            .when(
                pl.col(col).cast(pl.String).str.to_uppercase().is_in(["FALSE", "F", "0"])
            )
            .then(False)
            .otherwise(None)
            .alias(col)
        )

    df = df.with_columns(
        pl.col("clone").map_elements(_canonical_clone, return_dtype=pl.String)
    )

    convention = id_convention()
    df = df.with_columns(apply_convention_expr("tumour_id"))

    dupe_count = df.group_by(["tumour_id", "clone"]).len().filter(pl.col("len") > 1)
    if dupe_count.height:
        raise ValueError(
            f"clone_info: {dupe_count.height} (tumour_id, clone) pairs have >1 row.\n"
            f"First offenders: {dupe_count.head(5).to_dicts()}"
        )

    df = df.select(
        [
            "tumour_id",
            "clone",
            "is_tree_clone",
            "is_seeding_clone",
            "is_clonal",
            "is_shared",
            "is_primary_clone",
            "is_met_clone",
        ]
    ).sort(["tumour_id", "clone"])

    out = canonical_output_path(TABLE)
    df.write_parquet(out)
    log(
        f"wrote {out} ({df.height} rows × {df.width} cols) "
        f"[id_convention={convention}, non_tree_clones={df.filter(~pl.col('is_tree_clone')).height}]"
    )
    return df


def tumours_by_cohort(cohort: str) -> list[str]:
    """Return the tumour_ids that belong to a given analysis cohort.

    - ``primary`` — any tumour with at least one clone flagged as
      ``is_primary_clone`` (i.e. the primary tumour was sampled and analysed).
    - ``mets``    — any tumour with at least one clone flagged as
      ``is_met_clone`` (i.e. a metastatic sample was analysed).

    Reads the already-built ``clone_info`` parquet; raises if missing.
    """
    if cohort not in ("primary", "mets"):
        raise ValueError(f"cohort must be 'primary' or 'mets', got {cohort!r}")
    p = canonical_output_path(TABLE)
    if not p.is_file():
        raise FileNotFoundError(
            f"clone_info parquet missing: {p}. Run `pixi run build` first."
        )
    df = pl.read_parquet(p)
    flag_col = "is_primary_clone" if cohort == "primary" else "is_met_clone"
    return (
        df.filter(pl.col(flag_col) == True)
          .select("tumour_id")
          .unique()
          .sort("tumour_id")
          .to_series()
          .to_list()
    )


if __name__ == "__main__":
    build()
