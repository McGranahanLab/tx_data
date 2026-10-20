"""Build conipher_trees → canonical Parquet(s).

Reads the aggregated Tx842 ``tx842_trees.RDS``, extracts the two tables the
ALPACA figure pipeline consumes (CCF per clone/region, corrected clonality
per clone/region), and emits:

    data/conipher_ccf.parquet         # (tumour_id, clone, region, ccf)
    data/conipher_clonality.parquet   # (tumour_id, clone, region, clonality)

Both are long format, one row per (tumour × clone × region). CCF values are
kept as integer "percentage" (0..~105) to match the source — the paper's
R scripts use raw min/max/mean on these values and we don't want to make
rescaling an invisible side-effect.

Rationale for shelling out to R
-------------------------------
``tx_data``'s pixi env is Python-only by design. The aggregated RDS is an
unnamed list of 15-slot per-patient bundles — too deeply nested for
pyreadr's single-object import to handle cleanly. We delegate the RDS → TSV
dump to a tiny R helper (``scripts/_dump_conipher_trees.R``) and keep all
canonicalisation on the Python side.

If ``Rscript`` is not on PATH we fail with a clear message pointing at the
``tx800_alpaca_replica`` pixi env, which already carries ``r-base``.
"""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path

import polars as pl

from tx_data.builds._base import canonical_output_path, log
from tx_data.ids import apply_convention_expr
from tx_data.paths import REPO_ROOT, id_convention
from tx_data.sources import resolve_source

TABLE_CCF        = "conipher_ccf"
TABLE_CLONALITY  = "conipher_clonality"

_DUMPER_R = REPO_ROOT / "scripts" / "_dump_conipher_trees.R"
_TX800_LOCAL_RSCRIPT = (
    REPO_ROOT.parent / "tx800_alpaca_replica" / ".pixi" / "envs" / "default" / "bin" / "Rscript"
)


def _resolve_rscript() -> str:
    rscript = shutil.which("Rscript")
    if rscript:
        return rscript
    if _TX800_LOCAL_RSCRIPT.is_file():
        return str(_TX800_LOCAL_RSCRIPT)
    raise FileNotFoundError(
        "Rscript not on PATH and tx800_alpaca_replica's pixi R isn't available either. "
        "Install an R env with base R (`r-base`), or run this build from a shell that "
        "already has Rscript."
    )


def _run_dumper(rds_path: Path, out_dir: Path) -> tuple[Path, Path]:
    rscript = _resolve_rscript()
    cmd = [
        rscript, str(_DUMPER_R),
        f"--input={rds_path}",
        f"--output-dir={out_dir}",
    ]
    log(f"conipher_trees: invoking {' '.join(cmd)}")
    res = subprocess.run(cmd, check=False, capture_output=True, text=True)
    if res.stdout:
        for line in res.stdout.rstrip().splitlines():
            log(line)
    if res.returncode != 0:
        raise RuntimeError(
            f"conipher_trees dumper failed (exit {res.returncode}):\n"
            f"stdout:\n{res.stdout}\nstderr:\n{res.stderr}"
        )
    return out_dir / "conipher_ccf.tsv", out_dir / "conipher_clonality.tsv"


def _canonicalise(df: pl.DataFrame, value_col: str, value_dtype: pl.DataType) -> pl.DataFrame:
    """Canonicalise tumour_id + clone; rename `source_*` columns to canonical form."""
    out = df.rename({"source_tumour_id": "tumour_id", "source_clone": "clone"})
    out = out.with_columns(apply_convention_expr("tumour_id"))
    # clone: "1","2",... → "clone1","clone2",...
    out = out.with_columns(
        pl.when(pl.col("clone").cast(pl.String).str.starts_with("clone"))
        .then(pl.col("clone").cast(pl.String))
        .otherwise(pl.lit("clone") + pl.col("clone").cast(pl.String))
        .alias("clone")
    )
    out = out.with_columns(pl.col("value").cast(value_dtype).alias(value_col)).drop("value")
    # Primary key uniqueness guard.
    dupes = out.group_by(["tumour_id", "clone", "region"]).len().filter(pl.col("len") > 1)
    if dupes.height:
        raise ValueError(
            f"{value_col}: {dupes.height} duplicate (tumour_id, clone, region) rows.\n"
            f"First offenders: {dupes.head(5).to_dicts()}"
        )
    return out.select(["tumour_id", "clone", "region", value_col]).sort(
        ["tumour_id", "clone", "region"]
    )


def build() -> tuple[pl.DataFrame, pl.DataFrame]:
    rds_path = resolve_source("conipher_trees")
    log(f"conipher_trees: reading {rds_path}")

    with tempfile.TemporaryDirectory(prefix="tx_data_conipher_") as tmpdir:
        ccf_tsv, cln_tsv = _run_dumper(rds_path, Path(tmpdir))

        ccf_raw = pl.read_csv(ccf_tsv, separator="\t")
        cln_raw = pl.read_csv(cln_tsv, separator="\t")

    ccf = _canonicalise(ccf_raw, value_col="ccf", value_dtype=pl.Int32)
    cln = _canonicalise(cln_raw, value_col="clonality", value_dtype=pl.String)

    # Allowed-values guard for clonality.
    allowed = {"clonal", "subclonal", "absent"}
    seen = set(cln["clonality"].unique().to_list()) - {None}
    unknown = seen - allowed
    if unknown:
        raise ValueError(
            f"conipher_clonality: unexpected values {sorted(unknown)} (allowed: {sorted(allowed)})"
        )

    ccf_out = canonical_output_path(TABLE_CCF)
    cln_out = canonical_output_path(TABLE_CLONALITY)
    ccf.write_parquet(ccf_out)
    cln.write_parquet(cln_out)
    convention = id_convention()
    log(
        f"wrote {ccf_out} ({ccf.height} rows × {ccf.width} cols) "
        f"[id_convention={convention}, tumours={ccf.n_unique('tumour_id')}, "
        f"clonality_values={sorted(seen)}]"
    )
    log(f"wrote {cln_out} ({cln.height} rows × {cln.width} cols)")
    return ccf, cln


if __name__ == "__main__":
    build()
