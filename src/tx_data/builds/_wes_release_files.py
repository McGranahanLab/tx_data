"""Shared reading for the per-tumour WES sources in the tx842 release tree (alias `tx842_release`).

Each table is one file per tumour under `<patient>/<folder>/tumour_<n>/` (`mutation_clusters` for
SCoutput, `mutation_trees` for treeTable / allTrees / clusterInfo; Piotr, 2026-11-24). Tumours are
discovered from the folders; `tumour_id` comes from the path (`LTX0001` + `tumour_1` → `LTX0001-Tumour1`).
CONIPHER writes SAMPLE R-style (`LTX0001_SU_T1.R1..<hash>`); `dash_sample` converts it to the form used
everywhere else in tx_data (`LTX0001_SU_T1-R1--<hash>`).
"""

from __future__ import annotations

import re
from collections.abc import Callable
from pathlib import Path

import polars as pl

from tx_data.builds._base import Issues, canonical_output_path, log
from tx_data.sources import resolve_source, source_entry

TUMOUR_DIR = re.compile(r"tumour_(\d+)$")
ID_COLS = {"tumour_id": pl.String, "patient_id": pl.String}


def tumour_id(patient: str, n: str) -> str:
    return f"{patient}-Tumour{n}"


def discover(table: str, issues: Issues | None = None) -> list[tuple[str, str, Path]]:
    """`[(patient_id, n, path), ...]` for every `<patient>/<folder>/tumour_<n>/` holding the table's file;
    tumour folders without it are recorded in `issues` (skipped)."""
    pattern = source_entry(table)["relative_path_pattern"]
    depth = len(Path(pattern).parts)                 # {patient_id}/<folder>/tumour_{n}/<file>
    folder = Path(pattern).parts[1]
    root = resolve_source(table, patient_id="__probe__", tumour_n="1").parents[depth - 1]
    if not root.is_dir():
        log(f"{table}: release root does not exist: {root}")
        return []
    found, empty = [], []
    for tdir in sorted(root.glob(f"*/{folder}/tumour_*")):
        m = TUMOUR_DIR.search(tdir.name)
        if not m or tdir.parts[-3].startswith("_"):
            continue
        patient, n = tdir.parts[-3], m.group(1)
        f = resolve_source(table, patient_id=patient, tumour_n=n)
        (found if f.is_file() else empty).append((patient, n, f))
    for p, n, f in empty:
        if issues is not None:
            issues.add(tumour_id(p, n), "skipped", f"no {f.name} in the tumour folder", f)
        else:
            log(f"warning: {table}: no {f.name} for {tumour_id(p, n)}")
    return found


def dash_sample(col: str = "sample") -> pl.Expr:
    """R-style `LTX0001_SU_T1.R1..<hash>` -> `LTX0001_SU_T1-R1--<hash>`; dash-form names unchanged."""
    return pl.col(col).str.replace_all(r"\.\.", "--").str.replace_all(r"\.", "-").alias(col)


def build_per_tumour(table: str, schema: dict, read_one: Callable[[str, str, Path], pl.DataFrame],
                     sort: tuple[str, ...]) -> pl.DataFrame:
    """Read every discovered tumour with `read_one(patient, n, path)` (frame without the id columns),
    attach tumour_id / patient_id, check the schema, write the canonical parquet. A tumour whose file is
    missing or can't be read is left out and recorded in data/_issues/<table>.tsv."""
    issues = Issues(table)
    files = discover(table, issues)
    frames = []
    for p, n, f in files:
        try:
            frames.append(read_one(p, n, f).with_columns(pl.lit(tumour_id(p, n)).alias("tumour_id"),
                                                         pl.lit(p).alias("patient_id")))
        except Exception as e:
            issues.add(tumour_id(p, n), "skipped", f"unreadable: {type(e).__name__}: {e}", f)
    log(f"{table}: read {len(frames)} WES tumours")
    if frames:
        df = pl.concat(frames, how="diagonal_relaxed").select([pl.col(c).cast(t) for c, t in schema.items()])
    else:
        df = pl.DataFrame(schema=schema)
    df = df.sort(*sort, nulls_last=True)
    out = canonical_output_path(table)
    df.write_parquet(out)
    log(f"wrote {out} ({df.height} rows × {df.width} cols)")
    issues.write()
    return df
