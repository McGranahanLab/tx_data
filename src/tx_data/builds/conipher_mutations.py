"""Build conipher_mutations (WES) → canonical Parquet.

CONIPHER mutations of each WES tumour from the tx842 release tree (Piotr, 2026-11-24):
`<patient>/mutation_clusters/tumour_<n>/<patient>_<n>.SCoutput.CLEAN.tsv` (one row per mutation × sample)
joined with `<patient>/mutation_trees/tumour_<n>/treeTable.tsv` on (mutation_id, SAMPLE) for the tree cluster.
WES counterpart of `wgs_conipher_mutations`, same columns except the WGS id columns; input of the WGD
pipeline's step 00 (PGDD input) through WGD_pipeline's `config.wes.toml`.

Tumours are discovered from the `mutation_clusters/tumour_<n>/` folders. `tumour_id` comes from the path
(`LTX0001` + `tumour_1` → `LTX0001-Tumour1`); the file's CASE_ID must agree (`LTX0001_1`).
`tree_cluster` is null for mutations treeTable leaves out, and for a whole tumour without a treeTable
(logged). `sample_id` is the clone_proportions sample (`LTX0001_SU_T1-R1--<hash>`) matched ignoring
punctuation, as SCoutput may write it R-style (`LTX0001_SU_T1.R1..<hash>`); null when unmatched or when
clone_proportions isn't built.
"""

from __future__ import annotations

import re
from pathlib import Path

import polars as pl

from tx_data.builds._base import canonical_output_path, log
from tx_data.builds._wgs_tumour_files import read_tsv, require
from tx_data.normalize import canonical_chr_expr
from tx_data.sources import resolve_source

TABLE = "conipher_mutations"
TREE_TABLE_PATTERN = "relative_path_pattern_tree_table"

SOURCE = {  # SCoutput column -> canonical (name, type)
    "SAMPLE": ("sample", pl.String), "CHR": ("chr", pl.String), "POS": ("pos", pl.Int64),
    "REF": ("ref", pl.String), "ALT": ("alt", pl.String), "mutation_id": ("mutation_id", pl.String),
    "CLUSTER": ("cluster", pl.Int64), "MUT_COPY": ("mut_copy", pl.Float64),
    "COPY_NUMBER_A": ("cn_a", pl.Float64), "COPY_NUMBER_B": ("cn_b", pl.Float64),
}
SCHEMA = {"tumour_id": pl.String, "patient_id": pl.String, "sample": pl.String, "sample_id": pl.String,
          "chr": pl.String, "pos": pl.Int64, "ref": pl.String, "alt": pl.String, "mutation_id": pl.String,
          "cluster": pl.Int64, "tree_cluster": pl.Int64, "mut_copy": pl.Float64, "cn_a": pl.Float64,
          "cn_b": pl.Float64, "ploidy": pl.Float64}
TUMOUR_DIR = re.compile(r"tumour_(\d+)$")


def discover() -> list[tuple[str, str, Path]]:
    """`[(patient_id, n, scoutput_path), ...]` for every `<patient>/mutation_clusters/tumour_<n>/` with a file."""
    probe = resolve_source(TABLE, patient_id="__probe__", tumour_n="1")
    root = probe.parents[3]                        # release_tx842/
    if not root.is_dir():
        log(f"{TABLE}: release root does not exist: {root}")
        return []
    found, empty = [], []
    for tdir in sorted(root.glob("*/mutation_clusters/tumour_*")):
        m = TUMOUR_DIR.search(tdir.name)
        if not m or tdir.parts[-3].startswith("_"):
            continue
        patient, n = tdir.parts[-3], m.group(1)
        f = resolve_source(TABLE, patient_id=patient, tumour_n=n)
        (found if f.is_file() else empty).append((patient, n, f))
    if empty:
        log(f"warning: {TABLE}: no SCoutput.CLEAN.tsv in {len(empty)} tumour folder(s): "
            f"{[f'{p}_{n}' for p, n, _ in empty]}")
    return found


def _norm(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9]", "", s)


def _sample_ids() -> dict[str, dict[str, str]]:
    """tumour_id -> {punctuation-free sample: clone_proportions sample_id}; empty if not built."""
    cp = canonical_output_path("clone_proportions")
    if not cp.is_file():
        log(f"warning: {TABLE}: clone_proportions not built, sample_id left null")
        return {}
    out: dict[str, dict[str, str]] = {}
    for t, s in pl.read_parquet(cp, columns=["tumour_id", "sample_id"]).unique().iter_rows():
        out.setdefault(t, {})[_norm(s)] = s
    return out


def read_one(patient: str, n: str, path: Path, sample_ids: dict[str, str]) -> pl.DataFrame:
    muts = read_tsv(path)
    require(muts, SOURCE, path)
    case = f"{patient}_{n}"
    if "CASE_ID" in muts.columns:
        bad = set(muts["CASE_ID"].drop_nulls()) - {case}
        if bad:
            raise ValueError(f"{path}: CASE_ID {sorted(bad)} does not match the folder ({case})")
    ploidy = pl.col("PLOIDY").cast(pl.Float64) if "PLOIDY" in muts.columns else pl.lit(None, pl.Float64)
    df = muts.select(*[pl.col(s).cast(t).alias(n_) for s, (n_, t) in SOURCE.items()], ploidy.alias("ploidy"))
    tt_path = resolve_source(TABLE, pattern_key=TREE_TABLE_PATTERN, patient_id=patient, tumour_n=n)
    if tt_path.is_file():
        tt = read_tsv(tt_path)
        require(tt, ("mutation_id", "SAMPLE", "treeCLUSTER"), tt_path)
        tt = tt.select(pl.col("mutation_id"), pl.col("SAMPLE").alias("sample"),
                       pl.col("treeCLUSTER").cast(pl.Int64).alias("tree_cluster")).unique(["mutation_id", "sample"])
        df = df.join(tt, on=["mutation_id", "sample"], how="left")
    else:
        df = df.with_columns(pl.lit(None, pl.Int64).alias("tree_cluster"))
    return df.with_columns(
        canonical_chr_expr("chr"),
        pl.lit(f"{patient}-Tumour{n}").alias("tumour_id"),
        pl.lit(patient).alias("patient_id"),
        pl.col("sample").map_elements(lambda s: sample_ids.get(_norm(s)), return_dtype=pl.String).alias("sample_id"),
    )


def build() -> pl.DataFrame:
    files = discover()
    ids = _sample_ids()
    frames, no_tree = [], []
    for patient, n, path in files:
        if not resolve_source(TABLE, pattern_key=TREE_TABLE_PATTERN, patient_id=patient, tumour_n=n).is_file():
            no_tree.append(f"{patient}_{n}")
        frames.append(read_one(patient, n, path, ids.get(f"{patient}-Tumour{n}", {})))
    log(f"{TABLE}: read {len(frames)} WES tumours")
    if no_tree:
        log(f"warning: {TABLE}: no treeTable.tsv for {len(no_tree)} tumour(s), tree_cluster null: {no_tree}")
    if frames:
        df = pl.concat(frames, how="diagonal_relaxed").select([pl.col(c).cast(t) for c, t in SCHEMA.items()])
    else:
        df = pl.DataFrame(schema=SCHEMA)
    if ids and df.height:
        unmatched = df.filter(pl.col("sample_id").is_null())["sample"].unique().to_list()
        if unmatched:
            log(f"warning: {TABLE}: {len(unmatched)} sample(s) not in clone_proportions: {unmatched[:10]}")
    df = df.sort("tumour_id", "sample", "chr", "pos", nulls_last=True)
    out = canonical_output_path(TABLE)
    df.write_parquet(out)
    log(f"wrote {out} ({df.height} rows × {df.width} cols)")
    return df


if __name__ == "__main__":
    build()
