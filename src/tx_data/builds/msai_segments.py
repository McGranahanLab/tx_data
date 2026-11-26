"""Build msai_segments (WES) → canonical Parquet.

Fractional allele-specific CN per segment per sample for the WES cohort, from the ALPACA input
`ALPACA/input/<tumour_id>/ALPACA_input_table.csv` (Piotr, 2026-11-24). WES counterpart of
`wgs_consensus_segments` and input of WGD step 00 (`segments` in WGD_pipeline's `config.wes.toml`).

Why the ALPACA input: WES copy number is ASCAT, then Refphase, and Refphase detects MSAI itself, so
its output is already MSAI-corrected and there are no separate MSAI_adjusted files (WGS runs Battenberg,
which doesn't, hence the post-hoc MSAI correction there). The ALPACA input carries the same Refphase
solution that feeds CONIPHER. The name `msai_segments` is kept for symmetry with the WGS input.

`ALPACA/input` and `ALPACA/output` are symbolic links that the ALPACA pipeline does not update, so they
can point at different ALPACA versions (Piotr, 2026-11-24). The build logs where both resolve to.

`segment` (`<chr>_<start>_<end>`) is split into chr / start / end. cpnA / cpnB are phased alleles,
kept as cn_a / cn_b; cn_major / cn_minor are their max / min (what step 00 reads).
"""

from __future__ import annotations

import os
from pathlib import Path

import polars as pl

from tx_data.builds._base import Issues, canonical_output_path, log
from tx_data.builds._wes_release_files import dash_sample
from tx_data.builds._wgs_tumour_files import r_int
from tx_data.normalize import canonical_chr_expr
from tx_data.sources import resolve_source

TABLE = "msai_segments"
SOURCE_COLS = ("segment", "sample", "cpnA", "cpnB")
SCHEMA = {"tumour_id": pl.String, "patient_id": pl.String, "sample": pl.String, "chr": pl.String,
          "start": pl.Int64, "end": pl.Int64, "cn_a": pl.Float64, "cn_b": pl.Float64, "cn_major": pl.Float64,
          "cn_minor": pl.Float64, "cn_total": pl.Float64}


def log_alpaca_links(input_root: Path) -> None:
    """Log where ALPACA/input and ALPACA/output resolve to (unsynchronised symlinks)."""
    for d in (input_root, input_root.parent / "output"):
        kind = f"symlink -> {os.readlink(d)}" if d.is_symlink() else "directory"
        log(f"{TABLE}: {d} is a {kind} (resolves to {d.resolve()})")


def discover(issues: Issues | None = None) -> list[tuple[str, Path]]:
    probe = resolve_source(TABLE, tumour_id="__probe__")
    root = probe.parent.parent                          # ALPACA/input/
    if not root.is_dir():
        log(f"{TABLE}: ALPACA input root does not exist: {root}")
        return []
    log_alpaca_links(root)
    found = [(d.name, d / probe.name) for d in sorted(root.iterdir()) if d.is_dir() and (d / probe.name).is_file()]
    missing = [d.name for d in sorted(root.iterdir()) if d.is_dir() and not (d / probe.name).is_file()]
    for t in missing:
        if issues is not None:
            issues.add(t, "skipped", f"ALPACA input folder has no {probe.name}", root / t / probe.name)
        else:
            log(f"warning: {TABLE}: no {probe.name} for {t}")
    return found


def read_one(tumour_id: str, path: Path) -> pl.DataFrame:
    df = pl.read_csv(path, infer_schema_length=0, null_values=["NA", "NaN", ""])
    missing = [c for c in SOURCE_COLS if c not in df.columns]
    if missing:
        raise ValueError(f"{path}: missing required columns {missing} (has {df.columns})")
    if "tumour_id" in df.columns:
        bad = set(df["tumour_id"].drop_nulls()) - {tumour_id}
        if bad:
            raise ValueError(f"{path}: tumour_id {sorted(bad)} does not match the folder ({tumour_id})")
    df = df.with_columns(pl.col("segment").str.split("_").list.to_struct(fields=["_chr", "_start", "_end"])
                           .alias("_seg")).unnest("_seg")
    df = r_int(df, ("_start", "_end"), path)
    a, b = pl.col("cpnA").cast(pl.Float64), pl.col("cpnB").cast(pl.Float64)
    return df.select(
        pl.lit(tumour_id).alias("tumour_id"), pl.lit(tumour_id.split("-")[0]).alias("patient_id"),
        pl.col("sample"), pl.col("_chr").alias("chr"), pl.col("_start").alias("start"),
        pl.col("_end").alias("end"), a.alias("cn_a"), b.alias("cn_b"),
        pl.max_horizontal(a, b).alias("cn_major"), pl.min_horizontal(a, b).alias("cn_minor"),
        (a + b).alias("cn_total"),
    ).with_columns(dash_sample(), canonical_chr_expr("chr"))


def build() -> pl.DataFrame:
    issues = Issues(TABLE)
    files = discover(issues)
    frames = []
    for t, f in files:
        try:
            frames.append(read_one(t, f))
        except Exception as e:
            issues.add(t, "skipped", f"unreadable: {type(e).__name__}: {e}", f)
    log(f"{TABLE}: read {len(frames)} WES tumours")
    df = (pl.concat(frames, how="diagonal_relaxed").select([pl.col(c).cast(t) for c, t in SCHEMA.items()])
          if frames else pl.DataFrame(schema=SCHEMA))
    df = df.sort("tumour_id", "sample", "chr", "start", nulls_last=True)
    out = canonical_output_path(TABLE)
    df.write_parquet(out)
    log(f"wrote {out} ({df.height} rows × {df.width} cols)")
    issues.write()
    return df


if __name__ == "__main__":
    build()
