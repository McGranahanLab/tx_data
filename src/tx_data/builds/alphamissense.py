"""Build alphamissense → canonical Parquet.

Reads per-tumour AlphaMissense-annotated muttable files. Confirmed by
User 2026-09-02 that the real upstream is NOT raw VEP output —
it's a muttable-shaped TSV with `alpha_missense` (score) and `am_class`
columns appended. That means:

- No `Location` column to split; `chr` and `pos` are direct.
- No `Allele` column; the alt allele is muttable's `var`.
- No `Consequence` / `Gene` columns (they weren't used downstream anyway).
- Files are per-mutation × per-sample (long form, like muttable), so
  we still dedup to per-(tumour, chr, pos, alt) — the same variant
  appears once per sample it was called in, all rows carrying the
  identical AM score.

File layout:
  <root>/output/annotated_muttables/tx842/<tumour_id>_muttable_alpha.tsv
where `<tumour_id>` in the filename is the muttable spelling
(`LTX0001_tumour1`), which the loader normalises back to canonical form.
"""

from __future__ import annotations

import re
from pathlib import Path

import polars as pl

from tx_data.builds._base import canonical_output_path, log
from tx_data.normalize import canonical_chr_expr, canonical_tumour_id
from tx_data.sources import resolve_source, source_entry

TABLE = "alphamissense"

_NULL_MARKERS = ["NA", "", "-", "nan"]

# Real upstream emits DeepMind's original labels ("Likely benign" /
# "Likely pathogenic"), the VEP plugin (which our downstream code assumes)
# renames them to "benign" / "pathogenic". `not_classified` means AM had no
# opinion for this variant — coerced to null so downstream can distinguish
# "AM said benign" from "AM had no answer".
_AM_CLASS_CANONICAL = {
    "Likely benign": "benign",
    "Likely pathogenic": "pathogenic",
    "ambiguous": "ambiguous",
    "benign": "benign",
    "pathogenic": "pathogenic",
    "not_classified": None,
}


def _read_am(p: Path) -> pl.DataFrame:
    """Read one AlphaMissense-annotated muttable TSV.

    Only pulls the columns we actually need. `chr` is forced to String
    so X / Y / MT don't crash the auto-inferrer; `alpha_missense` is
    forced to Float64 so a chunk of all-null rows early in the file
    doesn't trap us into an Int/Utf8 inference.
    """
    return pl.read_csv(
        source=p,
        separator="\t",
        null_values=_NULL_MARKERS,
        infer_schema_length=10000,
        schema_overrides={
            "chr": pl.String,
            "pos": pl.Int64,
            "var": pl.String,
            "alpha_missense": pl.Float64,
            "am_class": pl.String,
        },
        columns=["chr", "pos", "var", "alpha_missense", "am_class"],
    )


def _discover_am_files() -> list[tuple[str, Path]]:
    """Return ``[(tumour_id_from_filename, path), …]`` for every AM file on disk.

    Uses ``relative_path_pattern`` from ``_sources.yml`` to derive both the
    directory to scan and the filename shape to parse — so a future rename
    only needs the yml change.
    """
    entry = source_entry(TABLE)
    pattern = entry["relative_path_pattern"]
    probe = resolve_source(TABLE, tumour_id="__probe__")
    scan_dir = probe.parent
    filename_template = Path(pattern).name
    prefix, suffix = filename_template.split("{tumour_id}", 1)
    fname_re = re.compile(re.escape(prefix) + r"(?P<tid>.+?)" + re.escape(suffix) + r"$")

    if not scan_dir.is_dir():
        log(f"AM scan dir does not exist: {scan_dir}")
        return []

    found: list[tuple[str, Path]] = []
    for p in sorted(scan_dir.iterdir()):
        if not p.is_file():
            continue
        m = fname_re.match(p.name)
        if not m:
            continue
        found.append((m.group("tid"), p))
    return found


def build() -> pl.DataFrame:
    discovered = _discover_am_files()
    log(f"discovered {len(discovered)} AM files")

    frames = []
    for tid_from_filename, p in discovered:
        tid_canonical = canonical_tumour_id(tid_from_filename)
        df = _read_am(p).with_columns(pl.lit(tid_canonical).alias("tumour_id"))
        frames.append(df)

    raw = pl.concat(frames, how="vertical_relaxed") if frames else pl.DataFrame()

    if raw.is_empty():
        log("no AM rows to write")
        out = canonical_output_path(TABLE)
        raw.write_parquet(out)
        return raw

    # Rename to the canonical column names downstream expects, coerce
    # `am_class` from whatever variant the upstream emitted to the canonical
    # {benign, ambiguous, pathogenic, null} set, and normalise `chr` (strip
    # any `chr` prefix so it matches muttable / driver_list / alpaca).
    parsed = raw.rename({"alpha_missense": "am_pathogenicity", "var": "alt"}).with_columns(
        [
            pl.col("am_class").replace(_AM_CLASS_CANONICAL, default=None),
            canonical_chr_expr("chr"),
        ]
    )

    # Dedup to per-variant. Real files are long-form (per mutation × sample),
    # so the same (tumour, chr, pos, alt) may appear many times with the
    # identical AM score; `max` collapses agreeing values and ignores nulls.
    scored = parsed.group_by(["tumour_id", "chr", "pos", "alt"]).agg(
        [
            pl.col("am_pathogenicity").max().alias("am_pathogenicity"),
            pl.col("am_class").drop_nulls().first().alias("am_class"),
        ]
    )

    out = canonical_output_path(TABLE)
    scored.write_parquet(out)
    log(f"wrote {out} ({scored.height} rows × {scored.width} cols)")
    return scored


if __name__ == "__main__":
    build()
