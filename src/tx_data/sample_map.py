"""Sample-hash → tumour/lesion lookup backed by an explicit mapping TSV.

The release mapping file is the authoritative crosswalk between an opaque
`sample_hash` (the suffix after `--` in a full `sample_name_hash` string)
and its canonical identifiers:

    sample_hash   patient_id   patient_tumour     lesion    histology_group
    9f1a26d3h2l3  LTX0287      LTX0287_tumour1    lesion1   LUSC

Build scripts that would otherwise regex-parse name tokens should route
every lookup through here so that any new release-format quirk is a
one-line fix in the TSV, not a brittle rebuild of a pattern.

Normal samples (buffy coats etc.) are intentionally **absent** from this
mapping — they don't correspond to a tumour. Builds that need to flag
normals should do so independently (e.g. the `_N<n>` substring check in
`kallisto.py`).
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from functools import lru_cache

import polars as pl

from tx_data.sources import resolve_source

SOURCE_NAME = "sample_hash_map"

# Stem = everything before `--` in a `sample_name_hash`. For tumour samples
# that follow the `<LTX|CRUK><digits>_<TYPE>_T<n>-R<n>` convention (reliable
# for kallisto RNA data per User 2026-10-21), the stem decomposes further.
_TUMOUR_STEM_RE = re.compile(
    r"^(?P<patient_id>(?:LTX|CRUK)\d+)"
    r"_[A-Z]+"
    r"_T(?P<tumour_ordinal>\d+)"
    r"-R(?P<region_ordinal>\d+)$"
)


@dataclass(frozen=True)
class SampleRecord:
    patient_id: str
    patient_tumour: str
    lesion: str
    histology_group: str


def extract_hash_suffix(sample_name_hash: str) -> str:
    """Return the suffix after `--` in a full `sample_name_hash` string.

    The hash suffix is the DNA join key against `sample_hash_map`.
    """
    if "--" not in sample_name_hash:
        raise ValueError(
            f"sample_name_hash missing '--' hash delimiter: {sample_name_hash!r}"
        )
    return sample_name_hash.rsplit("--", 1)[-1]


def extract_sample_stem(sample_name_hash: str) -> str:
    """Return the stem (everything before `--`) of a `sample_name_hash`.

    DNA and RNA assays of the same biopsy share the same stem but have
    different hashes, so the stem is the stable RNA↔DNA join key.
    """
    if "--" not in sample_name_hash:
        raise ValueError(
            f"sample_name_hash missing '--' hash delimiter: {sample_name_hash!r}"
        )
    return sample_name_hash.rsplit("--", 1)[0]


@dataclass(frozen=True)
class StemParts:
    patient_id: str
    tumour_ordinal: int
    region_ordinal: int


def parse_tumour_stem(stem: str) -> StemParts | None:
    """Parse a tumour-sample stem into `(patient_id, tumour_ordinal, region_ordinal)`.

    Returns `None` for stems that don't match the tumour pattern — callers
    should combine with their own normal-detection logic (e.g. `_N<n>`).
    """
    m = _TUMOUR_STEM_RE.match(stem)
    if m is None:
        return None
    return StemParts(
        patient_id=m.group("patient_id"),
        tumour_ordinal=int(m.group("tumour_ordinal")),
        region_ordinal=int(m.group("region_ordinal")),
    )


@lru_cache(maxsize=1)
def _load_map() -> dict[str, SampleRecord]:
    src = resolve_source(SOURCE_NAME)
    df = pl.read_csv(src, separator="\t", null_values=["NA", "", "nan"])
    expected = {
        "sample_hash", "patient_id", "patient_tumour", "lesion", "histology_group",
    }
    missing = expected - set(df.columns)
    if missing:
        raise ValueError(
            f"{src}: missing expected columns {sorted(missing)}; got {df.columns!r}"
        )

    dupes = df.group_by("sample_hash").len().filter(pl.col("len") > 1)
    if dupes.height:
        raise ValueError(
            f"{src}: {dupes.height} duplicate sample_hash row(s). "
            f"First offenders: {dupes.head(5).to_dicts()}"
        )

    out: dict[str, SampleRecord] = {}
    for row in df.iter_rows(named=True):
        out[row["sample_hash"]] = SampleRecord(
            patient_id=row["patient_id"],
            patient_tumour=row["patient_tumour"],
            lesion=row["lesion"],
            histology_group=row["histology_group"],
        )
    return out


@lru_cache(maxsize=1)
def _tumour_histology_index() -> dict[str, str]:
    """`patient_tumour → histology_group` index, used for samples not in the
    mapping by hash (e.g. RNA/kallisto samples, which share a `patient_tumour`
    with a DNA sample but carry a different hash)."""
    out: dict[str, str] = {}
    for rec in _load_map().values():
        prev = out.get(rec.patient_tumour)
        if prev is None:
            out[rec.patient_tumour] = rec.histology_group
        elif prev != rec.histology_group:
            raise ValueError(
                f"sample_hash_map: tumour {rec.patient_tumour!r} carries "
                f"conflicting histology_group values: {prev!r} vs "
                f"{rec.histology_group!r}"
            )
    return out


def lookup_tumour_histology(patient_tumour: str) -> str | None:
    """Return the histology group recorded against a `patient_tumour`.

    `patient_tumour` must be in the mapping's spelling (`LTX0001_tumour1`),
    not canonical form — callers should pass the pre-normalised value.
    """
    return _tumour_histology_index().get(patient_tumour)


def lookup(sample_name_hash: str) -> SampleRecord | None:
    """Return the mapping record for a full `sample_name_hash`, or `None`.

    `None` means the hash suffix isn't in the mapping — typically a
    non-tumour sample (normal, buffy coat) that the release maps elsewhere.
    """
    return _load_map().get(extract_hash_suffix(sample_name_hash))


def lookup_strict(sample_name_hash: str) -> SampleRecord:
    """Like `lookup`, but raises if the hash is absent."""
    rec = lookup(sample_name_hash)
    if rec is None:
        raise KeyError(
            f"sample_name_hash not in {SOURCE_NAME}: {sample_name_hash!r}"
        )
    return rec


def reset_cache() -> None:
    """Clear the cached mapping — call between tests that swap the fixture."""
    _load_map.cache_clear()
    _tumour_histology_index.cache_clear()
