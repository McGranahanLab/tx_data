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

from dataclasses import dataclass
from functools import lru_cache

import polars as pl

from tx_data.sources import resolve_source

SOURCE_NAME = "sample_hash_map"


@dataclass(frozen=True)
class SampleRecord:
    patient_id: str
    patient_tumour: str
    lesion: str
    histology_group: str


def extract_hash_suffix(sample_name_hash: str) -> str:
    """Return the suffix after `--` in a full `sample_name_hash` string.

    The only guaranteed structure of a sample identifier is
    `<LTX|CRUK><digits>_..._--<hash>`. Everything before `--` is
    presentational and varies between releases; the hash suffix is the
    stable join key against the mapping table.
    """
    if "--" not in sample_name_hash:
        raise ValueError(
            f"sample_name_hash missing '--' hash delimiter: {sample_name_hash!r}"
        )
    return sample_name_hash.rsplit("--", 1)[-1]


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
