"""Canonical-form normalisers for cross-file identifiers.

- ``tumour_id`` canonical form: ``LTX0001-Tumour1``  (hyphen, capital T).
- ``clone`` canonical form: ``clone1`` (string with 'clone' prefix; NaN → None).
- ``chr`` canonical form: no ``chr`` prefix (``"17"``, ``"X"``, ``"MT"``),
  always String. Applies to muttable.chr, alphamissense.chr,
  alpaca.segment_chr, driver_list.chr_hg19/hg38.

Rationale documented in ``catalog/wgd_calls.yml`` under ``id_normalisation``.
"""

from __future__ import annotations

import math
import re
from typing import Iterable

import polars as pl

_MUTTABLE_TUMOUR_RE = re.compile(r"^(LTX\d+)_tumour(\d+)$", re.IGNORECASE)
_WGD_TUMOUR_RE = re.compile(r"^(LTX\d+)-Tumour(\d+)$")


def canonical_tumour_id(value: str | None) -> str | None:
    """Coerce any known tumour_id spelling into canonical `LTX0001-Tumour1`."""
    if value is None:
        return None
    if _WGD_TUMOUR_RE.match(value):
        return value
    m = _MUTTABLE_TUMOUR_RE.match(value)
    if m:
        return f"{m.group(1)}-Tumour{int(m.group(2))}"
    raise ValueError(f"Unrecognised tumour_id form: {value!r}")


def canonical_tumour_id_expr(col: str) -> pl.Expr:
    """Polars expression version — vectorised coercion."""
    return (
        pl.when(pl.col(col).str.contains(r"^LTX\d+-Tumour\d+$"))
        .then(pl.col(col))
        .otherwise(
            pl.col(col).str.replace(r"^(LTX\d+)_tumour(\d+)$", r"${1}-Tumour${2}")
        )
        .alias(col)
    )


def canonical_clone(value: float | int | str | None) -> str | None:
    """Coerce a mutation_cluster/clone value to canonical `clone{N}` or None.

    NaN and None → None (mutation was excluded from clustering).
    """
    if value is None:
        return None
    if isinstance(value, float) and math.isnan(value):
        return None
    if isinstance(value, str):
        if value.startswith("clone"):
            return value
        try:
            n = int(float(value))
        except ValueError as exc:
            raise ValueError(f"Cannot coerce clone value {value!r}") from exc
        return f"clone{n}"
    if isinstance(value, (int, float)):
        return f"clone{int(value)}"
    raise TypeError(f"Unsupported clone type: {type(value).__name__}")


def canonical_clone_expr(col: str) -> pl.Expr:
    """Polars expression version — vectorised coercion.

    Assumes ``col`` is a numeric mutation_cluster column (Float or Int) that
    may contain nulls / NaNs.
    """
    return (
        pl.when(pl.col(col).is_null() | pl.col(col).is_nan())
        .then(None)
        .otherwise(pl.format("clone{}", pl.col(col).cast(pl.Int64)))
        .alias(col)
    )


def canonical_chr(value) -> str | None:
    """Strip a leading ``chr`` prefix and coerce to string.

    Accepts int (``17``), string (``"chr17"``, ``"17"``, ``"X"``), or None.
    Returns None for null. `MT` and `chrM` both collapse to `MT`.
    """
    if value is None:
        return None
    if isinstance(value, float) and math.isnan(value):
        return None
    s = str(value)
    if s.startswith("chr"):
        s = s[3:]
    if s == "M":
        s = "MT"
    return s


def canonical_chr_expr(col: str) -> pl.Expr:
    """Polars expression version — vectorised prefix strip + ``M`` → ``MT``."""
    return (
        pl.col(col)
        .cast(pl.Utf8)
        .str.strip_prefix("chr")
        .replace({"M": "MT"})
        .alias(col)
    )


def patient_from_tumour_id(tumour_id: str) -> str:
    """`LTX0001-Tumour1` → `LTX0001`. Assumes canonical form."""
    m = _WGD_TUMOUR_RE.match(tumour_id)
    if not m:
        raise ValueError(f"Not canonical tumour_id: {tumour_id!r}")
    return m.group(1)


def tumour_ordinal_from_tumour_id(tumour_id: str) -> int:
    """`LTX0001-Tumour2` → 2. Assumes canonical form."""
    m = _WGD_TUMOUR_RE.match(tumour_id)
    if not m:
        raise ValueError(f"Not canonical tumour_id: {tumour_id!r}")
    return int(m.group(2))


def build_tumour_id(patient_id: str, ordinal: int) -> str:
    """`LTX0001`, 2 → `LTX0001-Tumour2`."""
    return f"{patient_id}-Tumour{ordinal}"


def assert_all_canonical(tumour_ids: Iterable[str]) -> None:
    """Guard: fail loudly if any tumour_id isn't in canonical form."""
    bad = [t for t in tumour_ids if not _WGD_TUMOUR_RE.match(t)]
    if bad:
        raise ValueError(
            f"Non-canonical tumour_ids: {bad[:5]!r} (showing first 5 of {len(bad)})"
        )
