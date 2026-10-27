"""Tumour-id flavour handling: parse, translate, emit.

Multiple ID conventions appear across the Tx842 ecosystem. We auto-detect
the input flavour rather than requiring each builder to declare it up
front — a new source file often lands in a subtly different form, and
brittle assumptions have already burned us once (User 2026-09-23).

Recognised forms:

  Patient prefix  ×  Suffix pattern

    LTX<digits>   +  '-Tumour<n>'    -> canonical LTX (wgd_calls, alpaca, ...)
    LTX<digits>   +  '_tumour<n>'    -> muttable spelling (case-insensitive)
    LTX<digits>   +  '_<n>'          -> cloneInfo spelling (new 2026-09-23)
    LTX<digits>   +  '-Cluster<n>'   -> some clinical assets
    CRUK<digits>  +  '_<n>'          -> CRUK canonical
    CRUK<digits>  +  '-Cluster<n>'   -> CRUK cluster split (Zenodo tree_data)
    CRUK<digits>  +  '-Tumour<n>'    -> hypothetical, kept for symmetry

All variants encode the same (patient, tumour_ordinal) grain — the choice
of separator and word is cosmetic. The output form is dictated by
``tx_data.paths.id_convention()``:

  * ``'ltx'``  → ``LTX0001-Tumour1``
  * ``'cruk'`` → ``CRUK0001_1``

If input and target patient prefixes differ, the ``id_crosswalk`` table
supplies the patient-level translation. Suffix ordinal is preserved
verbatim.
"""

from __future__ import annotations

import re
from functools import lru_cache
from pathlib import Path

import polars as pl

from tx_data.builds._base import canonical_output_path
from tx_data.paths import id_convention

# Ordered so we try the specific forms first ('-Tumour', '-Cluster',
# '_tumour') before the bare-underscore form (which is a superset of the
# suffix character).
_PATTERNS: tuple[re.Pattern[str], ...] = (
    re.compile(r"^(LTX\d+)-Tumour(\d+)$"),
    re.compile(r"^(LTX\d+)-Cluster(\d+)$"),
    re.compile(r"^(LTX\d+)_tumour(\d+)$", re.IGNORECASE),
    re.compile(r"^(LTX\d+)_(\d+)$"),
    re.compile(r"^(CRUK\d+)-Tumour(\d+)$"),
    re.compile(r"^(CRUK\d+)-Cluster(\d+)$"),
    re.compile(r"^(CRUK\d+)_(\d+)$"),
)


def parse_tumour_id(value: str) -> tuple[str, int]:
    """Return ``(patient_id, tumour_ordinal)`` from any recognised form.

    Raises ``ValueError`` if the input doesn't match any known pattern.
    """
    for regex in _PATTERNS:
        m = regex.match(value)
        if m:
            return m.group(1), int(m.group(2))
    raise ValueError(f"Unrecognised tumour_id form: {value!r}")


@lru_cache(maxsize=1)
def _crosswalk_maps() -> tuple[dict[str, str], dict[str, str]]:
    """Load ``id_crosswalk`` and return ``(cruk_to_ltx, ltx_to_cruk)`` dicts."""
    parquet = Path(canonical_output_path("id_crosswalk"))
    if not parquet.is_file():
        raise FileNotFoundError(
            f"id_crosswalk parquet not found at {parquet}. "
            "Run `python scripts/build_id_crosswalk.py` first."
        )
    df = pl.read_parquet(parquet)
    if not {"ltx_patient_id", "cruk_patient_id"}.issubset(df.columns):
        raise ValueError(
            f"id_crosswalk parquet missing expected columns; got {df.columns!r}"
        )
    cruk_to_ltx_map = dict(zip(df["cruk_patient_id"], df["ltx_patient_id"]))
    ltx_to_cruk_map = dict(zip(df["ltx_patient_id"], df["cruk_patient_id"]))
    if len(cruk_to_ltx_map) != df.height or len(ltx_to_cruk_map) != df.height:
        raise ValueError(
            "id_crosswalk contains duplicate patient IDs; rebuild is needed."
        )
    return cruk_to_ltx_map, ltx_to_cruk_map


def _emit(patient: str, ordinal: int, target: str) -> str:
    """Emit canonical form for the given target convention."""
    if target == "ltx":
        return f"{patient}-Tumour{ordinal}"
    if target == "cruk":
        return f"{patient}_{ordinal}"
    raise ValueError(f"Unsupported target id_convention: {target!r}")


def apply_convention_scalar(value: str) -> str:
    """Translate any recognised tumour_id form into the configured convention.

    Auto-detects the input flavour, walks the crosswalk only if the patient
    prefix family (LTX vs CRUK) differs from the target.
    """
    patient, ordinal = parse_tumour_id(value)
    target = id_convention()

    if patient.startswith("LTX") and target == "cruk":
        _, ltx_to_cruk_map = _crosswalk_maps()
        if patient not in ltx_to_cruk_map:
            raise KeyError(f"LTX patient {patient!r} not in id_crosswalk")
        patient = ltx_to_cruk_map[patient]
    elif patient.startswith("CRUK") and target == "ltx":
        cruk_to_ltx_map, _ = _crosswalk_maps()
        if patient not in cruk_to_ltx_map:
            raise KeyError(f"CRUK patient {patient!r} not in id_crosswalk")
        patient = cruk_to_ltx_map[patient]
    # else: patient already in target family — just re-emit canonical.

    return _emit(patient, ordinal, target)


# WGS tumour label: the Battenberg inventory's `PATIENT` column and the
# Battenberg output folder names. It identifies a tumour, not a patient:
# single-tumour patients carry the bare patient id (`LTX0001`), multi-tumour
# patients get a `_tumour<n>` suffix (`LTX0002_tumour1`; `_tumour_2` also
# accepted). Unlike WES, single tumours have no suffix at all.
_WGS_TUMOUR_LABEL = re.compile(r"^((?:LTX|CRUK)\d+)(?:_tumour_?(\d+))?$", re.IGNORECASE)


def parse_wgs_tumour_label(label: str) -> tuple[str, int]:
    """Return ``(patient_id, tumour_ordinal)`` from a WGS tumour label.

    A bare patient id means the patient's only tumour, ordinal 1.
    Raises ``ValueError`` on anything else.
    """
    m = _WGS_TUMOUR_LABEL.match(label)
    if m is None:
        raise ValueError(f"Unrecognised WGS tumour label: {label!r}")
    return m.group(1).upper(), int(m.group(2) or 1)


def wgs_tumour_id(label: str) -> str:
    """WGS tumour label → canonical tumour_id (``LTX0001`` → ``LTX0001-Tumour1``)."""
    patient, ordinal = parse_wgs_tumour_label(label)
    return apply_convention_scalar(f"{patient}_{ordinal}")


def apply_convention_expr(col: str = "tumour_id") -> pl.Expr:
    """Polars expression: coerce a ``tumour_id`` column into the configured form."""
    return pl.col(col).map_elements(apply_convention_scalar, return_dtype=pl.String)


# ---- Backwards-compatible scalar helpers ------------------------------------


def cruk_to_ltx(value: str) -> str:
    """Force CRUK-form input to canonical LTX. Ignores ``id_convention``."""
    patient, ordinal = parse_tumour_id(value)
    if not patient.startswith("CRUK"):
        raise ValueError(f"cruk_to_ltx got non-CRUK input: {value!r}")
    cruk_to_ltx_map, _ = _crosswalk_maps()
    if patient not in cruk_to_ltx_map:
        raise KeyError(f"CRUK patient {patient!r} not in id_crosswalk")
    return _emit(cruk_to_ltx_map[patient], ordinal, "ltx")


def ltx_to_cruk(value: str) -> str:
    """Force LTX-form input to canonical CRUK. Ignores ``id_convention``."""
    patient, ordinal = parse_tumour_id(value)
    if not patient.startswith("LTX"):
        raise ValueError(f"ltx_to_cruk got non-LTX input: {value!r}")
    _, ltx_to_cruk_map = _crosswalk_maps()
    if patient not in ltx_to_cruk_map:
        raise KeyError(f"LTX patient {patient!r} not in id_crosswalk")
    return _emit(ltx_to_cruk_map[patient], ordinal, "cruk")
