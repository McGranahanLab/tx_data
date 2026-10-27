"""Path resolution for tx_data.

Every source location is declared in ``catalog/_sources.yml`` as an
absolute Nemo path (``/nemo/...``). At build time it is re-rooted under
``nemo_root()``, which picks (in order):
  1. Environment variable ``TX_NEMO_ROOT``.
  2. ``nemo_root`` field in ``config/paths.yml``.
  3. Repo-local ``nemo_mock/``.

``nemo_mock/`` mirrors ``/nemo`` one-to-one, so the WES release tree, the
WGS resources tree and every alt root resolve with the same rule:
``/nemo/<rest>`` → ``<nemo_root>/<rest>``.

``data_root()`` is kept for callers that want the WES release root
(``tracerx_842``); it is derived from ``nemo_root()``.
"""

from __future__ import annotations

import os
from pathlib import Path, PurePosixPath
from typing import Mapping

import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_PATHS_YML = REPO_ROOT / "config" / "paths.yml"
DEFAULT_MOCK_ROOT = REPO_ROOT / "nemo_mock"
NEMO_PREFIX = PurePosixPath("/nemo")
DEFAULT_PROD_ROOT = Path(NEMO_PREFIX)
WES_RELEASE_ROOT = PurePosixPath(
    "/nemo/project/proj-tracerx-lung/tctProjects/lungTx/Tx842/release/tracerx_842"
)


def _check_legacy_env() -> None:
    # TX_DATA_ROOT used to point at the WES release root only. Alt roots and
    # WGS roots can't be derived from it, so a stale setting would silently
    # mix prod and mock files. Fail loudly instead.
    if os.environ.get("TX_DATA_ROOT"):
        raise RuntimeError(
            "TX_DATA_ROOT is no longer supported. Set TX_NEMO_ROOT to the Nemo "
            "mount (`/nemo` on the HPC) or use `pixi run build prod|sandbox`."
        )


def nemo_root() -> Path:
    """Resolve the directory that stands in for ``/nemo``."""
    _check_legacy_env()
    env = os.environ.get("TX_NEMO_ROOT")
    if env:
        return Path(env).expanduser()

    if DEFAULT_PATHS_YML.is_file():
        cfg = yaml.safe_load(DEFAULT_PATHS_YML.read_text()) or {}
        val = cfg.get("nemo_root")
        if val:
            p = Path(val).expanduser()
            return p if p.is_absolute() else (REPO_ROOT / p)

    return DEFAULT_MOCK_ROOT


def nemo_path(abs_path: str | PurePosixPath) -> Path:
    """Re-root an absolute ``/nemo/...`` path under ``nemo_root()``."""
    p = PurePosixPath(abs_path)
    try:
        rest = p.relative_to(NEMO_PREFIX)
    except ValueError:
        raise ValueError(f"Source path {str(p)!r} is not under {NEMO_PREFIX}") from None
    return nemo_root() / rest


def data_root() -> Path:
    """WES release root (``.../tracerx_842``) under the current ``nemo_root()``."""
    return nemo_path(WES_RELEASE_ROOT)


def id_convention() -> str:
    """Global tumour-id convention: 'ltx' (default) or 'cruk'.

    Precedence:
      1. ``TX_ID_CONVENTION`` env var.
      2. ``id_convention`` field in ``config/paths.yml``.
      3. ``'ltx'`` (tx_data legacy).
    """
    env = os.environ.get("TX_ID_CONVENTION")
    if env:
        v = env.strip().lower()
    elif DEFAULT_PATHS_YML.is_file():
        cfg = yaml.safe_load(DEFAULT_PATHS_YML.read_text()) or {}
        v = str(cfg.get("id_convention", "ltx")).strip().lower()
    else:
        v = "ltx"

    if v not in ("ltx", "cruk"):
        raise ValueError(
            f"id_convention must be 'ltx' or 'cruk' (got {v!r}). "
            "Set TX_ID_CONVENTION or edit config/paths.yml."
        )
    return v


def data_root_is_mock() -> bool:
    """True when ``nemo_root()`` is the repo-local mock tree (nemo_mock/)."""
    try:
        return nemo_root().resolve() == DEFAULT_MOCK_ROOT.resolve()
    except OSError:
        return False


def build_root_for_mode(mode: str | None) -> Path | None:
    """Return the ``nemo_root`` requested by the build command.

    `None` preserves the existing `nemo_root()` resolution order.
    "sandbox" forces the repo-local mock tree.
    "prod" forces the real ``/nemo`` mount.
    """
    if mode is None:
        return None
    if mode == "sandbox":
        return DEFAULT_MOCK_ROOT
    if mode == "prod":
        return DEFAULT_PROD_ROOT
    raise ValueError(f"Unknown build mode: {mode!r}")


def _root_for(
    source_root: str | None,
    cohort_root: str | None,
    alt_source_roots: Mapping[str, str],
) -> Path:
    """Resolve a table's root directory.

    ``source_root`` None → the table's cohort root. ``in_repo`` → REPO_ROOT,
    for cohort-independent static reference data committed alongside the
    code. Any other value is an ``alt_source_roots`` alias.
    """
    if source_root == "in_repo":
        return REPO_ROOT
    if source_root is None:
        if cohort_root is None:
            raise ValueError("Table has no `source_root` and no cohort root")
        return nemo_path(cohort_root)
    if source_root not in alt_source_roots:
        raise KeyError(f"alt_source_root alias {source_root!r} not in _sources.yml")
    return nemo_path(alt_source_roots[source_root])


def resolve_source_path(
    relative_path: str,
    *,
    source_root: str | None = None,
    cohort_root: str | None = None,
    alt_source_roots: Mapping[str, str] | None = None,
) -> Path:
    """Resolve a per-table `relative_path` against the right root."""
    return _root_for(source_root, cohort_root, alt_source_roots or {}) / relative_path


def resolve_source_pattern(
    relative_path_pattern: str,
    *,
    source_root: str | None = None,
    cohort_root: str | None = None,
    alt_source_roots: Mapping[str, str] | None = None,
    **fmt: str,
) -> Path:
    """Same as `resolve_source_path` for patterns with `{placeholders}`."""
    return resolve_source_path(
        relative_path_pattern.format(**fmt),
        source_root=source_root,
        cohort_root=cohort_root,
        alt_source_roots=alt_source_roots,
    )
