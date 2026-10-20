"""Path resolution for tx_data.

`data_root()` picks (in order):
  1. Environment variable ``TX_DATA_ROOT``.
  2. ``data_root`` field in ``config/paths.yml``.
  3. Repo-local ``nemo_mock/``.

Per-table source paths live in ``catalog/_sources.yml`` and are resolved
against either the main ``data_root()`` or one of the ``alt_source_roots``
declared there.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Mapping

import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_PATHS_YML = REPO_ROOT / "config" / "paths.yml"
DEFAULT_MOCK_ROOT = REPO_ROOT / "nemo_mock"
DEFAULT_PROD_ROOT = Path(
    "/nemo/project/proj-tracerx-lung/tctProjects/lungTx/Tx842/release/tracerx_842"
)



def data_root() -> Path:
    """Resolve the cohort source root."""
    env = os.environ.get("TX_DATA_ROOT")
    if env:
        return Path(env).expanduser()

    if DEFAULT_PATHS_YML.is_file():
        cfg = yaml.safe_load(DEFAULT_PATHS_YML.read_text()) or {}
        val = cfg.get("data_root")
        if val:
            p = Path(val).expanduser()
            return p if p.is_absolute() else (REPO_ROOT / p)

    return DEFAULT_MOCK_ROOT


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
    """True when `data_root()` is the repo-local mock tree (nemo_mock/)."""
    try:
        return data_root().resolve() == DEFAULT_MOCK_ROOT.resolve()
    except OSError:
        return False


def build_root_for_mode(mode: str | None) -> Path | None:
    """Return the source root requested by the build command.

    `None` preserves the existing `data_root()` resolution order.
    "sandbox" forces the repo-local mock tree.
    "prod" forces the HPC release path.
    """
    if mode is None:
        return None
    if mode == "sandbox":
        return DEFAULT_MOCK_ROOT
    if mode == "prod":
        return DEFAULT_PROD_ROOT
    raise ValueError(f"Unknown build mode: {mode!r}")


def _alt_root(alias: str, alt_source_roots: Mapping[str, str]) -> Path:
    """Resolve an alt-source-root alias to an absolute path.

    On the HPC the alias maps to an absolute path from ``_sources.yml``.
    Locally we substitute the mock tree (``nemo_mock/alt/<alias>/``) so
    the same table can be built either way.
    """
    root = data_root()
    if root == DEFAULT_MOCK_ROOT or (root.is_dir() and (root / "alt" / alias).exists()):
        return root / "alt" / alias

    if alias not in alt_source_roots:
        raise KeyError(f"alt_source_root alias {alias!r} not in _sources.yml")
    return Path(alt_source_roots[alias]).expanduser()


def resolve_source_path(
    relative_path: str,
    *,
    source_root: str | None = None,
    alt_source_roots: Mapping[str, str] | None = None,
) -> Path:
    """Resolve a per-table `relative_path` against the right root."""
    if source_root is None:
        return data_root() / relative_path
    return _alt_root(source_root, alt_source_roots or {}) / relative_path


def resolve_source_pattern(
    relative_path_pattern: str,
    *,
    source_root: str | None = None,
    alt_source_roots: Mapping[str, str] | None = None,
    **fmt: str,
) -> Path:
    """Same as `resolve_source_path` for patterns with `{placeholders}`."""
    return resolve_source_path(
        relative_path_pattern.format(**fmt),
        source_root=source_root,
        alt_source_roots=alt_source_roots,
    )
