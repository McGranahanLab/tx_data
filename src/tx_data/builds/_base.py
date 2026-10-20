"""Common helpers for build scripts."""

from __future__ import annotations

from pathlib import Path

from tx_data.paths import REPO_ROOT

DATA_DIR = REPO_ROOT / "data"


def canonical_output_path(name: str) -> Path:
    """Where a canonical Parquet for a table lives."""
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    return DATA_DIR / f"{name}.parquet"


def log(msg: str) -> None:
    """Tiny logger — one function so build output is easy to grep."""
    print(f"[tx_data.build] {msg}", flush=True)
