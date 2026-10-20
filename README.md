# tx_data

Semantic data layer + curated store for TRACERx-scale cohort analysis.

## Idea

Scattered files on the HPC (`/Nemo`) are converted, via idempotent scripts, into a canonical Parquet-based store described by a YAML catalog. A DuckDB file attaches over the Parquet for fast SQL access. The result is designed for both human analysts and LLM assistants.

## Path abstraction

Code never hardcodes `/Nemo`. Paths resolve via `tx_data.paths.data_root()`, which reads:

1. Environment variable `TX_DATA_ROOT` (highest priority)
2. `config/paths.yml` `data_root:` (default: `nemo_mock`)

`nemo_mock/` mirrors the `/Nemo` directory structure with synthetic data, so the same code runs locally and on the cluster:

```bash
# Local dev — uses nemo_mock/
python scripts/build_clinical.py

# Explicit build modes for the full build
pixi build sandbox
pixi build prod

# Production on HPC
export TX_DATA_ROOT=/Nemo
python scripts/build_clinical.py
```

## Layout

- `nemo_mock/` — synthetic mock data, mirrors `/Nemo` structure
- `catalog/` — YAML per table: columns, semantics, NULL rules, provenance
- `scripts/` — idempotent source → canonical Parquet conversion scripts
- `data/` — derived Parquet (gitignored, regenerable)
- `db/` — DuckDB file (gitignored, regenerable)
- `src/tx_data/` — Python helpers (path resolver, catalog loader)
- `tests/` — smoke + validation tests
- `config/paths.yml` — data-root config

## Getting started

```bash
pip install -e .[dev]
pytest
```
