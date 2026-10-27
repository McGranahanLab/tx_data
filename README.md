# tx_data

Semantic data layer + curated store for TRACERx-scale cohort analysis.

## Idea

Scattered files on the HPC (`/Nemo`) are converted, via idempotent scripts, into a canonical Parquet-based store described by a YAML catalog. A DuckDB file attaches over the Parquet for fast SQL access. The result is designed for both human analysts and LLM assistants.

## Cohorts

The lab has two sequencing cohorts, and every table belongs to exactly one:

| `cohort` | Assay | Main Nemo root |
|---|---|---|
| `wes` | multi-region whole-exome (TRACERx 842) | `/nemo/project/proj-tracerx-lung/tctProjects/lungTx/Tx842/release/tracerx_842` |
| `wgs` | whole-genome, mostly single-sample | `/nemo/project/proj-tracerx-lung/tctProjects/CN-CCF/wgs_resources` |

Static reference tables (hg19) use `cohort: reference`. The cohort is
declared in `catalog/_sources.yml` (per source) and in each
`catalog/<table>.yml`; tests keep the two in sync. All tables built so far
are WES. Do not join WES and WGS tables without doing so explicitly.

## Path abstraction

Code never hardcodes `/nemo`. Every source root in `catalog/_sources.yml` is
an absolute `/nemo/...` path, re-rooted under `tx_data.paths.nemo_root()`:

1. Environment variable `TX_NEMO_ROOT` (highest priority)
2. `config/paths.yml` `nemo_root:` (default: `nemo_mock`)

`nemo_mock/` mirrors `/nemo` one-to-one with synthetic data
(`/nemo/<rest>` ↔ `nemo_mock/<rest>`), so the WES tree, the WGS tree and
alt roots elsewhere on Nemo all resolve with the same rule, locally and on
the cluster:

```bash
# Local dev — uses nemo_mock/
python scripts/build_clinical.py

# Explicit build modes for the full build
pixi run build sandbox
pixi run build prod

# Production on HPC
export TX_NEMO_ROOT=/nemo
python scripts/build_clinical.py
```

`TX_DATA_ROOT` (the old WES-only override) is retired: setting it raises an
error so a stale value can't mix mock and real files.

## Layout

- `nemo_mock/` — synthetic mock data, mirrors `/nemo` (`nemo_mock/_cohort_manifest.json` is mock-only metadata for the WES mock cohort)
- `catalog/` — YAML per table: columns, semantics, NULL rules, provenance
- `scripts/` — idempotent source → canonical Parquet conversion scripts
- `data/` — derived Parquet (gitignored, regenerable)
- `db/` — DuckDB file (gitignored, regenerable)
- `src/tx_data/` — Python helpers (path resolver, catalog loader)
- `tests/` — smoke + validation tests
- `config/paths.yml` — `nemo_root` + id-convention config

## Getting started

```bash
pip install -e .[dev]
pytest
```
