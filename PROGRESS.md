# tx_data — implementation progress

Session-safe log so any future session can pick up cleanly. Update as you go.

## Goal

Semantic data layer over `/Nemo` TRACERx cohort files: YAML catalog + Parquet
canonical store + synthetic `nemo_mock/` for local dev. Downstream projects
(currently `mut_essential_wgd`) consume via `import tx_data`.

## Design (locked)

- **Path resolution**: `tx_data.paths.data_root()` reads `TX_DATA_ROOT` env var
  first, then falls back to `config/paths.yml`. Default `nemo_mock/`.
- **Alt source roots**: `catalog/_sources.yml` `alt_source_roots:` map lets a
  table declare `source_root: <alias>` to live outside the main release tree
  (used for AlphaMissense at `CN-CCF/alphamissense/…`).
- **Canonical IDs**:
  - `tumour_id` canonical form: `LTX0001-Tumour1` (matches wgd_calls).
  - `clone` canonical form: string `clone1`, `clone10` (matches wgd_calls).
  - Loaders coerce muttable's `LTX0001_tumour1` / `1.0` on import.
- **Grain rules** (see catalog per-table):
  - muttable: per (mutation × sample). Filter by `is_present`.
  - wgd_calls: per (tumour × clone). Filter by `status == 'resolved'`.
  - alphamissense: raw per (variant × transcript × sample); dedup to per variant.
  - clinical: per patient; loader unpivots two-lesion columns to per-tumour.

## Status

| Step | Status | Files |
|---|---|---|
| Catalog: muttable | ✅ | `catalog/muttable.yml` |
| Catalog: wgd_calls | ✅ | `catalog/wgd_calls.yml` |
| Catalog: alphamissense | ✅ | `catalog/alphamissense.yml` |
| Catalog: clinical | ✅ | `catalog/clinical.yml` |
| Sources registry | ✅ | `catalog/_sources.yml` |
| Path resolver | ✅ | `src/tx_data/paths.py` |
| Catalog loader | ✅ | `src/tx_data/catalog.py` |
| Sources loader | ✅ | `src/tx_data/sources.py` |
| Normalise helpers | ✅ | `src/tx_data/normalize.py` |
| Build: muttable | ✅ | `src/tx_data/builds/muttable.py` |
| Build: wgd_calls | ✅ | `src/tx_data/builds/wgd_calls.py` |
| Build: alphamissense | ✅ | `src/tx_data/builds/alphamissense.py` |
| Build: clinical | ✅ | `src/tx_data/builds/clinical.py` |
| Mock generator | ✅ | `scripts/generate_mocks.py` |
| Tests | ✅ | `tests/test_pipeline.py` (6 passing) |

## How to run (pixi)

```bash
cd projects/tx_data
pixi install       # once
pixi run mocks     # regenerate nemo_mock/
pixi run build     # → data/*.parquet
pixi run test      # 12 tests, ~0.6s
pixi run all       # all three, chained
```

Cluster / HPC: pixi works there too; for the final production run we'll wrap
this same env in an Apptainer container.

Mock cohort: 10 patients, 11 tumours. Covers every wgd_calls `class`
(no_wgd, clonal_wgd, mut_supported, ploidy_only), sub_class (single, parallel,
sequential), status (resolved, unresolved, needs_follow_up), a multi-tumour
patient (LTX0008), and a tumour with `number_of_gds_at_clone == 2` at MRCA.

Legend: ✅ done · ⏳ in progress or todo.

## Layout to build

```
src/tx_data/
├── __init__.py
├── paths.py          # data_root(), alt_root(alias)
├── catalog.py        # load_catalog(name) → dict; load_all() → dict[str,dict]
├── sources.py        # load_sources() → dict; resolve_source(name)
├── normalize.py      # canonical tumour_id / clone helpers
└── builds/
    ├── __init__.py
    ├── _base.py      # build_table() common runner
    ├── muttable.py
    ├── wgd_calls.py
    ├── alphamissense.py
    └── clinical.py

scripts/
├── build_muttable.py
├── build_wgd_calls.py
├── build_alphamissense.py
├── build_clinical.py
├── build_all.py
└── generate_mocks.py
```

## Resume playbook

1. Check this file's status table → find first `⏳` row.
2. Read the corresponding catalog YAML (`catalog/<name>.yml`) for column
   spec + normalisation rules + join keys.
3. Look at the previous build script for the pattern (same shape across
   tables: read source per _sources.yml → apply catalog columns + coercions
   → write Parquet to `data/<name>.parquet`).
4. Update this file's status row when done. Commit.

## Cross-references

- **Downstream project**: `../mut_essential_wgd/` (see its own PROGRESS.md).
- **Analysis plan**: `../mut_essential_wgd/notes/analysis_plan.md`.
- **Deferred items**: `../mut_essential_wgd/notes/TODO.md`.
