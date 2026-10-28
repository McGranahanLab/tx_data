# tx_data — implementation progress

Session-safe log so any future session can pick up cleanly. Update as you go.

## Goal

Semantic data layer over `/Nemo` TRACERx cohort files: YAML catalog + Parquet
canonical store + synthetic `nemo_mock/` for local dev. Downstream projects
(currently `mut_essential_wgd`) consume via `import tx_data`.

## Design (locked)

- **Cohorts** (2026-10-27): two sequencing cohorts, `wes` (multi-region exome,
  Tx842 release `lungTx/Tx842/release/tracerx_842`) and `wgs` (whole genome,
  mostly `CN-CCF/wgs_resources`). `_sources.yml` `cohorts:` declares each root;
  every source entry and every `catalog/<table>.yml` carries `cohort:`
  (`wes` / `wgs` / `reference`). All tables so far are `wes`; hg19 tables are
  `reference`. The build manifest records `table_cohorts`. Genome builds:
  WES hg19 (moving to hg38 soon), WGS hg38. WGS tables are named `wgs_<table>`;
  WGS mocks come from `scripts/generate_mocks_wgs.py` (own seed, own
  `nemo_mock/_cohort_manifest_wgs.json`). First WGS table:
  `wgs_battenberg_inventory` (2026-10-27).
- **WGS tumour labels** (2026-10-27): WGS sources (inventory `PATIENT`, Battenberg
  folder names) label tumours as bare `LTX0001` (single tumour) or
  `LTX0002_tumour<n>`. Canonical WGS tables standardise to WES form via
  `tx_data.ids.wgs_tumour_id` (`LTX0001-Tumour1`) and keep the raw string as
  `wgs_tumour_label` for locating files.
- **WGS CN tables** (2026-10-27): `wgs_consensus_segments` (fractional sample-level
  CN from Battenberg `tree_input/`) and `wgs_clone_proportions` (ALPACA input
  cp_table). Both discover tumours from the built inventory; sample names in
  other sources are matched via `tx_data.wgs_samples.match_sample`.
  2026-10-28: segments now come from the MSAI-adjusted consensus file, and
  `wgs_tree_clones` (ALPACA `tree_paths.json`) defines single-clone tumours
  (`n_clones_tree = 1`). Pending: ALPACA output (integer clone CN), remaining
  tree_input / alpaca input files.
- **Path resolution**: every root in `_sources.yml` is an absolute `/nemo/...`
  path, re-rooted under `tx_data.paths.nemo_root()` (`TX_NEMO_ROOT` env var,
  then `config/paths.yml` `nemo_root:`, default `nemo_mock/`). `nemo_mock/`
  mirrors `/nemo` 1:1, so WES, WGS and alt roots all mock the same way.
  `TX_DATA_ROOT` is retired and raises if set. `data_root()` still returns the
  WES release root.
- **Alt source roots**: `catalog/_sources.yml` `alt_source_roots:` map lets a
  table declare `source_root: <alias>` to live outside its cohort's main tree
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
├── paths.py          # nemo_root(), nemo_path(), data_root() (WES)
├── catalog.py        # load_catalog(name) → dict; load_all() → dict[str,dict]
├── sources.py        # load_sources(); resolve_source(name); source_cohort(); cohort_root()
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
