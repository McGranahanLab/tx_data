# tx_data — agent operating manual

You are working on `tx_data`: a semantic data layer + curated store for TRACERx-scale cohort analysis. Your job is to ingest scattered files under `/Nemo` (mirrored locally as `nemo_mock/`), produce canonical Parquet, describe each table in a YAML catalog, and iterate on the assistant that queries this store.

## Session-start ritual

At the start of every session, read in this order:

1. This file (`agent/SYSTEM.md`) — always
2. `labbook/STATE.md` — current situation
3. `labbook/YYYY-MM-DD.md` for today and yesterday — recent raw context
4. `agent/prompts/<relevant>.md` on demand for the specific task
5. `catalog/*.yml` on demand for the tables involved

## Golden rules

- **Never hardcode `/Nemo` paths.** Resolve via `tx_data.paths.data_root()`.
- **`nemo_mock/` mirrors `/Nemo` exactly.** Same directory structure, synthetic values only.
- **ETL scripts are idempotent.** Re-running produces byte-identical Parquet (modulo the sidecar timestamp).
- **The catalog is the truth.** Every canonical table has a YAML entry with per-column semantics, NULL rules, and provenance.
- **Never propagate source column names.** Rename to snake_case during ETL.
- **All dates in UTC, ISO 8601.**
- **Never invent field names or values.** If a semantic is unclear, ask.
- **All non-trivial work goes in the lab book.** See `labbook/`.

## Where things live

- Operating rules → `agent/SYSTEM.md` (this file)
- Detail prompts → `agent/prompts/*.md`
  - Coding conventions → `agent/prompts/conventions.md`
  - Semantic-layer format → `agent/prompts/data-model.md`
  - ETL script shape → `agent/prompts/etl-pattern.md`
  - Data governance → `agent/prompts/privacy.md`
  - Commit / review checklist → `agent/prompts/review.md`
- Current project state → `labbook/STATE.md`
- Architecture decisions → `labbook/DECISIONS.md`
- Daily raw log → `labbook/YYYY-MM-DD.md`
- Table catalogs → `catalog/<table>.yml`
- ETL scripts → `scripts/build_<table>.py`
- Mock data → `nemo_mock/<mirrored-path>/...`
- Python helpers → `src/tx_data/`

## Format rules

- ISO 8601 dates in filenames and body (`2026-07-30`, `2026-07-30 13:30 UTC`).
- Never write "today", "yesterday", "next week" in lab-book entries — use absolute dates. Relative time rots.
- Every non-trivial claim links to a commit hash or file path.
- No markdown tables in the lab book — use bullet lists. Tables rot.

## Write access

- **`labbook/YYYY-MM-DD.md`** — write freely.
- **`labbook/STATE.md`** — update at end of each work session; changes visible in git diff.
- **`labbook/DECISIONS.md`** — propose new entries in the daily log; the human moves them here after confirming.
- **`agent/*.md`** — propose edits; the human reviews. Policy files, gated.
- **`catalog/*.yml`** — write freely for new tables; propose diffs for existing tables.
- **`scripts/*.py`** — write freely.
- **`nemo_mock/**`** — write freely.

## Data governance

Real data lives on `/Nemo` and never leaves the HPC. Local dev uses `nemo_mock/` only. Full rules in `agent/prompts/privacy.md`.
