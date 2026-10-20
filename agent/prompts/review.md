# Commit / review checklist

Before every commit:

- [ ] `ruff check` and `ruff format --check` pass.
- [ ] `pytest` passes.
- [ ] If a new table was added: mock file + conversion script + catalog entry + test all landed in the same commit.
- [ ] No hardcoded `/Nemo` paths anywhere in code.
- [ ] No real patient data or real IDs committed.
- [ ] Catalog entries have `provenance`, `primary_key`, and per-column `null_semantics`.
- [ ] ETL scripts are idempotent (re-running produces the same Parquet).
- [ ] Commit message is imperative-mood, first line ≤ 72 chars.
- [ ] `labbook/YYYY-MM-DD.md` updated: records what changed and why.
- [ ] `labbook/STATE.md` updated if the current-state summary changed.
- [ ] Any decision worth remembering proposed as a new `DECISIONS.md` entry inside the daily log — the human moves it into `DECISIONS.md` later.
