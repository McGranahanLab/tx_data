# Coding conventions

## Python

- Python ≥ 3.10, type hints where sensible.
- Line length 100.
- `from __future__ import annotations` at the top of new modules.
- Prefer `pathlib.Path` over string paths.
- Prefer `polars` over `pandas` for new ETL scripts; `pandas` is fine if the script is genuinely trivial.
- Use `duckdb` for anything SQL-shaped.
- Format + lint with `ruff`.
- Tests via `pytest`.

## Naming

- snake_case for identifiers and column names.
- Table names in the catalog: singular noun (`clinical`, `variant`, not `clinicals` / `variants`).
- Script names: `build_<table>.py`.
- Catalog files: `<table>.yml`.
- No abbreviations in column names unless universally known (`chr`, `id`, `gc`).

## Commit messages

- Imperative mood, present tense ("Add clinical table", not "Added" or "Adds").
- First line ≤ 72 chars.
- Body wraps at 100 chars.
- Reference the table or module changed if it's not obvious from the first line.

## Testing

- Every catalog entry gets a test that loads it and checks basic invariants (row count > 0, primary key unique, no unexpected NULLs).
- Every ETL script is tested end-to-end against `nemo_mock/`.
- Fast: the full suite should complete in seconds on mock data.
