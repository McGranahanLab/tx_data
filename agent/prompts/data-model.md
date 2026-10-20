# Semantic-layer format

Each canonical table has one YAML file under `catalog/`. Files starting with `_` are templates and ignored by `catalog.load_all()`.

## Required fields

- `name` — canonical table name (snake_case, singular)
- `description` — one-paragraph plain-English summary
- `provenance` — object:
  - `source_paths` — list of relative paths under `data_root`, or a glob
  - `source_format` — `csv`, `tsv`, `xlsx`, `rds`, `vcf`, `parquet`, `custom`
  - `pipeline_version` — upstream pipeline version, if known
  - `notes` — anything else worth capturing
- `primary_key` — column name, or list of columns for composite keys
- `columns` — list, each with:
  - `name`
  - `type` — DuckDB type (`INTEGER`, `VARCHAR`, `DOUBLE`, `BOOLEAN`, `DATE`, `TIMESTAMP`, `LIST<...>`, etc.)
  - `description` — plain-English
  - `allowed_values` — for enums; omit otherwise
  - `null_semantics` — one of: `not_measured`, `not_applicable`, `unknown`, `measured_as_absent`, `no_nulls_allowed`
  - `notes` — edge cases, historical baggage, gotchas
- `relationships` — optional list of foreign-key links to other tables

## Optional fields

- `row_count_expected` — rough expectation for sanity checks
- `ontology_refs` — links to external ontologies (OncoTree, HGNC, SO, HPO)
- `version` — version of this catalog entry itself

## Style rules

- Descriptions target a smart reader who has never seen this dataset. No lab jargon without expansion.
- `notes` capture the "if you don't know this, you will get it wrong" points. This is where the ancestry-style gotchas live (`"Asian"` not a value; use `East_Asian` / `South_Asian`).
- Provenance is not optional. If we can't say where a column came from, we can't trust it.
- Enum values in `allowed_values` are snake_case in the canonical store, even if the source uses other casing. Document the mapping in `notes`.

A working example lives at `catalog/_template.yml` (added with the first real table).
