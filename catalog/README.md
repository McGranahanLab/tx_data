# Catalog

One YAML file per canonical table. Files starting with `_` are templates/examples and ignored by `catalog.load_all()`.

## Fields per table

- `name` — canonical table name (snake_case)
- `description` — plain-English summary
- `provenance` — source path(s) + version + notes on how it was produced
- `primary_key` — column(s) uniquely identifying a row
- `columns` — list, each with:
  - `name`
  - `type` (SQL/DuckDB type)
  - `description`
  - `allowed_values` (if enum)
  - `null_semantics` (what NULL means: not_measured / not_applicable / unknown / measured_as_absent)
  - `notes` (edge cases, gotchas)
- `relationships` — foreign keys to other tables

A `_template.yml` will be added with the first real table.
