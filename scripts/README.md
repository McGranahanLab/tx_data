# Conversion scripts

One script per canonical table. Each script:

1. Resolves source paths via `tx_data.paths.data_root()` (never hardcodes `/Nemo`)
2. Reads source file(s)
3. Cleans / renames / casts
4. Writes canonical Parquet to `data/<table>.parquet`
5. Writes a metadata sidecar to `data/_metadata/<table>.json` with source hashes + timestamp + script version

Scripts must be **idempotent** — re-running should produce byte-identical Parquet output (modulo the timestamp in the sidecar).
