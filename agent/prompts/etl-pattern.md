# Standard ETL script shape

One script per canonical table at `scripts/build_<table>.py`.

## Required behaviour

1. Resolves source paths via `tx_data.paths.data_root()`. Never hardcodes `/Nemo`.
2. Reads source file(s).
3. Cleans + renames + casts to the catalog schema.
4. Writes canonical Parquet to `data/<table>.parquet`.
5. Writes a metadata sidecar to `data/_metadata/<table>.json` with:
   - `source_files` — list of `{path, sha256}` per source file
   - `script_path` — repo-relative path to this script
   - `built_at` — UTC ISO 8601 timestamp
   - `row_count`, `column_count`
6. Prints a one-line summary on success:
   `built <table>: N rows, M cols -> data/<table>.parquet`

## Idempotency

Re-running on unchanged inputs must produce byte-identical Parquet (modulo the timestamp in the sidecar). Concretely:

- Sort output rows by primary key before write.
- Fix explicit column ordering.
- Seed any sampling / RNG deterministically.
- Fixed Parquet compression (`zstd`, level 3) and row-group size.

## Skeleton

```python
"""Build canonical <table>.parquet from source under data_root."""
from __future__ import annotations

import polars as pl

from tx_data.paths import data_root, derived_dir

SOURCE_REL = "path/under/data_root/<file>.csv"
TABLE = "<table>"
PRIMARY_KEY = ["patient_id"]


def build() -> pl.DataFrame:
    src = data_root() / SOURCE_REL
    df = pl.read_csv(src)
    # rename, cast, clean
    return df.sort(PRIMARY_KEY)


def main() -> None:
    df = build()
    out = derived_dir() / f"{TABLE}.parquet"
    out.parent.mkdir(parents=True, exist_ok=True)
    df.write_parquet(out, compression="zstd", compression_level=3)
    # write sidecar (see helper — to be added)
    print(f"built {TABLE}: {df.height} rows, {df.width} cols -> {out}")


if __name__ == "__main__":
    main()
```

A `write_sidecar()` helper will land in `src/tx_data/etl.py` alongside the first real ETL script, so every script writes sidecars consistently.

## Testing

- Unit test the pure `build()` function against mock inputs.
- Integration test the full `main()` end-to-end on `nemo_mock/`.
