"""Build wgs_clone_proportions → canonical Parquet.

WGS counterpart of `clone_proportions`: ALPACA input clone composition
tables, one per tumour at `wgs_resources/alpaca/input/<label>/cp_table.csv`
(same wide layout as WES: N sample columns + trailing `clone` column).

Tumours are discovered from the built `wgs_battenberg_inventory`; the
folder name is the WGS tumour label (`LTX0001` / `LTX0002_tumour1`) and is
standardised to WES-style `tumour_id`. Sample columns are matched to the
inventory `sample_name_hash` within the tumour; the source header is kept
as `sample_id`.

A per-tumour clone count (and so the single-clone flag) is a downstream
decision — e.g. whether a clone with proportion 0 in every sample counts —
so this table keeps every row as on source.
"""

from __future__ import annotations

import polars as pl

from tx_data import wgs_samples
from tx_data.builds._base import canonical_output_path, log
from tx_data.builds.clone_proportions import _read_cp_table
from tx_data.sources import resolve_source

TABLE = "wgs_clone_proportions"

SCHEMA = {
    "tumour_id": pl.String,
    "patient_id": pl.String,
    "wgs_tumour_label": pl.String,
    "clone": pl.String,
    "sample_id": pl.String,
    "sample_name_hash": pl.String,
    "proportion": pl.Float64,
}


def build() -> pl.DataFrame:
    wgs_samples.clear_cache()
    tumours = wgs_samples.tumours()

    frames: list[pl.DataFrame] = []
    absent: list[str] = []
    for label in tumours.get_column("wgs_tumour_label").to_list():
        path = resolve_source(TABLE, wgs_tumour_label=label)
        if not path.is_file():
            absent.append(label)
            continue
        # WES reader melts to (tumour_id, clone, sample_id, proportion); its
        # tumour_id slot carries the label here and is replaced below.
        frames.append(
            _read_cp_table(label, path).rename({"tumour_id": "wgs_tumour_label"})
        )
    log(f"read cp_table for {len(frames)}/{tumours.height} WGS tumours")
    if absent:
        log(f"warning: no cp_table.csv for {len(absent)} tumour(s): {absent}")

    if not frames:
        df = pl.DataFrame(schema=SCHEMA)
    else:
        df = pl.concat(frames, how="vertical")
        df = df.with_columns(
            pl.struct("wgs_tumour_label", "sample_id")
            .map_elements(
                lambda r: wgs_samples.match_sample(r["wgs_tumour_label"], r["sample_id"]),
                return_dtype=pl.String,
            )
            .alias("sample_name_hash")
        )
        unmatched = df.filter(pl.col("sample_name_hash").is_null()).select(
            "wgs_tumour_label", "sample_id").unique()
        if unmatched.height:
            log(
                f"warning: {unmatched.height} cp_table sample column(s) not matched to the "
                f"inventory; first: {unmatched.head(5).rows()}"
            )
        df = df.join(tumours, on="wgs_tumour_label", how="left").select(list(SCHEMA))

    df = df.sort("tumour_id", "sample_id", "clone")
    out = canonical_output_path(TABLE)
    df.write_parquet(out)
    log(f"wrote {out} ({df.height} rows × {df.width} cols)")
    return df


if __name__ == "__main__":
    build()
