"""Build purity_ploidy → canonical Parquet.

Source is a flat CSV keyed by `sample_name_hash`. The build script parses
the hash into canonical `tumour_id` + `region` so the table joins directly
against every other per-tumour table while preserving the original hash
for lossless joins against muttable + kallisto.
"""

from __future__ import annotations

import polars as pl

from tx_data.builds._base import canonical_output_path, log
from tx_data.ids import apply_convention_scalar
from tx_data.sources import resolve_source

TABLE = "purity_ploidy"

_SAMPLE_RE = (
    r"^(?P<prefix>LTX|CRUK)"
    r"(?P<patient_digits>\d+)"
    r"_SU_T(?P<tumour_ordinal>\d+)"
    r"-R(?P<region>\d+)"
    r"--(?P<hash>[a-z0-9]+)$"
)


def _derive_tumour_id(sample_hash: str) -> str:
    import re

    m = re.match(_SAMPLE_RE, sample_hash)
    if m is None:
        raise ValueError(
            f"sample_name_hash does not match expected pattern "
            f"'<LTX|CRUK>#_SU_T#-R#--<hash>': {sample_hash!r}"
        )
    patient = f"{m.group('prefix')}{m.group('patient_digits')}"
    ordinal = int(m.group("tumour_ordinal"))
    return apply_convention_scalar(f"{patient}-Tumour{ordinal}")


def _derive_region(sample_hash: str) -> str:
    import re

    m = re.match(_SAMPLE_RE, sample_hash)
    if m is None:
        raise ValueError(f"sample_name_hash parse failed: {sample_hash!r}")
    return f"R{int(m.group('region'))}"


def build() -> pl.DataFrame:
    src = resolve_source(TABLE)
    log(f"read {src}")

    df = pl.read_csv(src, null_values=["NA", "", "nan"], infer_schema_length=10000)

    expected = {"sample_name_hash", "purity", "ploidy"}
    missing = expected - set(df.columns)
    if missing:
        raise ValueError(
            f"purity_ploidy source missing expected columns {sorted(missing)}; "
            f"got {df.columns!r}"
        )

    df = df.with_columns(
        [
            pl.col("sample_name_hash")
            .map_elements(_derive_tumour_id, return_dtype=pl.String)
            .alias("tumour_id"),
            pl.col("sample_name_hash")
            .map_elements(_derive_region, return_dtype=pl.String)
            .alias("region"),
            pl.col("purity").cast(pl.Float64),
            pl.col("ploidy").cast(pl.Float64),
        ]
    )

    dupes = df.group_by("sample_name_hash").len().filter(pl.col("len") > 1)
    if dupes.height:
        raise ValueError(
            f"purity_ploidy: {dupes.height} duplicate sample_name_hash rows. "
            f"First offenders: {dupes.head(5).to_dicts()}"
        )

    tr_dupes = df.group_by(["tumour_id", "region"]).len().filter(pl.col("len") > 1)
    if tr_dupes.height:
        raise ValueError(
            f"purity_ploidy: {tr_dupes.height} duplicate (tumour_id, region) rows. "
            f"First offenders: {tr_dupes.head(5).to_dicts()}"
        )

    out_df = df.select(
        ["tumour_id", "region", "sample_name_hash", "purity", "ploidy"]
    ).sort(["tumour_id", "region"])

    out = canonical_output_path(TABLE)
    out_df.write_parquet(out)
    log(f"wrote {out} ({out_df.height} rows × {out_df.width} cols)")
    return out_df


if __name__ == "__main__":
    build()
