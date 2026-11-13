"""Generate the alpaca_corrected / alpaca_wgd_ratios_corrected mocks.

Runs darwins_analysis `alpaca_post_hoc_correction` (apc/correct.py) on the built mock tables
(`alpaca`, `clone_proportions`, `alpaca_wgd_ratios`, `wgd_calls`) and writes its output where the
release run is PROPOSED to put it (catalog/alpaca_corrected.yml):

    ALPACA/output/post_hoc_corrected/alpaca_corrected.parquet
    ALPACA/output/post_hoc_corrected/wgd_ratios_corrected.tsv

    APC_SRC=/path/to/darwins_analysis/alpaca_post_hoc_correction \
        python scripts/generate_mock_alpaca_corrected.py [--data-dir data]

Re-run after the ALPACA, clone-proportion or wgd_calls mocks change. Deterministic. With the
placeholder wgd_calls no mock clone meets the correction rule, so the corrected CN equals ALPACA's.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import polars as pl

REPO_ROOT = Path(__file__).resolve().parents[1]
OUT = (REPO_ROOT / "nemo_mock" / "project" / "proj-tracerx-lung" / "tctProjects" / "lungTx" / "Tx842"
       / "release" / "tracerx_842" / "ALPACA" / "output" / "post_hoc_corrected")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-dir", type=Path, default=REPO_ROOT / "data")
    args = ap.parse_args()
    src = Path(os.environ.get("APC_SRC", REPO_ROOT.parent / "darwins_analysis" / "alpaca_post_hoc_correction"))
    sys.path.insert(0, str(src))
    from apc import correct  # noqa: PLC0415

    d = args.data_dir
    res = correct.run(pl.read_parquet(d / "alpaca.parquet"), pl.read_parquet(d / "clone_proportions.parquet"),
                      pl.read_parquet(d / "wgd_calls.parquet"), pl.read_parquet(d / "alpaca_wgd_ratios.parquet"))
    OUT.mkdir(parents=True, exist_ok=True)
    res["alpaca_corrected"].write_parquet(OUT / "alpaca_corrected.parquet")
    res["wgd_ratios_corrected"].write_csv(OUT / "wgd_ratios_corrected.tsv", separator="\t", null_value="NA")
    print(f"wrote {OUT} ({int(res['alpaca_corrected']['corrected'].sum())} corrected rows, "
          f"{res['corrections'].height} clones)")


if __name__ == "__main__":
    main()
