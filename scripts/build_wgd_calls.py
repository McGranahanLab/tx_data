"""Entry point: build canonical wgd_calls (WES WGD calls from WGD_pipeline) Parquet.

  python scripts/build_wgd_calls.py [--source-dir <WGD_pipeline WES run, default output/wes>/wgd/05_table]

Without --source-dir: $TX_WES_WGD_CALLS_DIR, else the _sources.yml (mock) location.
"""
import argparse
from pathlib import Path

from tx_data.builds.wgd_calls import build

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--source-dir", type=Path)
    build(ap.parse_args().source_dir)
