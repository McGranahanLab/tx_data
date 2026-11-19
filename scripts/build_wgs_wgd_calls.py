"""Entry point: build canonical wgs_wgd_calls Parquet.

  python scripts/build_wgs_wgd_calls.py [--source-dir <WGD_WGS_pipeline run>/wgd/05_table]

Without --source-dir: $TX_WGS_WGD_CALLS_DIR, else the _sources.yml (mock) location.
"""
import argparse
from pathlib import Path

from tx_data.builds.wgs_wgd_calls import build

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--source-dir", type=Path)
    build(ap.parse_args().source_dir)
