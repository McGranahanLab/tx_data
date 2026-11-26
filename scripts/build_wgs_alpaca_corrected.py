"""Entry point: build canonical wgs_alpaca_corrected + wgs_alpaca_wgd_ratios_corrected Parquet.

  python scripts/build_wgs_alpaca_corrected.py [--source-dir <WGD_pipeline WGS run>/correction]

Without --source-dir: $TX_WGS_ALPACA_CORRECTED_DIR, else the _sources.yml (mock) location.
"""
import argparse
from pathlib import Path

from tx_data.builds.alpaca_corrected import build_wgs

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--source-dir", type=Path)
    build_wgs(ap.parse_args().source_dir)
