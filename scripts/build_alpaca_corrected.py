"""Entry point: build canonical alpaca_corrected + alpaca_wgd_ratios_corrected Parquet.

  python scripts/build_alpaca_corrected.py [--source-dir <correction --out dir>]

Without --source-dir: $TX_ALPACA_CORRECTED_DIR, else the _sources.yml (mock) location.
"""
import argparse
from pathlib import Path

from tx_data.builds.alpaca_corrected import build

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--source-dir", type=Path)
    build(ap.parse_args().source_dir)
