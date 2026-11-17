"""Entry point: build canonical wgs_alpaca_wgd_ratios Parquet (needs wgs_battenberg_inventory built)."""
from tx_data.builds.wgs_alpaca_wgd_ratios import build

if __name__ == "__main__":
    build()
