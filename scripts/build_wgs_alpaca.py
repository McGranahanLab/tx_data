"""Entry point: build canonical wgs_alpaca Parquet (needs wgs_battenberg_inventory built)."""
from tx_data.builds.wgs_alpaca import build

if __name__ == "__main__":
    build()
