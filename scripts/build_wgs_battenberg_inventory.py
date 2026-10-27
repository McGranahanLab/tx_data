"""Entry point: build canonical wgs_battenberg_inventory Parquet."""
from tx_data.builds.wgs_battenberg_inventory import build

if __name__ == "__main__":
    build()
