"""Entry point: build canonical wgs_clone_proportions Parquet (needs wgs_battenberg_inventory built)."""
from tx_data.builds.wgs_clone_proportions import build

if __name__ == "__main__":
    build()
