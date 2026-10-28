"""Entry point: build canonical wgs_tree_clones Parquet (needs wgs_battenberg_inventory built)."""
from tx_data.builds.wgs_tree_clones import build

if __name__ == "__main__":
    build()
