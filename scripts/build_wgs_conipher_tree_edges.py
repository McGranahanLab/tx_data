"""Entry point: build canonical wgs_conipher_tree_edges Parquet (needs wgs_battenberg_inventory built)."""
from tx_data.builds.wgs_conipher_tree_edges import build

if __name__ == "__main__":
    build()
