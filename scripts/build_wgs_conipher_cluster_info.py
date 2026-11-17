"""Entry point: build canonical wgs_conipher_cluster_info Parquet (needs wgs_battenberg_inventory built)."""
from tx_data.builds.wgs_conipher_cluster_info import build

if __name__ == "__main__":
    build()
