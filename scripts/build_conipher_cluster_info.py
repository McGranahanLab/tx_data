"""Entry point: build canonical conipher_cluster_info (WES) Parquet."""
from tx_data.builds.conipher_cluster_info import build

if __name__ == "__main__":
    build()
