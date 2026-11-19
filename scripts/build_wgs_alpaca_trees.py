"""Entry point: build canonical wgs_alpaca_trees Parquet."""
from tx_data.builds.wgs_alpaca_trees import build

if __name__ == "__main__":
    build()
