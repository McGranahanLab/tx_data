"""Entry point: build canonical alpaca_trees Parquet."""
from tx_data.builds.alpaca_trees import build

if __name__ == "__main__":
    build()
