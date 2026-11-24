"""Entry point: build canonical alpaca (WES ALPACA clone CN) Parquet."""
from tx_data.builds.alpaca import build

if __name__ == "__main__":
    build()
