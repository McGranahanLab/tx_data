"""Entry point: build canonical wgd_calls Parquet."""
from tx_data.builds.wgd_calls import build

if __name__ == "__main__":
    build()
