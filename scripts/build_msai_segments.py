"""Entry point: build canonical msai_segments (WES, from the ALPACA input) Parquet."""
from tx_data.builds.msai_segments import build

if __name__ == "__main__":
    build()
