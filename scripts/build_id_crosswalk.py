"""Entry point: build canonical id_crosswalk Parquet."""
from tx_data.builds.id_crosswalk import build

if __name__ == "__main__":
    build()
