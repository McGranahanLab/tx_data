"""Entry point: build canonical cohort_metadata Parquet."""
from tx_data.builds.cohort_metadata import build

if __name__ == "__main__":
    build()
