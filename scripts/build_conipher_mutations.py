"""Entry point: build canonical conipher_mutations (WES) Parquet (sample_id needs clone_proportions built)."""
from tx_data.builds.conipher_mutations import build

if __name__ == "__main__":
    build()
