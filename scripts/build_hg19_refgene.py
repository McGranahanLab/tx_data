"""Entry point: build canonical hg19_refgene Parquet from UCSC files."""
from tx_data.builds.hg19_refgene import build

if __name__ == "__main__":
    build()
