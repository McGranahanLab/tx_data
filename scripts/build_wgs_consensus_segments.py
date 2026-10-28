"""Entry point: build canonical wgs_consensus_segments Parquet (needs wgs_battenberg_inventory built)."""
from tx_data.builds.wgs_consensus_segments import build

if __name__ == "__main__":
    build()
