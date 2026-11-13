"""Entry point: build canonical alpaca_corrected + alpaca_wgd_ratios_corrected Parquet."""
from tx_data.builds.alpaca_corrected import build

if __name__ == "__main__":
    build()
