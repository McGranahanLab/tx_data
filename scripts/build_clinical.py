"""Entry point: build canonical clinical Parquets (per-patient + per-tumour)."""
from tx_data.builds.clinical import build

if __name__ == "__main__":
    build()
