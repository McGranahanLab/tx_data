"""Entry point: build canonical seeding_table Parquet.

Depends on id_crosswalk when the source tumour_id family (LTX vs CRUK)
differs from the target ``id_convention``.
"""
from tx_data.builds.seeding_table import build

if __name__ == "__main__":
    build()
