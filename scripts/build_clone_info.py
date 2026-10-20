"""Entry point: build canonical clone_info Parquet.

Depends on id_crosswalk (patient-level CRUK ↔ LTX). If the crosswalk
parquet is missing this build will fail with a clear message.
"""
from tx_data.builds.clone_info import build

if __name__ == "__main__":
    build()
