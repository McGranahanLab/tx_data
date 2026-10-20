"""Standalone build entrypoint for ``conipher_trees`` (RDS → two parquets).

Pattern mirrors ``build_wgd_calls.py`` etc. Shells out to R under the hood.
"""

from tx_data.builds.conipher_trees import build

if __name__ == "__main__":
    build()
