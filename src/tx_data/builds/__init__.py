"""Per-table build modules: source (TSV/etc.) → canonical Parquet."""

from tx_data.builds import (
    alpaca,
    alphamissense,
    clinical,
    clone_info,
    clone_proportions,
    conipher_trees,
    driver_list,
    hg19_refgene,
    id_crosswalk,
    kallisto,
    muttable,
    purity_ploidy,
    seeding_table,
    wgd_calls,
)

__all__ = [
    "muttable",
    "wgd_calls",
    "alphamissense",
    "clinical",
    "driver_list",
    "alpaca",
    "kallisto",
    "clone_proportions",
    "id_crosswalk",
    "clone_info",
    "seeding_table",
    "hg19_refgene",
    "conipher_trees",
    "purity_ploidy",
]

BUILDERS = {
    "muttable": muttable.build,
    "wgd_calls": wgd_calls.build,
    "alphamissense": alphamissense.build,
    "clinical": clinical.build,
    "driver_list": driver_list.build,
    "alpaca": alpaca.build,
    "kallisto": kallisto.build,
    "clone_proportions": clone_proportions.build,
    "id_crosswalk": id_crosswalk.build,
    "clone_info": clone_info.build,
    "seeding_table": seeding_table.build,
    "hg19_refgene": hg19_refgene.build,
    "conipher_trees": conipher_trees.build,
    "purity_ploidy": purity_ploidy.build,
}
