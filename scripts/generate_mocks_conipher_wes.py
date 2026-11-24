"""Write the WES CONIPHER mutation mocks into nemo_mock's tx842 release tree (conipher_mutations, 2026-11-24).

Per tumour of the WES cohort manifest (nemo_mock/_cohort_manifest.json):
  release_tx842/<patient>/mutation_clusters/tumour_<n>/<patient>_<n>.SCoutput.CLEAN.tsv
  release_tx842/<patient>/mutation_trees/tumour_<n>/{treeTable.tsv, allTrees.txt, clusterInfo.txt}
                                                                            (tumours with a tree only)

Clusters are the manifest's tree clones plus its non-tree clusters (left out of treeTable, NA in
clusterInfo). allTrees.txt: the manifest tree first; LTX0004 gets a second, alternative tree. SAMPLE is
written R-style (`LTX0001_SU_T1.R1..<hash>`) and mutation_id without a patient prefix, as on Nemo, so the
build's conversion to the dash form is exercised. Truncal mutations sit at 2 copies on 2|2 in tumours with a
clonal GD. Deterministic.

  python scripts/generate_mocks_conipher_wes.py
"""

from __future__ import annotations

import json
import random

from tx_data.paths import DEFAULT_MOCK_ROOT

MANIFEST = DEFAULT_MOCK_ROOT / "_cohort_manifest.json"
RELEASE = DEFAULT_MOCK_ROOT / "project/proj-tracerx-lung/tracerx/_RELEASE/release_tx842"
N_MUTS_PER_CLUSTER = 20
ALT_TREE = {"LTX0004"}                    # allTrees.txt with 2 trees (n_trees = 2)
COLS = ("CASE_ID", "SAMPLE", "CHR", "POS", "REF", "ALT", "REF_COUNT", "VAR_COUNT", "DEPTH", "CLUSTER",
        "CCF_PHYLO", "CCF_OBS", "MUT_COPY", "COPY_NUMBER_A", "COPY_NUMBER_B", "ACF", "PLOIDY", "mutation_id")


def r_style(sample_name_hash: str) -> str:
    """`LTX0001_SU_T1-R1--abc` -> `LTX0001_SU_T1.R1..abc` (R make.names)."""
    return sample_name_hash.replace("-", ".")


def main() -> None:
    rng = random.Random(20261124)
    tumours = json.loads(MANIFEST.read_text())["tumours"]
    for t in tumours:
        patient, n = t["patient_id"], t["tumour_ordinal"]
        case = f"{patient}_{n}"
        tree = [int(c["name"][5:]) for c in t["clones"]] if t["has_tree"] else []
        clusters = (tree or [1]) + list(t["non_tree_clusters"])
        clonal_gd = bool(t["clones"]) and t["clones"][0]["gds_at_clone"] > 0
        samples = [r_style(s["sample_name_hash"]) for s in t["samples"]]
        scout, tt = ["\t".join(COLS)], ["\t".join(("mutation_id", "SAMPLE", "treeCLUSTER", "CCF"))]
        for cl in clusters:
            for _ in range(N_MUTS_PER_CLUSTER):
                chrom, pos = str(rng.randint(1, 22)), rng.randint(1_000_000, 50_000_000)
                ref, alt = rng.sample("ACGT", 2)
                mid = f":{chrom}:{pos}:{ref}:{alt}"   # real form: no patient prefix
                for s in samples:
                    cn = 2 if clonal_gd else 1
                    mut_copy = 2 if clonal_gd and cl == 1 else 1
                    ccf = 1.0 if cl == 1 else round(rng.uniform(0.2, 0.8), 2)
                    depth = rng.randint(60, 200)
                    var = round(depth * 0.25 * ccf)
                    scout.append("\t".join(map(str, (case, s, chrom, pos, ref, alt, depth - var, var, depth, cl, ccf,
                                                     ccf, mut_copy, cn, cn, 0.5, 2.0 * cn, mid))))
                    if cl in tree:
                        tt.append(f"{mid}\t{s}\t{cl}\t{round(100 * ccf, 1)}")
        mc = RELEASE / patient / "mutation_clusters" / f"tumour_{n}"
        mc.mkdir(parents=True, exist_ok=True)
        (mc / f"{case}.SCoutput.CLEAN.tsv").write_text("\n".join(scout) + "\n")
        if tree:
            mt = RELEASE / patient / "mutation_trees" / f"tumour_{n}"
            mt.mkdir(parents=True, exist_ok=True)
            (mt / "treeTable.tsv").write_text("\n".join(tt) + "\n")
            edges = [(int(c["parent"][5:]), int(c["name"][5:])) for c in t["clones"] if c["parent"] != "diploid"]
            trees = [edges or [(1, 1)]]
            if patient in ALT_TREE and edges[-1][0] != 1:
                trees.append(edges[:-1] + [(1, edges[-1][1])])     # last clone re-attached to the trunk
            (mt / "allTrees.txt").write_text("".join(
                f"# tree {k}\n" + "".join(f"{a}\t{b}\n" for a, b in tr) for k, tr in enumerate(trees, 1)))
            info = ["\t".join(("clusterID", "SAMPLE", "truncal", "treeClust", "meanCCF"))]
            for cl in clusters:
                for s in samples:
                    off = cl not in tree
                    ccf = "NA" if off else (98.0 if cl == 1 else round(rng.uniform(20, 80), 1))
                    info.append("\t".join(map(str, (cl, s, "TRUE" if cl == 1 else "FALSE", "NA" if off else "TRUE",
                                                     ccf))))
            (mt / "clusterInfo.txt").write_text("\n".join(info) + "\n")
    print(f"wrote WES CONIPHER mutation mocks for {len(tumours)} tumours under {RELEASE}")


if __name__ == "__main__":
    main()
