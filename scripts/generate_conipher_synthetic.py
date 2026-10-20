"""Generate synthetic CONIPHER benchmark cases into ``synthetic/conipher/``.

Unlike ``generate_mocks.py`` (which mirrors /Nemo file layouts with values that
are only schema-plausible), these cases are simulated *from* a known clone tree,
so read counts carry the phylogenetic signal CONIPHER is meant to recover:

    expected VAF = purity * CCF * multiplicity / (purity * (cn_a + cn_b) + 2 * (1 - purity))
    var_count    ~ Binomial(depth, expected VAF)

Each case directory holds:

  - ``input.tsv``            CONIPHER input table (same columns as CONIPHER's
                             ``inst/extdata/input_table.tsv``). CLUSTER is the true
                             clone; CCF_OBS / CCF_PHYLO / MUT_COPY are estimated
                             from the simulated counts, so the table also works for
                             tree-building-only runs.
  - ``truth_tree.tsv``       true clone → parent edges, mutation counts, per-region CCF
  - ``truth_mutations.tsv``  true clone and multiplicity per mutation

Stdlib only, fixed seed: re-running produces byte-identical files.
These cases are fully synthetic; they are not derived from any patient data.
"""

from __future__ import annotations

import json
import random
from dataclasses import dataclass
from pathlib import Path

from tx_data.paths import REPO_ROOT

OUT_ROOT = REPO_ROOT / "synthetic" / "conipher"
SEED = 20261020
GENERATOR_VERSION = 1

# Chromosomes carrying a clonal single-copy gain (cn_a=2, cn_b=1) in cases that
# set `gained_fraction`. All other mutations sit on diploid segments.
GAINED_CHROMS = ("5", "8")
ALL_CHROMS = tuple(str(c) for c in range(1, 23))
OTHER_CHROMS = tuple(c for c in ALL_CHROMS if c not in GAINED_CHROMS)
BASES = "ACGT"

INPUT_COLUMNS = [
    "CASE_ID", "SAMPLE", "CHR", "POS", "REF", "ALT",
    "REF_COUNT", "VAR_COUNT", "DEPTH",
    "CLUSTER", "CCF_PHYLO", "CCF_OBS", "MUT_COPY",
    "COPY_NUMBER_A", "COPY_NUMBER_B", "ACF", "PLOIDY",
]  # fmt: skip


@dataclass(frozen=True)
class Clone:
    name: str
    parent: str  # "diploid" for the trunk
    n_mutations: int
    ccf: tuple[float, ...]  # true cancer cell fraction, one per region


@dataclass(frozen=True)
class Case:
    case_id: str
    description: str
    purity: tuple[float, ...]  # one per region
    clones: tuple[Clone, ...]
    mean_depth: int = 150
    gained_fraction: float = 0.0  # fraction of trunk mutations on gained segments

    @property
    def n_regions(self) -> int:
        return len(self.purity)

    @property
    def ploidy(self) -> float:
        # Gained chromosomes 5 + 8 are ~12% of the genome; only matters when present.
        return 2.12 if self.gained_fraction > 0 else 2.0


CASES = (
    Case(
        "SYN001",
        "Linear tree, 2 regions, diploid.",
        purity=(0.7, 0.55),
        clones=(
            Clone("clone1", "diploid", 30, (1.0, 1.0)),
            Clone("clone2", "clone1", 20, (0.6, 0.4)),
            Clone("clone3", "clone2", 20, (0.3, 0.15)),
        ),
    ),
    Case(
        "SYN002",
        "Branching tree, 3 regions, diploid; clone4 is private to R1.",
        purity=(0.8, 0.6, 0.5),
        clones=(
            Clone("clone1", "diploid", 30, (1.0, 1.0, 1.0)),
            Clone("clone2", "clone1", 20, (0.7, 0.5, 0.0)),
            Clone("clone3", "clone1", 20, (0.2, 0.4, 0.9)),
            Clone("clone4", "clone2", 15, (0.6, 0.0, 0.0)),
        ),
    ),
    Case(
        "SYN003",
        "Branching tree, 2 regions, clonal gain on chr5/chr8, low purity in R1.",
        purity=(0.35, 0.8),
        gained_fraction=0.5,
        clones=(
            Clone("clone1", "diploid", 40, (1.0, 1.0)),
            Clone("clone2", "clone1", 20, (0.8, 0.3)),
            Clone("clone3", "clone1", 20, (0.0, 0.5)),
        ),
    ),
)


def validate(case: Case) -> None:
    """Fail loudly if a case violates the tree's sum rule or shape."""
    names = {c.name for c in case.clones}
    roots = [c for c in case.clones if c.parent == "diploid"]
    assert len(roots) == 1, f"{case.case_id}: need exactly one trunk"
    assert all(ccf == 1.0 for ccf in roots[0].ccf), f"{case.case_id}: trunk CCF must be 1"
    for clone in case.clones:
        assert len(clone.ccf) == case.n_regions, f"{case.case_id}/{clone.name}: CCF length"
        assert clone.parent == "diploid" or clone.parent in names
    for parent in case.clones:
        children = [c for c in case.clones if c.parent == parent.name]
        for r in range(case.n_regions):
            total = sum(c.ccf[r] for c in children)
            assert total <= parent.ccf[r] + 1e-9, (
                f"{case.case_id}: children of {parent.name} exceed its CCF in R{r + 1}"
            )


def _binomial(rng: random.Random, n: int, p: float) -> int:
    return sum(rng.random() < p for _ in range(n))


def _positions(rng: random.Random, chroms: tuple[str, ...], n: int, used: set) -> list:
    """Draw `n` unused sites. Chromosomes without a mutation yet are filled first:
    CONIPHER's plotting (`fun.chrom.length`) fails if any of chr1..max has none."""
    out = []
    while len(out) < n:
        covered = {chrom for chrom, _ in used}
        empty = [c for c in chroms if c not in covered]
        key = (empty[0] if empty else rng.choice(chroms), rng.randint(1_000_000, 100_000_000))
        if key not in used:
            used.add(key)
            out.append(key)
    return out


def simulate(case: Case, rng: random.Random) -> tuple[list[dict], list[dict]]:
    """Return (input rows, truth mutation rows) for one case."""
    validate(case)
    used: set = set()
    mutations = []
    for idx, clone in enumerate(case.clones, start=1):
        n_gained = round(clone.n_mutations * case.gained_fraction) if idx == 1 else 0
        diploid_chroms = OTHER_CHROMS if case.gained_fraction > 0 else ALL_CHROMS
        sites = _positions(rng, GAINED_CHROMS, n_gained, used)
        sites += _positions(rng, diploid_chroms, clone.n_mutations - n_gained, used)
        for chrom, pos in sites:
            gained = chrom in GAINED_CHROMS
            ref = rng.choice(BASES)
            alt = rng.choice(BASES.replace(ref, ""))
            # Trunk mutations on a gained segment that arose before the gain sit on 2 copies.
            multiplicity = 2 if gained and rng.random() < 0.5 else 1
            mutations.append(
                {"clone": clone, "cluster": idx, "chr": chrom, "pos": pos, "ref": ref,
                 "alt": alt, "cn_a": 2 if gained else 1, "cn_b": 1, "multiplicity": multiplicity}
            )  # fmt: skip
    mutations.sort(key=lambda m: (int(m["chr"]), m["pos"]))

    rows, truth = [], []
    for m in mutations:
        truth.append(
            {"CASE_ID": case.case_id, "CHR": m["chr"], "POS": m["pos"], "REF": m["ref"],
             "ALT": m["alt"], "CLONE": m["clone"].name, "CLUSTER": m["cluster"],
             "MULTIPLICITY": m["multiplicity"]}
        )  # fmt: skip
        cn_total = m["cn_a"] + m["cn_b"]
        for r in range(case.n_regions):
            purity = case.purity[r]
            denom = purity * cn_total + 2 * (1 - purity)
            vaf_expected = purity * m["clone"].ccf[r] * m["multiplicity"] / denom
            depth = max(30, round(rng.gauss(case.mean_depth, case.mean_depth * 0.25)))
            var_count = _binomial(rng, depth, vaf_expected)
            vaf = var_count / depth
            mut_copy = vaf * denom / purity
            ccf_obs = mut_copy / m["multiplicity"]
            rows.append(
                {
                    "CASE_ID": case.case_id, "SAMPLE": f"{case.case_id}_R{r + 1}",
                    "CHR": m["chr"], "POS": m["pos"], "REF": m["ref"], "ALT": m["alt"],
                    "REF_COUNT": depth - var_count, "VAR_COUNT": var_count, "DEPTH": depth,
                    "CLUSTER": m["cluster"], "CCF_PHYLO": round(min(1.0, ccf_obs), 4),
                    "CCF_OBS": round(ccf_obs, 4), "MUT_COPY": round(mut_copy, 4),
                    "COPY_NUMBER_A": m["cn_a"], "COPY_NUMBER_B": m["cn_b"],
                    "ACF": purity, "PLOIDY": case.ploidy,
                }
            )  # fmt: skip
    rows.sort(key=lambda row: (row["SAMPLE"], int(row["CHR"]), row["POS"]))
    return rows, truth


def _write_tsv(path: Path, rows: list[dict], columns: list[str]) -> None:
    lines = ["\t".join(columns)]
    lines += ["\t".join(str(row[c]) for c in columns) for row in rows]
    path.write_text("\n".join(lines) + "\n")


def write_case(case: Case, out_root: Path = OUT_ROOT) -> Path:
    # One RNG per case, seeded from the case ID, so adding or editing a case never
    # changes the others.
    rng = random.Random(f"{SEED}:{case.case_id}")
    rows, truth = simulate(case, rng)
    case_dir = out_root / case.case_id
    case_dir.mkdir(parents=True, exist_ok=True)
    _write_tsv(case_dir / "input.tsv", rows, INPUT_COLUMNS)
    _write_tsv(case_dir / "truth_mutations.tsv", truth, list(truth[0]))
    regions = [f"CCF_R{r + 1}" for r in range(case.n_regions)]
    tree_rows = [
        {"CLONE": c.name, "CLUSTER": i, "PARENT": c.parent, "N_MUTATIONS": c.n_mutations,
         **dict(zip(regions, c.ccf))}
        for i, c in enumerate(case.clones, start=1)
    ]  # fmt: skip
    _write_tsv(case_dir / "truth_tree.tsv", tree_rows, list(tree_rows[0]))
    print(f"wrote {case_dir}  ({len(truth)} mutations x {case.n_regions} regions)")
    return case_dir


def write_manifest(out_root: Path = OUT_ROOT) -> None:
    manifest = {
        "generator": "scripts/generate_conipher_synthetic.py",
        "generator_version": GENERATOR_VERSION,
        "seed": SEED,
        "cases": {
            c.case_id: {"description": c.description, "purity": list(c.purity),
                        "ploidy": c.ploidy, "gained_fraction": c.gained_fraction}
            for c in CASES
        },
    }  # fmt: skip
    (out_root / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


def main() -> None:
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    for case in CASES:
        write_case(case)
    write_manifest()


if __name__ == "__main__":
    main()
