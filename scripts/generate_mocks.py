"""Generate synthetic mock cohort into ``nemo_mock/``.

Produces four coherent tables spanning every meaningful edge case:

  - 10 tumours (~9 patients, one with a second tumour for multi-tumour coverage)
  - Every wgd_calls `class`: no_wgd, clonal_wgd, mut_supported, ploidy_only
  - A `sub_class` split (parallel / sequential / single subclonal GD)
  - Every status: mostly `resolved`, one `unresolved`, one `needs_follow_up`
  - Both High and Low quality
  - A subclonal-WGD tumour where the WGD event is on an internal clone (so
    tree traversal has to run) plus a parallel-WGD tumour with two independent
    doublings
  - Every gene set: essential-only, driver-only, non-essential (drawn from
    fixed lists so we can assert on them in tests)
  - Every variant class: SNV missense, SNV synonymous, SNV nonsense, indel

Files written under ``nemo_mock/``:
  - ``tx842_mutation_table.tsv.gz``  (long format, one row per mutation×sample)
  - ``WGD/release_tx842/WGD_calls.tsv``
  - ``tx842_clinical.tsv``
  - ``alt/cn_ccf_alphamissense/output/annotated_muttables/tx842/<tumour_muttable_form>_muttable_alpha.tsv``

Deterministic (fixed seed) so mocks re-generate identically.
"""

from __future__ import annotations

import gzip
import random
from dataclasses import dataclass, field
from pathlib import Path

from tx_data.paths import REPO_ROOT

MOCK_ROOT = REPO_ROOT / "nemo_mock"
AM_ROOT = MOCK_ROOT / "alt" / "cn_ccf_alphamissense" / "output" / "annotated_muttables" / "tx842"
WGD_DIR = MOCK_ROOT / "WGD" / "release_tx842"
ALPACA_DIR = MOCK_ROOT / "ALPACA" / "output" / "cohort_results"
CLONE_PROP_ROOT = MOCK_ROOT / "ALPACA" / "input"
KALLISTO_DIR = MOCK_ROOT / "alt" / "tx842_mets_rnaseq" / "_aggregate" / "kallisto_mat"

RNG = random.Random(20260831)

# --- Fixed gene panel (chr:pos ranges chosen to look hg19-ish) ---------------
# Each gene gets a fake 1kb coding region on chr10 / chr17 etc. Keep it tiny so
# mock TSVs stay small.
GENES = [
    # (symbol, chr, start_pos, ensembl_id, category)
    ("TP53",     "17",  7577000, "ENSG00000141510", "driver"),
    ("KRAS",     "12", 25398000, "ENSG00000133703", "driver"),
    ("PDCD11",   "10", 105200000, "ENSG00000148843", "nonessential"),
    ("GSTO2",    "10", 106039000, "ENSG00000065621", "nonessential"),
    ("ARHGAP19", "10",  99003000, "ENSG00000213390", "nonessential"),
    ("EGFR",     "7",  55086000, "ENSG00000146648", "driver"),
    ("MYC",      "8", 128748000, "ENSG00000136997", "driver"),
    ("POLR2A",   "17",  7389000, "ENSG00000181222", "essential"),
    ("RPL5",     "1",  93297000, "ENSG00000122406", "essential"),
    ("EIF3A",    "10", 120790000, "ENSG00000107581", "essential"),
    ("RPS3",     "11", 75111000, "ENSG00000149273", "essential"),
    ("MRPS2",    "9", 139036000, "ENSG00000122140", "essential"),
    ("HSPE1",    "2", 198365000, "ENSG00000115541", "essential"),
    ("SNRPD1",   "18", 19195000, "ENSG00000167088", "essential"),
    ("PSMA1",    "11", 14528000, "ENSG00000129084", "essential"),
]

CODON_TABLE = {
    "missense":   ("Ctt", "Gtt", "L", "V"),
    "synonymous": ("Ctt", "Ctc", "L", "L"),
    "nonsense":   ("Cga", "Tga", "R", "*"),
}


# --- Tumour cohort definition -------------------------------------------------

@dataclass
class Clone:
    name: str
    parent: str
    gds_at_clone: int = 0


@dataclass
class TumourSpec:
    patient_id: str
    tumour_ordinal: int
    histology_group: str
    class_: str
    sub_class: str | None
    status: str
    quality: str
    n_samples: int
    clones: list[Clone]
    total_gds: int = field(init=False)

    def __post_init__(self) -> None:
        self.total_gds = sum(c.gds_at_clone for c in self.clones)

    @property
    def tumour_id_muttable(self) -> str:
        return f"{self.patient_id}_tumour{self.tumour_ordinal}"

    @property
    def tumour_id_canonical(self) -> str:
        return f"{self.patient_id}-Tumour{self.tumour_ordinal}"


def _linear_tree(n_clones: int) -> list[Clone]:
    """clone1 (root) → clone2 → clone3 → ..."""
    clones = [Clone("clone1", "diploid")]
    for i in range(2, n_clones + 1):
        clones.append(Clone(f"clone{i}", f"clone{i - 1}"))
    return clones


def _branching_tree() -> list[Clone]:
    """
        clone1
       /      \\
    clone2    clone3
       |         |
    clone4    clone5
    """
    return [
        Clone("clone1", "diploid"),
        Clone("clone2", "clone1"),
        Clone("clone3", "clone1"),
        Clone("clone4", "clone2"),
        Clone("clone5", "clone3"),
    ]


def cohort() -> list[TumourSpec]:
    specs: list[TumourSpec] = []

    # 1. No-WGD LUAD.
    specs.append(TumourSpec("LTX0001", 1, "LUAD", "no_wgd", None,
                            "resolved", "High", 3, _linear_tree(3)))

    # 2. Clonal-WGD LUAD (single WGD at MRCA).
    tree = _linear_tree(4)
    tree[0].gds_at_clone = 1
    specs.append(TumourSpec("LTX0002", 1, "LUAD", "clonal_wgd", None,
                            "resolved", "High", 3, tree))

    # 3. Clonal-WGD LUSC.
    tree = _linear_tree(3)
    tree[0].gds_at_clone = 1
    specs.append(TumourSpec("LTX0003", 1, "LUSC", "clonal_wgd", None,
                            "resolved", "Low", 2, tree))

    # 4. Subclonal WGD (mut_supported), single event on an internal clone.
    tree = _branching_tree()
    tree[3].gds_at_clone = 1  # clone4 acquires the WGD
    specs.append(TumourSpec("LTX0004", 1, "LUAD", "mut_supported",
                            "single subclonal GD", "resolved", "High", 4, tree))

    # 5. Parallel subclonal WGDs (mut_supported, parallel only).
    tree = _branching_tree()
    tree[3].gds_at_clone = 1  # clone4
    tree[4].gds_at_clone = 1  # clone5 — independent branch
    specs.append(TumourSpec("LTX0005", 1, "LUAD", "mut_supported",
                            "parallel only", "resolved", "High", 4, tree))

    # 6. Sequential subclonal WGDs (mut_supported, sequential only).
    tree = _linear_tree(4)
    tree[1].gds_at_clone = 1  # clone2
    tree[2].gds_at_clone = 1  # clone3 — sequential on the same path
    specs.append(TumourSpec("LTX0006", 1, "LUAD", "mut_supported",
                            "sequential only", "resolved", "High", 3, tree))

    # 7. ploidy_only subclonal WGD (weaker call).
    tree = _branching_tree()
    tree[3].gds_at_clone = 1
    specs.append(TumourSpec("LTX0007", 1, "LUAD", "ploidy_only",
                            "single subclonal GD", "resolved", "Low", 3, tree))

    # 8. Multi-tumour patient: Tumour1 is clonal_wgd, Tumour2 is no_wgd.
    tree = _linear_tree(3)
    tree[0].gds_at_clone = 1
    specs.append(TumourSpec("LTX0008", 1, "LUAD", "clonal_wgd", None,
                            "resolved", "High", 2, tree))
    specs.append(TumourSpec("LTX0008", 2, "LUSC", "no_wgd", None,
                            "resolved", "High", 2, _linear_tree(2)))

    # 9. Unresolved tumour — should be filtered out downstream.
    specs.append(TumourSpec("LTX0009", 1, "other", "no_wgd", None,
                            "unresolved", "Low", 2, _linear_tree(2)))

    # 10. needs_follow_up tumour — also filtered out.
    tree = _linear_tree(3)
    tree[0].gds_at_clone = 2  # double-clonal at MRCA (exercises number_of_gds > 1)
    specs.append(TumourSpec("LTX0010", 1, "LUAD", "clonal_wgd", None,
                            "needs_follow_up", "High", 2, tree))

    return specs


# --- Muttable ----------------------------------------------------------------

def _assign_mutations(spec: TumourSpec) -> list[dict]:
    """Return a list of mutation-clone assignments for this tumour.

    Every non-root clone gets 2–4 mutations across the fixed gene panel.
    Root gets more (trunk mutations).
    """
    muts: list[dict] = []
    for i, clone in enumerate(spec.clones):
        n = 6 if clone.parent == "diploid" else RNG.randint(2, 4)
        for _ in range(n):
            g = RNG.choice(GENES)
            # Pick a variant class.
            klass = RNG.choices(
                ["missense", "synonymous", "nonsense", "indel"],
                weights=[0.55, 0.25, 0.10, 0.10],
            )[0]
            pos = g[2] + RNG.randint(0, 990)
            if klass == "indel":
                # Simple 1bp deletion.
                start = pos
                stop = pos
                pos_left = pos - 1
                ref = RNG.choice(["A", "C", "G", "T"])
                var = "-"
                is_snv, is_indel, is_mnv = False, True, False
                func = "exonic" if g[4] != "nonessential" else "ncRNA_intronic"
                aachange = "-"
            else:
                start = stop = pos_left = pos
                ref_c, alt_c, ref_aa, alt_aa = CODON_TABLE[klass]
                ref = ref_c[0].upper()
                var = alt_c[0].upper()
                is_snv, is_indel, is_mnv = True, False, False
                func = "exonic"
                aachange = f"{g[0]}:NM_XXX:exon1:c.{ref}1{var}:p.{ref_aa}1{alt_aa}"
            muts.append(
                dict(
                    clone_name=clone.name,
                    chr=g[1], start=start, stop=stop, pos_left=pos_left,
                    ref=ref, var=var,
                    is_snv=is_snv, is_indel=is_indel, is_mnv=is_mnv,
                    gene=g[0], gene_ensembl=g[3], func=func, aachange=aachange,
                    variant_class=klass,
                )
            )
    return muts


# Cache per-tumour mutation lists so muttable and AlphaMissense share them.
_MUTATIONS_CACHE: dict[str, list[dict]] = {}


def _mutations_for(spec: TumourSpec) -> list[dict]:
    key = spec.tumour_id_canonical
    if key not in _MUTATIONS_CACHE:
        _MUTATIONS_CACHE[key] = _assign_mutations(spec)
    return _MUTATIONS_CACHE[key]


def _write_muttable(specs: list[TumourSpec]) -> None:
    out = MOCK_ROOT / "tx842_mutation_table.tsv.gz"
    out.parent.mkdir(parents=True, exist_ok=True)

    # Full header from the catalog spec.
    header = [
        "patient_id", "patient_tumour", "chr", "start", "stop",
        "sample_name_hash", "ref", "var", "pos",
        "var_count", "ref_count", "cov", "vaf",
        "max_var_count", "min_depth", "max_vaf", "is_present",
        "gl_sample_name_hash", "is_snv", "is_indel", "is_mnv",
        "varscan", "mutect", "vaf02",
        "varscan_any", "mutect_any", "vaf02_any",
        "gene", "func_refgene", "aachange_refgene",
        "snp137", "cosmic90_coding", "cadd13gt10",
        "ljb23_sift", "ljb23_pp2hvar", "gnomad211_af",
        "CALLED", "PASS",
        "pan_cancer_driver", "lung_driver", "luad_driver", "lusc_driver",
        "oncogene", "tumor_suppressor", "cosmic_count",
        "is_deleterious", "driver_category", "is_driver_mut",
        "tumour_id", "ffpe", "ffpe_confidence", "mutation_id",
        "cn_a", "cn_b", "major_cn", "minor_cn", "mutcpn",
        "phylo_ccf", "obs_ccf", "clone_region_proportion",
        "mean_cluster_region_ccf", "cluster_region_clonality",
        "is_trunk", "timing", "mutation_cluster",
        "clean_cluster", "tree_cluster", "tree_cn_removed_cluster",
        "is_cn_fail",
        "ith_state", "ith_state_ambiguity1",
        "is_homogen", "is_artefact_sig_fail", "is_homogen_ambiguity1",
        "rna_var_count", "rna_vaf", "rna_cov",
        "is_expressed", "is_expressed_any",
    ]

    rows = []
    for spec in specs:
        muts = _mutations_for(spec)
        sample_ids = [
            f"{spec.patient_id}_SU_T{spec.tumour_ordinal}-R{i + 1}--{RNG.randbytes(6).hex()}"
            for i in range(spec.n_samples)
        ]
        gl = f"{spec.patient_id}_BS_GL--{RNG.randbytes(6).hex()}"

        for m in muts:
            cluster_idx = int(m["clone_name"].removeprefix("clone"))
            mutation_id = f"{spec.patient_id}_{spec.tumour_ordinal}:{m['chr']}:{m['pos_left']}:{m['ref']}:{m['var']}"
            is_trunk = m["clone_name"] == "clone1"
            for s_idx, sid in enumerate(sample_ids):
                # Sample-level presence: 90% present in the "primary" sample per clone,
                # 60% detected in others via rediscovery.
                is_present = RNG.random() < (0.9 if s_idx == 0 else 0.6)
                vaf = round(RNG.uniform(0.05, 0.45), 4) if is_present else 0.0
                cov = RNG.randint(60, 350)
                var_count = int(vaf * cov)
                ref_count = cov - var_count

                # Local CN — coarse but plausible: WGD-carrying clones get higher CN.
                # Use `mutation_cluster` to look up cumulative WGD on the fly.
                cn_boost = _cumulative_gds(spec, m["clone_name"])
                major_cn = 2.0 + cn_boost
                minor_cn = 1.0 + max(0, cn_boost - 1)
                mutcpn = round(RNG.uniform(1.0, major_cn + minor_cn), 3)

                is_driver = m["gene"] in {"TP53", "KRAS", "EGFR", "MYC"}

                row = [
                    spec.patient_id, spec.tumour_id_muttable, m["chr"], m["start"], m["stop"],
                    sid, m["ref"], m["var"], m["pos_left"],
                    var_count, ref_count, cov, vaf,
                    max(var_count, 40), min(cov, 120), min(0.5, vaf + 0.05), is_present,
                    gl, m["is_snv"], m["is_indel"], m["is_mnv"],
                    True, RNG.random() < 0.7, vaf > 0.02,
                    True, True, True,
                    m["gene"], m["func"], m["aachange"],
                    "", "", m["variant_class"] == "nonsense",
                    "", "", "",
                    True,  # CALLED — set to True; we'll test False elsewhere in future
                    True,  # PASS — always True (matches Tx842)
                    is_driver, is_driver, is_driver and m["gene"] != "MYC",
                    is_driver and m["gene"] == "TP53",
                    m["gene"] in {"KRAS", "EGFR", "MYC"},
                    m["gene"] == "TP53",
                    0,
                    m["variant_class"] in ("nonsense", "indel"),
                    "1" if is_driver else "",
                    is_driver,
                    spec.tumour_ordinal, False, "high_conf",
                    mutation_id,
                    major_cn - minor_cn, minor_cn, major_cn, minor_cn, mutcpn,
                    round(RNG.uniform(0.6, 1.0), 3), round(RNG.uniform(0.6, 1.0), 3),
                    0.0, round(RNG.uniform(0.7, 1.0), 2),
                    "clonal" if is_trunk else "subclonal",
                    is_trunk, "early",
                    float(cluster_idx),           # mutation_cluster (float)
                    True, True, False, False,
                    1, 1, is_trunk, False, is_trunk,
                    "", "", "", False, False,
                ]
                rows.append(row)

    with gzip.open(out, "wt", newline="\n") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join(_fmt(x) for x in r) + "\n")

    print(f"wrote {out}  ({len(rows)} rows)")


def _cumulative_gds(spec: TumourSpec, clone_name: str) -> int:
    """Sum of gds_at_clone along the parent chain to root (inclusive)."""
    idx = {c.name: c for c in spec.clones}
    total = 0
    node: str | None = clone_name
    while node and node != "diploid":
        total += idx[node].gds_at_clone
        node = idx[node].parent
    return total


def _fmt(v: object) -> str:
    if v is None:
        return "NA"
    if isinstance(v, bool):
        return "True" if v else "False"
    if isinstance(v, float):
        return f"{v:g}"
    return str(v)


# --- WGD calls ---------------------------------------------------------------

def _write_wgd_calls(specs: list[TumourSpec]) -> None:
    WGD_DIR.mkdir(parents=True, exist_ok=True)
    out = WGD_DIR / "WGD_calls.tsv"

    header = [
        "tumour_id", "clone", "parent",
        "number_of_gds_at_clone", "total_gds_per_tumour",
        "quality", "status", "class", "sub_class",
    ]

    lines = ["\t".join(header)]
    for spec in specs:
        for c in spec.clones:
            lines.append("\t".join([
                spec.tumour_id_canonical,
                c.name,
                c.parent,
                str(c.gds_at_clone),
                str(spec.total_gds),
                spec.quality,
                spec.status,
                spec.class_ if spec.status == "resolved" else "",
                spec.sub_class or "",
            ]))

    out.write_text("\n".join(lines) + "\n")
    print(f"wrote {out}")


# --- Clinical ----------------------------------------------------------------

def _write_clinical(specs: list[TumourSpec]) -> None:
    patients: dict[str, list[TumourSpec]] = {}
    for s in specs:
        patients.setdefault(s.patient_id, []).append(s)

    header = [
        "REGTrialNo", "Patient_ID", "age", "sex", "ethnicity",
        "pre_surgery_FEV1", "registration_ECOG",
        "smoking_status", "smoking_status_group",
        "Number_cigs_perday_use", "Number_years_smoking_use",
        "pack_years_calculated", "pack_years_truncated",
        "age_start_smoking_survey", "age_start_smoking_exsmoker", "age_stop_smoking_exsmoker",
        "alcohol_units", "alcohol_frequency", "alcohol_frequency_over7units",
        "alcohol_years", "alcohol_years_truncated",
        "Surgery_type", "adjuvant_treatment_given", "adjuvant_treatment_YN",
        "Centrally_reviewed_histology_per_patient",
        "Lesion1_site_central.reviewed", "histology1_central.reviewed",
        "histology1_group_central.reviewed",
        "pTStage_v8_lesion1_central.reviewed", "pNStage_lesion1_central.reviewed",
        "pTNMStage_v8_lesion1_central.reviewed", "SizePath_lesion1_central.reviewed",
        "PathPleuInv_lesion1_central.reviewed", "margin_status_lesion1_central.reviewed",
        "Lesion2_site_central.reviewed", "histology2_central.reviewed",
        "histology2_group_central.reviewed",
        "pTStage_v8_lesion2_central.reviewed", "pNStage_lesion2_central.reviewed",
        "pTNMStage_v8_lesion2_central.reviewed", "SizePath_lesion2_central.reviewed",
        "PathPleuInv_lesion2_central.reviewed", "margin_status_lesion2_central.reviewed",
        "note_histology", "note_pathology_report",
        "new_primary_cancer", "new_primary_flag",
        "second_lesion_detected_in_baseline", "CTC_notes",
        "Relapse_cat_new", "Relapse_site_text",
        "Surgical_bed_recurrence", "Ipsilateral_lung", "Mediastinum",
        "Contralateral_lung", "Pleural_nodules_effusion",
        "Liver", "Adrenal", "Bone", "Brain",
        "Extrathoracic_lymph_nodes", "Other", "Brain_only",
        "extracranial_extrathoracic_single",
        "PrevCancerAge", "Previous_cancer",
        "HomeSmoke", "WorkSmoke30Yrs", "WorkArsenic", "WorkAsbestos",
        "WorkBenzene", "WorkBisphenolA", "WorkCH2O", "WorkChemicals",
        "WorkChromHex", "WorkDiesel", "WorkDioxins", "WorkMetals",
        "WorkNatFibres", "WorkPAHs", "WorkPBDEs", "WorkRadon",
        "WorkSolvents", "WorkVinylCl",
        "PathTissueFind",
    ]

    def _blank_row(pid: str) -> list[str]:
        return ["NA"] * len(header)

    lines = ["\t".join(header)]
    for pid, tumours in patients.items():
        row = _blank_row(pid)
        row[0] = f"A_{pid}"
        row[1] = pid
        row[2] = str(RNG.randint(45, 82))
        row[3] = RNG.choice(["Male", "Female"])
        row[4] = "White- British"
        row[7] = RNG.choice(["Never", "Ex-Smoker", "Current"])
        row[8] = row[7]
        row[11] = str(RNG.randint(0, 60))
        row[12] = row[11]
        row[21] = "Lobectomy"
        row[22] = "No adjuvant"
        row[23] = "No"
        row[24] = tumours[0].histology_group

        # Lesion1 always populated.
        t1 = tumours[0]
        row[25] = "Right Upper Lobe"
        row[26] = "Invasive adenocarcinoma" if t1.histology_group == "LUAD" else "Squamous cell carcinoma"
        row[27] = t1.histology_group
        row[28] = RNG.choice(["1a", "1b", "2a", "2b", "3a"])
        row[29] = "0"

        # Lesion2 populated iff multi-tumour.
        if len(tumours) > 1:
            t2 = tumours[1]
            row[34] = "Left Lower Lobe"
            row[35] = "Squamous cell carcinoma"
            row[36] = t2.histology_group
            row[37] = "1b"

        lines.append("\t".join(row))

    out = MOCK_ROOT / "tx842_clinical.tsv"
    out.write_text("\n".join(lines) + "\n")
    print(f"wrote {out}")


# --- AlphaMissense (muttable-shaped with AM columns appended) ----------------
#
# Real Tx842 layout (confirmed 2026-09-02): per-tumour TSVs whose columns are
# muttable's columns + `alpha_missense` (float score) + `am_class`
# (benign / ambiguous / pathogenic; NA for non-missense). The AM loader only
# reads `chr, pos, var, alpha_missense, am_class` — everything else can be a
# stub, but we emit the full muttable header shape so the mock's columns
# match what a real file would carry.

def _write_alphamissense(specs: list[TumourSpec]) -> None:
    AM_ROOT.mkdir(parents=True, exist_ok=True)

    # Muttable-shaped header (subset the loader needs, plus enough padding
    # to look realistic). Order doesn't matter for the loader — it reads by
    # name — but we mirror muttable's leading columns for readability.
    header = [
        "patient_id", "patient_tumour", "chr", "start", "stop",
        "ref", "var", "pos", "gene", "alpha_missense", "am_class",
    ]

    for spec in specs:
        muts = _mutations_for(spec)
        out = AM_ROOT / f"{spec.tumour_id_muttable}_muttable_alpha.tsv"

        lines = ["\t".join(header)]
        for m in muts:
            if m["variant_class"] == "missense":
                am_score = round(RNG.uniform(0.6, 0.99) if m["gene"] in {"TP53", "KRAS", "EGFR"}
                                 else RNG.uniform(0.02, 0.35), 4)
                # Emit raw upstream labels ("Likely benign" / "Likely
                # pathogenic") so the tx_data build's normalisation is
                # exercised on mocks too.
                am_class = ("Likely pathogenic" if am_score > 0.564
                            else "ambiguous" if am_score >= 0.34
                            else "Likely benign")
                am_score_s, am_class_s = str(am_score), am_class
            else:
                # AM doesn't score non-missense; upstream writes `not_classified`.
                am_score_s, am_class_s = "NA", "not_classified"

            lines.append("\t".join([
                spec.patient_id,
                spec.tumour_id_muttable,
                m["chr"],
                str(m["start"]),
                str(m["stop"]),
                m["ref"],
                m["var"],
                str(m["pos_left"]),
                m["gene"],
                am_score_s,
                am_class_s,
            ]))
        out.write_text("\n".join(lines) + "\n")
    print(f"wrote AM files for {len(specs)} tumours under {AM_ROOT}")


# --- Driver list ------------------------------------------------------------

def _write_driver_list() -> None:
    """Emit a synthetic TX842-shaped driver gene list.

    Categorises our fixed gene panel:
      - `TP53`, `KRAS`, `EGFR` → lung_mut_driver TRUE
      - `MYC` → mut_driver TRUE but not lung-specific
      - others → not drivers
    """
    out_dir = MOCK_ROOT / "drivers"
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / "tx842_driver_gene_list.tsv"

    header = [
        "Hugo_Symbol", "gene_id", "driver_gene", "driver_role", "CGC_TIER",
        "mut_driver", "lung_mut_driver", "mut_driver_source",
        "CN_driver", "lung_CN_driver", "CNA_type", "CN_driver_source",
        "chr_hg38", "start_hg38", "end_hg38", "cyto_hg38", "strand_hg38",
        "chr_hg19", "start_hg19", "end_hg19", "cyto_hg19",
        "SYNONYMS",
    ]

    lung_mut_drivers = {"TP53", "KRAS", "EGFR"}
    pan_mut_drivers = lung_mut_drivers | {"MYC"}
    lung_cn_drivers = {"MYC", "EGFR"}
    cn_drivers = lung_cn_drivers | {"KRAS"}

    lines = ["\t".join(header)]
    for sym, chrom, start, ensembl, category in GENES:
        is_lung_mut = sym in lung_mut_drivers
        is_mut = sym in pan_mut_drivers
        is_lung_cn = sym in lung_cn_drivers
        is_cn = sym in cn_drivers
        driver = is_mut or is_cn

        role = {
            "TP53": "tumor_suppressor",
            "KRAS": "oncogene",
            "EGFR": "oncogene",
            "MYC":  "oncogene",
        }.get(sym, "NA")

        cgc_tier = "1" if sym in {"TP53", "KRAS", "EGFR", "MYC"} else "NA"
        cna_type = "amp" if is_cn and sym != "TP53" else ("NA" if not is_cn else "del")

        end = start + 990
        lines.append("\t".join([
            sym, ensembl,
            "TRUE" if driver else "FALSE",
            role,
            cgc_tier,
            "TRUE" if is_mut else "FALSE",
            "TRUE" if is_lung_mut else "FALSE",
            "CGC_v104" if is_mut else "NA",
            "TRUE" if is_cn else "FALSE",
            "TRUE" if is_lung_cn else "FALSE",
            cna_type,
            "CGC_v104" if is_cn else "NA",
            f"chr{chrom}", str(start), str(end), f"{chrom}q00.0", "+",
            f"chr{chrom}", str(start), str(end), f"{chrom}q00.0",
            "NA",
        ]))

    out.write_text("\n".join(lines) + "\n")
    print(f"wrote {out}")


# --- Entry ------------------------------------------------------------------

# --- ALPACA (per-clone allele-specific CN) ----------------------------------
#
# For each tumour, emit rows for every clone × a fixed set of segments
# covering the mock GENES panel. Segment CN is realistic-ish: baseline 1|1;
# WGD-bearing clones get diploid (2|2) or amplified; a fraction of clones
# gets targeted LOH or homdel at TP53 / CDKN2A to exercise those paths.

_ALPACA_SEGMENT_HEADERS = ("tumour_id", "segment", "clone", "pred_CN_A", "pred_CN_B")


def _tumour_segments(spec: TumourSpec) -> list[tuple[str, int, int]]:
    """Return `[(chr, start, end), …]` for this tumour.

    Segments are a small superset of the mock gene panel — one segment per
    chromosome that appears in GENES, wide enough to cover every gene on that
    chromosome plus a buffer.
    """
    chrom_ranges: dict[str, tuple[int, int]] = {}
    for _sym, chrom, start, _ens, _cat in GENES:
        end = start + 990
        lo, hi = chrom_ranges.get(chrom, (start, end))
        chrom_ranges[chrom] = (min(lo, start), max(hi, end))

    segments: list[tuple[str, int, int]] = []
    for chrom, (lo, hi) in sorted(chrom_ranges.items()):
        # Add a small buffer so segment boundaries don't sit exactly on gene starts.
        segments.append((chrom, max(1, lo - 100), hi + 100))
    return segments


def _pick_alpaca_cn(spec: TumourSpec, clone: Clone, chrom: str) -> tuple[int, int]:
    """Return (pred_CN_A, pred_CN_B) for this clone × chromosome.

    Realistic-ish mock rules:
      - Baseline 1|1 diploid.
      - WGD in this clone's lineage → 2|2 (single WGD) or 4|4 (double).
      - TP53 (chr17) in a subclonal branch of a mut_supported tumour → LOH 1|0.
      - CDKN2A stand-in (chr7, using EGFR's chr) in one clone of LTX0005 → homdel 0|0.
      - Per-clone stochastic CN drift for non-root clones on non-hotspot
        chromosomes, so downstream SCNA analyses see non-trivial subclonal
        variation instead of a flat 1|1 landscape. Seeded per
        (tumour, clone, chrom) so runs stay reproducible.
    """
    # WGD baseline.
    cum_gd = _cumulative_gds(spec, clone.name)
    a = 1 + cum_gd
    b = 1 + cum_gd

    # Targeted edge cases so tests can assert on specific tumours.
    hotspot = False
    if spec.tumour_id_canonical == "LTX0004-Tumour1" and clone.name in ("clone2", "clone4") and chrom == "17":
        # TP53 LOH on the branch containing clone4.
        a, b = a + b, 0
        hotspot = True
    if spec.tumour_id_canonical == "LTX0005-Tumour1" and clone.name == "clone5" and chrom == "12":
        # Homdel on a small chromosome-12 segment (covering the mock KRAS region).
        a, b = 0, 0
        hotspot = True
    if spec.tumour_id_canonical == "LTX0002-Tumour1" and clone.name == "clone2" and chrom == "10":
        # Post-WGD LOH: 4→2|0 (loses one allele after doubling).
        a, b = a + b, 0
        hotspot = True

    # Stochastic per-clone SCNA drift so subclonal analyses see signal.
    # Root clones (parent="diploid") stay at baseline; only descendants drift.
    # Skips hotspot overrides so those stay deterministic.
    if not hotspot and clone.parent != "diploid":
        import hashlib
        import random as _random

        # `hash()` is salted per interpreter via PYTHONHASHSEED, so use a
        # stable digest. Keeps the SCNA-drift sub-RNG deterministic across
        # separate `python scripts/generate_mocks.py` invocations.
        seed_bytes = hashlib.blake2s(
            f"{spec.tumour_id_canonical}|{clone.name}|{chrom}".encode(),
            digest_size=4,
        ).digest()
        rng = _random.Random(int.from_bytes(seed_bytes, "big"))
        # 40% chance of some SCNA event on this (clone, chrom).
        if rng.random() < 0.40:
            event = rng.choices(
                ["gain_A", "gain_B", "loss_A", "loss_B", "loh_A", "loh_B"],
                weights=[0.25, 0.25, 0.20, 0.20, 0.05, 0.05],
            )[0]
            if event == "gain_A":
                a += 1
            elif event == "gain_B":
                b += 1
            elif event == "loss_A":
                a = max(0, a - 1)
            elif event == "loss_B":
                b = max(0, b - 1)
            elif event == "loh_A":
                a, b = a + b, 0
            elif event == "loh_B":
                a, b = 0, a + b
    return int(a), int(b)


def _write_alpaca(specs: list[TumourSpec]) -> None:
    ALPACA_DIR.mkdir(parents=True, exist_ok=True)
    out = ALPACA_DIR / "all_tumours_combined.csv"

    lines = [",".join(_ALPACA_SEGMENT_HEADERS)]
    for spec in specs:
        if spec.status != "resolved":
            # ALPACA is not run on unresolved tumours in real life; match that.
            continue
        segments = _tumour_segments(spec)
        for clone in spec.clones:
            for chrom, start, end in segments:
                a, b = _pick_alpaca_cn(spec, clone, chrom)
                lines.append(
                    ",".join(
                        [
                            spec.tumour_id_canonical,
                            f"{chrom}_{start}_{end}",
                            clone.name,
                            str(a),
                            str(b),
                        ]
                    )
                )
    out.write_text("\n".join(lines) + "\n")
    print(f"wrote {out}  ({len(lines) - 1} rows)")


# --- Kallisto (counts + TPM, wide parquet with a `gene_id` column) ---------
#
# One tumour region-sample column per (patient × tumour × region), plus one
# `_N01` normal for a subset of patients. Genes = the mock GENES panel'
# ensembl IDs, so lof_gene can join by ENSG.

def _write_kallisto(specs: list[TumourSpec]) -> None:
    """Emit two wide-format parquets: kallisto_counts + kallisto_tpm."""
    try:
        import polars as pl  # noqa: local import so mock gen still runs without polars for other outputs
    except ImportError:
        print("polars not available — skipping kallisto mock")
        return

    KALLISTO_DIR.mkdir(parents=True, exist_ok=True)

    gene_ids = [g[3] for g in GENES]  # ENSG…

    # Build sample list: for each resolved tumour, N regions; add a normal
    # sample for 60% of patients (deterministic via RNG).
    sample_columns: list[str] = []
    tumour_sample_map: dict[str, TumourSpec] = {}
    for spec in specs:
        if spec.status != "resolved":
            continue
        for r in range(1, spec.n_samples + 1):
            h = RNG.randbytes(16).hex()
            s = f"{spec.patient_id}_SU_T{spec.tumour_ordinal}-R{r}--{h}"
            sample_columns.append(s)
            tumour_sample_map[s] = spec

    # Normal samples (patient-level, one per selected patient).
    normal_samples_by_patient: dict[str, str] = {}
    patients = sorted({s.patient_id for s in specs if s.status == "resolved"})
    for pid in patients:
        if RNG.random() < 0.6:
            h = RNG.randbytes(16).hex()
            normal_samples_by_patient[pid] = f"{pid}_SU_N01--{h}"

    sample_columns.extend(sorted(normal_samples_by_patient.values()))

    # Emit realistic-ish TPM: baseline ~ log-normal centred on 10, with
    # per-tumour × per-gene knockouts for genes we know are homdel/LOH in
    # ALPACA mock (so the expression path can validate its own signal).
    tpm_rows = []
    count_rows = []
    for g_ens, g_meta in zip(gene_ids, GENES):
        sym, chrom, _start, _ensembl, _cat = g_meta
        tpm_row = {"gene_id": g_ens}
        count_row = {"gene_id": g_ens}
        for s in sample_columns:
            m = _SAMPLE_NAME_RE.match(s) if False else None
            # Baseline log-normal.
            baseline = max(0.0, RNG.gauss(10, 5))

            # Knockdown at TP53 (chr17) in tumour LTX0004-Tumour1 samples.
            if sym == "TP53" and s.startswith("LTX0004_"):
                baseline *= 0.15
            # Knockdown at KRAS (chr12) in LTX0005-Tumour1.
            if sym == "KRAS" and s.startswith("LTX0005_"):
                baseline *= 0.05

            tpm_row[s] = round(baseline, 4)
            count_row[s] = round(baseline * RNG.uniform(30, 80), 1)
        tpm_rows.append(tpm_row)
        count_rows.append(count_row)

    tpm_df = pl.DataFrame(tpm_rows)
    counts_df = pl.DataFrame(count_rows)

    tpm_out = KALLISTO_DIR / "kallisto_tpm.parquet"
    counts_out = KALLISTO_DIR / "kallisto_counts.parquet"
    tpm_df.write_parquet(tpm_out)
    counts_df.write_parquet(counts_out)
    print(f"wrote {tpm_out}  ({tpm_df.height} genes × {tpm_df.width - 1} samples)")
    print(f"wrote {counts_out}  ({counts_df.height} genes × {counts_df.width - 1} samples)")


# Placeholder (never matches); kept only so the block above doesn't hit an
# undefined name when we later add a real parser. Regex parsing is in
# tx_data/builds/kallisto.py; the mock just constructs plausible names.
_SAMPLE_NAME_RE = None


def _write_clone_proportions(specs: list[TumourSpec]) -> None:
    """Emit `ALPACA/input/<tumour_id>/cp_table.csv` for every resolved tumour.

    Layout matches the real source: N sample columns + trailing `clone`
    column. Values per sample column sum to 1. Ensures at least one clone
    per tumour is fully clonal (present) and at least one is absent
    (all-zero row) so downstream presence-threshold logic gets exercised.
    """
    for spec in specs:
        if spec.status != "resolved":
            continue

        tumour_dir = CLONE_PROP_ROOT / spec.tumour_id_canonical
        tumour_dir.mkdir(parents=True, exist_ok=True)

        # Sample IDs mirror the kallisto naming convention (`SU_T{n}-R{r}--<hash>`)
        # so the two mocks share a plausible cross-table shape.
        sample_ids: list[str] = []
        for r in range(1, spec.n_samples + 1):
            h = RNG.randbytes(6).hex()
            sample_ids.append(
                f"{spec.patient_id}_SU_T{spec.tumour_ordinal}-R{r}--{h}"
            )

        clone_names = [c.name for c in spec.clones]
        # proportions[sample_idx][clone_idx] → float
        proportions: dict[str, dict[str, float]] = {s: {c: 0.0 for c in clone_names} for s in sample_ids}

        # For each sample, pick a mixture of 1–3 clones (weighted to sum to 1).
        # First sample forces a fully-clonal state on the trunk (clone1) so
        # every mock tumour has at least one always-present clone.
        for idx, s in enumerate(sample_ids):
            if idx == 0:
                proportions[s][clone_names[0]] = 1.0
                continue
            k = RNG.choice([1, 2, 3])
            picks = RNG.sample(clone_names, k=min(k, len(clone_names)))
            weights = [RNG.uniform(0.1, 1.0) for _ in picks]
            tot = sum(weights)
            for c, w in zip(picks, weights):
                proportions[s][c] = round(w / tot, 4)
            # Adjust final entry so column sums to exactly 1.
            drift = 1.0 - sum(proportions[s].values())
            proportions[s][picks[-1]] = round(proportions[s][picks[-1]] + drift, 4)

        # Compose CSV: sample columns then `clone`.
        header = ",".join(sample_ids + ["clone"])
        rows = [header]
        for c in clone_names:
            cells = [f"{proportions[s][c]:.4g}" for s in sample_ids]
            rows.append(",".join(cells + [c]))

        out = tumour_dir / "cp_table.csv"
        out.write_text("\n".join(rows) + "\n")
    print(f"wrote {sum(1 for s in specs if s.status == 'resolved')} cp_table.csv files under {CLONE_PROP_ROOT}")


def _write_purity_ploidy(specs: list[TumourSpec]) -> None:
    """Emit `nemo_mock/tx842_purity_ploidy.csv`.

    One row per sequencing region for every resolved tumour. Uses the same
    `_SU_T{n}-R{r}--<hash>` convention as muttable / kallisto / clone_proportions.
    Hashes are drawn from the shared RNG, so they will NOT line up with the
    equivalent muttable rows — that's a mock limitation, not a convention
    violation. Downstream joins at mock level use (tumour_id, region) via the
    build script's derived columns, which do line up.
    """
    out = MOCK_ROOT / "tx842_purity_ploidy.csv"
    lines = ['"sample_name_hash","purity","ploidy"']
    for spec in specs:
        if spec.status != "resolved":
            continue
        for r in range(1, spec.n_samples + 1):
            h = RNG.randbytes(6).hex()
            sample_hash = f"{spec.patient_id}_SU_T{spec.tumour_ordinal}-R{r}--{h}"
            # Plausible mock values: purity ~ U(0.2, 0.95); ploidy ~ 2.0 for
            # no_wgd, ~4.0 for clonal_wgd, otherwise ~3.5. Add a touch of noise.
            if spec.class_ == "clonal_wgd":
                ploidy = round(RNG.uniform(3.7, 4.3), 2)
            elif spec.class_ == "no_wgd":
                ploidy = round(RNG.uniform(1.9, 2.3), 2)
            else:
                ploidy = round(RNG.uniform(2.8, 3.8), 2)
            purity = round(RNG.uniform(0.2, 0.95), 2)
            lines.append(f'"{sample_hash}",{purity},{ploidy}')
    out.write_text("\n".join(lines) + "\n")
    print(f"wrote {out}")


def _write_cohort_manifest(specs: list[TumourSpec]) -> None:
    """Emit `nemo_mock/_cohort_manifest.json` for cross-language mock gens.

    The R mock generator for `conipher_trees` (and any future R-side mock)
    reads this file instead of duplicating the Python-side `cohort()` logic.
    Only the fields R-side generators actually need are emitted — enough to
    re-derive clone trees, sample names, and tumour IDs in whatever
    convention the downstream mock requires.
    """
    import json

    manifest = {
        "schema_version": 1,
        "cohort_tag": "Tx842_mock",
        "random_seed": 20260831,
        "tumours": [
            {
                "patient_id": spec.patient_id,
                "tumour_ordinal": spec.tumour_ordinal,
                "tumour_id_canonical": spec.tumour_id_canonical,  # LTX0001-Tumour1
                "tumour_id_muttable": spec.tumour_id_muttable,    # LTX0001_tumour1
                "status": spec.status,
                "histology_group": spec.histology_group,
                "n_samples": spec.n_samples,
                "clones": [
                    {"name": c.name, "parent": c.parent, "gds_at_clone": c.gds_at_clone}
                    for c in spec.clones
                ],
            }
            for spec in specs
        ],
    }
    out = MOCK_ROOT / "_cohort_manifest.json"
    out.write_text(json.dumps(manifest, indent=2))
    print(f"wrote {out} ({len(manifest['tumours'])} tumours)")


def main() -> None:
    MOCK_ROOT.mkdir(parents=True, exist_ok=True)
    specs = cohort()
    print(f"cohort: {len(specs)} tumours across {len({s.patient_id for s in specs})} patients")

    _write_muttable(specs)
    _write_wgd_calls(specs)
    _write_clinical(specs)
    _write_alphamissense(specs)
    _write_driver_list()
    _write_alpaca(specs)
    _write_kallisto(specs)
    _write_clone_proportions(specs)
    _write_purity_ploidy(specs)
    _write_cohort_manifest(specs)
    print("done.")
    print(
        "note: `conipher_trees` mock (`nemo_mock/tx842_trees.RDS`) is produced by\n"
        "      `Rscript scripts/generate_mock_conipher_trees.R` — requires an R env\n"
        "      with the `jsonlite` package (e.g. the tx800_alpaca_replica pixi env)."
    )


if __name__ == "__main__":
    main()
