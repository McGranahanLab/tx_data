"""Generate the synthetic WGS mock cohort into ``nemo_mock/``.

Separate from ``generate_mocks.py`` (the WES mock). ``nemo_mock/`` mirrors
``/nemo``, so every WGS file lands at its real Nemo path under
``nemo_mock/``. The cohort spec below is the single source of truth for WGS
mocks; later generators (Battenberg subclones, MSAI, mosdepth, ...) read it
from ``nemo_mock/_cohort_manifest_wgs.json`` so all WGS mocks agree.

Cohort shape (11 tumours / 12 samples / 10 patients), chosen to exercise:
  - single-sample, single-clone tumours (the bulk of the real WGS cohort)
  - patients shared with the WES mock (LTX0001-LTX0004) with *different*
    sample hashes, plus WGS-only patients (LTX0101+)
  - one multi-region WGS tumour and one multi-tumour patient
  - one sample whose name breaks the usual `_SU_T<n>-R<n>` tokens (real
    WGS names have no closed catalogue)
  - the source's tumour-label convention in `PATIENT` and the Battenberg
    folder names: bare `LTX0001` for single-tumour patients,
    `LTX0104_tumour1` / `LTX0104_tumour2` for multi-tumour patients
  - purity from very low (0.09) to high (0.88); near-diploid and
    WGD-like ploidy (~3.6-3.9)

Copy-number model (consensus segments + ALPACA cp_table agree):
  - genome-wide hg38 segments (chr1-22, X), shared by a tumour's samples
  - clone1 (trunk) carries integer allele-specific states around a
    diploid or WGD baseline (from the tumour's ploidy); each further clone
    descends from clone1 and changes ~20% of segments
  - each sample's fractional CN = sum over clones of cp_table proportion ×
    clone state, plus Gaussian noise with sd 0.03 / purity per allele.
    So single-clone tumours sit near integers (worse at low purity) and
    multi-clone tumours have genuinely fractional segments.
  - LTX0102's segments file uses the legacy COPY_NUMBER_A / _B header.

Files written (mirrored /nemo paths under ``nemo_mock/``):
  - ``.../CN-CCF/wgs_resources/inventories/battenberg_output_inventory.tsv``
  - ``.../CN-CCF/battenberg_wgs/output/<label>/tree_input/<label>_copynumber_MSAI_adjusted_consensus_segments_tree_input.tsv``
    for multi-sample tumours; single-sample tumours can't be MSAI-adjusted and
    get ``<label>_copynumber_MSAI_NOT_adjusted_segments_tree_input.tsv`` instead
  - ``.../CN-CCF/wgs_resources/alpaca/input/<label>/cp_table.csv``
  - ``.../CN-CCF/wgs_resources/alpaca/input/<label>/tree_paths.json`` (ALPACA
    list-of-paths, trunk → leaf; every subclone is a child of clone1)
  - ``_cohort_manifest_wgs.json`` (mock-only metadata, not part of the mirror)

Deterministic (fixed seed) so mocks re-generate identically.
"""

from __future__ import annotations

import json
import random
from dataclasses import dataclass

from tx_data.paths import REPO_ROOT

MOCK_DIR = REPO_ROOT / "nemo_mock"
MOCK_CN_CCF = MOCK_DIR / "project" / "proj-tracerx-lung" / "tctProjects" / "CN-CCF"
INVENTORY = MOCK_CN_CCF / "wgs_resources" / "inventories" / "battenberg_output_inventory.tsv"
MOCK_BATTENBERG_OUT = MOCK_CN_CCF / "battenberg_wgs" / "output"
MOCK_ALPACA_INPUT = MOCK_CN_CCF / "wgs_resources" / "alpaca" / "input"

# Real Nemo locations the inventory points at (kept as /nemo strings, as on source).
BATTENBERG_OUT = "/nemo/project/proj-tracerx-lung/tctProjects/CN-CCF/battenberg_wgs/output"
MOSDEPTH_DIR = "/nemo/project/proj-tracerx-wgs/working/liul/mitochondria/data/mosdepth"

SEED = 20261027
RNG = random.Random(SEED)
# Separate stream for copy number so CN changes never reshuffle sample hashes.
CN_RNG = random.Random(SEED + 1)

# hg38 chromosome lengths (bp).
HG38_LENGTHS = {
    "1": 248956422, "2": 242193529, "3": 198295559, "4": 190214555,
    "5": 181538259, "6": 170805979, "7": 159345973, "8": 145138636,
    "9": 138394717, "10": 133797422, "11": 135086622, "12": 133275309,
    "13": 114364328, "14": 107043718, "15": 101991189, "16": 90338345,
    "17": 83257441, "18": 80373285, "19": 58617616, "20": 64444167,
    "21": 46709983, "22": 50818468, "X": 156040895,
}
NOISE_SD_AT_PURITY_1 = 0.03
LEGACY_HEADER_LABELS = {"LTX0102"}


@dataclass
class WgsSample:
    region: str
    purity: float
    # Replaces the usual `T<n>-R<n>` token to mimic non-standard names.
    name_token: str | None = None


@dataclass
class WgsTumour:
    patient_id: str
    tumour_ordinal: int
    ploidy: float
    n_clones: int
    samples: list[WgsSample]

    # Set by cohort(): True when the patient has more than one WGS tumour.
    multi_tumour_patient: bool = False

    @property
    def tumour_id_canonical(self) -> str:
        return f"{self.patient_id}-Tumour{self.tumour_ordinal}"

    @property
    def wgs_tumour_label(self) -> str:
        """Source `PATIENT` value / Battenberg folder name for this tumour."""
        if self.multi_tumour_patient:
            return f"{self.patient_id}_tumour{self.tumour_ordinal}"
        return self.patient_id


def cohort() -> list[WgsTumour]:
    def one(region_purity: float) -> list[WgsSample]:
        return [WgsSample("R1", region_purity)]

    specs = [
        # Shared with the WES mock; mirrors the real head row (rho ~0.26, psi ~3.86).
        WgsTumour("LTX0001", 1, 3.86072434911302, 1, one(0.25893)),
        WgsTumour("LTX0002", 1, 1.98, 1, one(0.71)),
        WgsTumour("LTX0003", 1, 2.10, 1, one(0.14)),       # low purity
        WgsTumour("LTX0004", 1, 3.20, 3, one(0.55)),       # subclonal
        # WGS-only patients.
        WgsTumour("LTX0101", 1, 2.00, 1, one(0.88)),       # high purity, diploid
        WgsTumour("LTX0102", 1, 3.60, 1, one(0.33)),       # WGD-like
        WgsTumour("LTX0103", 1, 2.90, 2,                   # multi-region WGS
                  [WgsSample("R1", 0.60), WgsSample("R2", 0.45)]),
        WgsTumour("LTX0104", 1, 2.30, 1, one(0.50)),       # multi-tumour patient
        WgsTumour("LTX0104", 2, 2.70, 2, one(0.40)),
        WgsTumour("LTX0105", 1, 2.40, 1, one(0.09)),       # very low purity
        WgsTumour("LTX0106", 1, 2.20, 1,                   # non-standard sample name
                  [WgsSample("R1", 0.47, name_token="LN1")]),
    ]
    per_patient: dict[str, int] = {}
    for t in specs:
        per_patient[t.patient_id] = per_patient.get(t.patient_id, 0) + 1
    for t in specs:
        t.multi_tumour_patient = per_patient[t.patient_id] > 1
    return specs


def _hash() -> str:
    return "".join(RNG.choice("0123456789abcdef") for _ in range(10))


def _sample_name(t: WgsTumour, s: WgsSample, h: str) -> str:
    token = s.name_token or f"T{t.tumour_ordinal}-{s.region}"
    return f"{t.patient_id}_SU_{token}--{h}"


def _write_inventory(specs: list[WgsTumour]) -> list[dict]:
    header = [
        "PATIENT", "TUMOUR_SAMPLE", "results_dir", "optimal_solution", "subclones",
        "purity", "ploidy", "msai", "mosdepth_summary", "coverage",
    ]
    rows: list[dict] = []
    for t in specs:
        for s in t.samples:
            name = _sample_name(t, s, _hash())
            # Solution label carries rounded psi/rho that differ slightly from
            # the precise purity/ploidy columns, as on the real file.
            psi = round(t.ploidy + RNG.uniform(-0.02, 0.02), 2)
            rho = round(s.purity + RNG.uniform(-0.01, 0.015), 2)
            solution = f"purpleDefault_psi{psi}_rho{rho}"
            rows.append({
                "PATIENT": t.wgs_tumour_label,
                "TUMOUR_SAMPLE": name,
                "results_dir": f"{BATTENBERG_OUT}/{t.wgs_tumour_label}/results/",
                "optimal_solution": solution,
                "subclones": f"{name}_{solution}_subclones.txt.gz",
                "purity": s.purity,
                "ploidy": t.ploidy,
                # Doubled slashes reproduce the real file's quirk.
                "msai": f"{BATTENBERG_OUT}//{t.wgs_tumour_label}/results//MSAI/multisample_MSAI.txt",
                "mosdepth_summary": f"{MOSDEPTH_DIR}/{name}.mosdepth.summary.txt",
                "coverage": round(RNG.uniform(60, 250), 2),
            })

    INVENTORY.parent.mkdir(parents=True, exist_ok=True)
    lines = ["\t".join(header)]
    lines += ["\t".join(str(r[c]) for c in header) for r in rows]
    INVENTORY.write_text("\n".join(lines) + "\n")
    print(f"wrote {INVENTORY} ({len(rows)} samples)")
    return rows


def _write_manifest(specs: list[WgsTumour], rows: list[dict], truth: dict[str, dict]) -> None:
    # Rows were emitted in spec order, one per sample.
    it = iter(rows)

    manifest = {
        "schema_version": 1,
        "cohort": "wgs",
        "cohort_tag": "TxWGS_mock",
        "genome_build": "hg38",
        "random_seed": SEED,
        "cn_noise_sd_at_purity_1": NOISE_SD_AT_PURITY_1,
        "tumours": [
            {
                "patient_id": t.patient_id,
                "tumour_ordinal": t.tumour_ordinal,
                "tumour_id_canonical": t.tumour_id_canonical,
                "wgs_tumour_label": t.wgs_tumour_label,
                "ploidy": t.ploidy,
                "n_clones": t.n_clones,
                **truth[t.wgs_tumour_label],
                "samples": [
                    {
                        "sample_name_hash": r["TUMOUR_SAMPLE"],
                        "region": s.region,
                        "standard_name": s.name_token is None,
                        "purity": r["purity"],
                        "optimal_solution": r["optimal_solution"],
                    }
                    for s, r in ((s, next(it)) for s in t.samples)
                ],
            }
            for t in specs
        ],
    }
    out = MOCK_DIR / "_cohort_manifest_wgs.json"
    out.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"wrote {out} ({len(specs)} tumours)")


def _segments_source_name(sample_name_hash: str) -> str:
    """Short SAMPLE form used in tree_input files: stem with '-' → '.'."""
    return sample_name_hash.split("--", 1)[0].replace("-", ".")


def _segments() -> list[tuple[str, int, int]]:
    segs: list[tuple[str, int, int]] = []
    for chrom, length in HG38_LENGTHS.items():
        n = CN_RNG.randint(2, 6)
        cuts = sorted(CN_RNG.sample(range(1_000_000, length - 1_000_000), n - 1))
        bounds = [10_001] + cuts + [length - 10_000]
        for a, b in zip(bounds[:-1], bounds[1:]):
            segs.append((chrom, a, b - 1))
    return segs


def _clone_states(t: WgsTumour, n_segs: int) -> dict[str, list[tuple[int, int]]]:
    """Integer (major, minor) per segment for every clone."""
    base = 2 if t.ploidy >= 3 else 1
    trunk: list[tuple[int, int]] = []
    for _ in range(n_segs):
        major, minor = base, base
        r = CN_RNG.random()
        if r < 0.20:
            major += 1                      # gain
        elif r < 0.35:
            minor = 0                       # LOH
        elif r < 0.42:
            major, minor = major + 2, minor  # amplification
        trunk.append((major, minor))
    states = {"clone1": trunk}
    for c in range(2, t.n_clones + 1):
        sub = []
        for major, minor in trunk:
            if CN_RNG.random() < 0.20:
                if minor > 0 and CN_RNG.random() < 0.5:
                    minor -= 1
                else:
                    major += 1
            sub.append((major, minor))
        states[f"clone{c}"] = sub
    return states


def _clone_proportions(t: WgsTumour, n_samples: int) -> list[dict[str, float]]:
    if t.n_clones == 1:
        return [{"clone1": 1.0} for _ in range(n_samples)]
    out = []
    for _ in range(n_samples):
        w = [CN_RNG.uniform(0.15, 1.0) for _ in range(t.n_clones)]
        tot = sum(w)
        props = [round(x / tot, 4) for x in w]
        props[-1] = round(1.0 - sum(props[:-1]), 4)
        out.append({f"clone{i + 1}": p for i, p in enumerate(props)})
    return out


def _write_cn(specs: list[WgsTumour], rows: list[dict]) -> dict[str, dict]:
    """Write consensus-segment and cp_table mocks; return per-label truth."""
    truth: dict[str, dict] = {}
    it = iter(rows)
    for t in specs:
        label = t.wgs_tumour_label
        samples = [(s, next(it)) for s in t.samples]
        segs = _segments()
        states = _clone_states(t, len(segs))
        props = _clone_proportions(t, len(samples))

        legacy = label in LEGACY_HEADER_LABELS
        cols = ("COPY_NUMBER_A", "COPY_NUMBER_B") if legacy else (
            "COPY_NUMBER_MAJOR", "COPY_NUMBER_MINOR")
        lines = ["\t".join(("SAMPLE", "CHR", "STARTPOS", "ENDPOS") + cols)]
        for (s, row), cp in zip(samples, props):
            sd = NOISE_SD_AT_PURITY_1 / s.purity
            src_name = _segments_source_name(row["TUMOUR_SAMPLE"])
            for i, (chrom, start, end) in enumerate(segs):
                major = sum(p * states[c][i][0] for c, p in cp.items())
                minor = sum(p * states[c][i][1] for c, p in cp.items())
                major = max(0.0, major + CN_RNG.gauss(0, sd))
                minor = max(0.0, minor + CN_RNG.gauss(0, sd))
                # Major is the larger allele by definition, as on real files.
                major, minor = max(major, minor), min(major, minor)
                lines.append(f"{src_name}\t{chrom}\t{start}\t{end}\t{major!r}\t{minor!r}")
        adjusted = len(samples) > 1
        out = MOCK_BATTENBERG_OUT / label / "tree_input" / (
            f"{label}_copynumber_MSAI_adjusted_consensus_segments_tree_input.tsv" if adjusted
            else f"{label}_copynumber_MSAI_NOT_adjusted_segments_tree_input.tsv")
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text("\n".join(lines) + "\n")

        # cp_table.csv: N sample columns (full sample names, as WES) + trailing clone.
        names = [row["TUMOUR_SAMPLE"] for _, row in samples]
        cp_lines = [",".join(names + ["clone"])]
        for c in sorted(states, key=lambda x: int(x[5:])):
            cp_lines.append(",".join([str(cp.get(c, 0.0)) for cp in props] + [c]))
        cp_out = MOCK_ALPACA_INPUT / label / "cp_table.csv"
        cp_out.parent.mkdir(parents=True, exist_ok=True)
        cp_out.write_text("\n".join(cp_lines) + "\n")

        # tree_paths.json: ALPACA list of trunk → leaf paths.
        subs = [c for c in states if c != "clone1"]
        paths = [["clone1", c] for c in subs] or [["clone1"]]
        (MOCK_ALPACA_INPUT / label / "tree_paths.json").write_text(json.dumps(paths) + "\n")

        truth[label] = {
            "n_segments": len(segs),
            "segments_header": "legacy_A_B" if legacy else "major_minor",
            "msai_adjusted": adjusted,
            "segments_sample_names": [_segments_source_name(n) for n in names],
            "clone_proportions": dict(zip(names, props)),
            "tree_paths": paths,
        }
    print(f"wrote MSAI-adjusted consensus segments + cp_table + tree_paths for {len(truth)} tumours")
    return truth


def main() -> None:
    specs = cohort()
    rows = _write_inventory(specs)
    truth = _write_cn(specs, rows)
    _write_manifest(specs, rows, truth)


if __name__ == "__main__":
    main()
