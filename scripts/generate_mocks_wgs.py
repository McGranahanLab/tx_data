"""Generate the synthetic WGS mock cohort into ``nemo_mock/``.

Separate from ``generate_mocks.py`` (the WES mock). ``nemo_mock/`` mirrors
``/nemo``, so every WGS file lands at its real Nemo path under
``nemo_mock/``. The cohort spec below is the single source of truth for WGS
mocks; later generators (Battenberg subclones, MSAI, mosdepth, ...) read it
from ``nemo_mock/_cohort_manifest_wgs.json`` so all WGS mocks agree.

Cohort shape (10 tumours / 11 samples / 9 patients), chosen to exercise:
  - single-sample, single-clone tumours (the bulk of the real WGS cohort)
  - patients shared with the WES mock (LTX0001-LTX0004) with *different*
    sample hashes, plus WGS-only patients (LTX0101+)
  - one multi-region WGS tumour and one multi-tumour patient
  - purity from very low (0.09) to high (0.88); near-diploid and
    WGD-like ploidy (~3.6-3.9)

Files written:
  - ``project/proj-tracerx-wgs/working/liul/TxWGS/inventory/2026-09-09_battenberg_output_inventory.tsv``
  - ``_cohort_manifest_wgs.json`` (mock-only metadata, not part of the mirror)

Deterministic (fixed seed) so mocks re-generate identically.
"""

from __future__ import annotations

import json
import random
from dataclasses import dataclass

from tx_data.paths import REPO_ROOT

MOCK_DIR = REPO_ROOT / "nemo_mock"
INVENTORY = (
    MOCK_DIR / "project" / "proj-tracerx-wgs" / "working" / "liul" / "TxWGS"
    / "inventory" / "2026-09-09_battenberg_output_inventory.tsv"
)

# Real Nemo locations the inventory points at (kept as /nemo strings, as on source).
BATTENBERG_OUT = "/nemo/project/proj-tracerx-lung/tctProjects/CN-CCF/battenberg_wgs/output"
MOSDEPTH_DIR = "/nemo/project/proj-tracerx-wgs/working/liul/mitochondria/data/mosdepth"

SEED = 20261027
RNG = random.Random(SEED)


@dataclass
class WgsSample:
    region: str
    purity: float


@dataclass
class WgsTumour:
    patient_id: str
    tumour_ordinal: int
    ploidy: float
    n_clones: int
    samples: list[WgsSample]

    @property
    def tumour_id_canonical(self) -> str:
        return f"{self.patient_id}-Tumour{self.tumour_ordinal}"


def cohort() -> list[WgsTumour]:
    def one(region_purity: float) -> list[WgsSample]:
        return [WgsSample("R1", region_purity)]

    return [
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
    ]


def _hash() -> str:
    return "".join(RNG.choice("0123456789abcdef") for _ in range(10))


def _sample_name(t: WgsTumour, s: WgsSample, h: str) -> str:
    return f"{t.patient_id}_SU_T{t.tumour_ordinal}-{s.region}--{h}"


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
                "PATIENT": t.patient_id,
                "TUMOUR_SAMPLE": name,
                "results_dir": f"{BATTENBERG_OUT}/{t.patient_id}/results/",
                "optimal_solution": solution,
                "subclones": f"{name}_{solution}_subclones.txt.gz",
                "purity": s.purity,
                "ploidy": t.ploidy,
                # Doubled slashes reproduce the real file's quirk.
                "msai": f"{BATTENBERG_OUT}//{t.patient_id}/results//MSAI/multisample_MSAI.txt",
                "mosdepth_summary": f"{MOSDEPTH_DIR}/{name}.mosdepth.summary.txt",
                "coverage": round(RNG.uniform(60, 250), 2),
            })

    INVENTORY.parent.mkdir(parents=True, exist_ok=True)
    lines = ["\t".join(header)]
    lines += ["\t".join(str(r[c]) for c in header) for r in rows]
    INVENTORY.write_text("\n".join(lines) + "\n")
    print(f"wrote {INVENTORY} ({len(rows)} samples)")
    return rows


def _write_manifest(specs: list[WgsTumour], rows: list[dict]) -> None:
    by_tumour: dict[tuple[str, int], list[dict]] = {}
    for r in rows:
        name = r["TUMOUR_SAMPLE"]
        ordinal = int(name.split("_SU_T", 1)[1].split("-", 1)[0])
        by_tumour.setdefault((r["PATIENT"], ordinal), []).append(r)

    manifest = {
        "schema_version": 1,
        "cohort": "wgs",
        "cohort_tag": "TxWGS_mock",
        "genome_build": "hg38",
        "random_seed": SEED,
        "tumours": [
            {
                "patient_id": t.patient_id,
                "tumour_ordinal": t.tumour_ordinal,
                "tumour_id_canonical": t.tumour_id_canonical,
                "ploidy": t.ploidy,
                "n_clones": t.n_clones,
                "samples": [
                    {
                        "sample_name_hash": r["TUMOUR_SAMPLE"],
                        "region": r["TUMOUR_SAMPLE"].split("-", 1)[1].split("--", 1)[0],
                        "purity": r["purity"],
                        "optimal_solution": r["optimal_solution"],
                    }
                    for r in by_tumour[(t.patient_id, t.tumour_ordinal)]
                ],
            }
            for t in specs
        ],
    }
    out = MOCK_DIR / "_cohort_manifest_wgs.json"
    out.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"wrote {out} ({len(specs)} tumours)")


def main() -> None:
    specs = cohort()
    rows = _write_inventory(specs)
    _write_manifest(specs, rows)


if __name__ == "__main__":
    main()
