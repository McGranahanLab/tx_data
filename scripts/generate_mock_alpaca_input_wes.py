"""Write the WES ALPACA input mocks: `ALPACA/input/<tumour_id>/ALPACA_input_table.csv` (msai_segments, 2026-11-24).

Per tumour with a cp_table.csv in the mock, each sample's fractional allele CN on each ALPACA segment is the
clone mixture of the mock ALPACA output (`all_tumours_combined.csv` pred_CN_A/B weighted by the cp_table
proportions), plus a little deterministic noise, so the input is consistent with the output and the clone
proportions. Columns as ALPACA's input table: tumour_id, segment, sample, cpnA, cpnB. Deterministic.

  python scripts/generate_mock_alpaca_input_wes.py
"""

from __future__ import annotations

import csv
import random
from collections import defaultdict

from tx_data.paths import DEFAULT_MOCK_ROOT

ALPACA = DEFAULT_MOCK_ROOT / "project/proj-tracerx-lung/tctProjects/lungTx/Tx842/release/tracerx_842/ALPACA"


def main() -> None:
    rng = random.Random(20261124)
    cn: dict[str, dict[str, dict[str, tuple[int, int]]]] = defaultdict(lambda: defaultdict(dict))
    with open(ALPACA / "output" / "cohort_results" / "all_tumours_combined.csv") as fh:
        for r in csv.DictReader(fh):
            cn[r["tumour_id"]][r["segment"]][r["clone"]] = (int(r["pred_CN_A"]), int(r["pred_CN_B"]))
    n = 0
    for tdir in sorted((ALPACA / "input").iterdir()):
        cp_path = tdir / "cp_table.csv"
        if not cp_path.is_file() or tdir.name not in cn:
            continue
        with open(cp_path) as fh:
            rows = list(csv.reader(fh))
        samples, props = rows[0][:-1], {r[-1]: [float(x) for x in r[:-1]] for r in rows[1:]}
        out = ["tumour_id,segment,sample,cpnA,cpnB"]
        for seg in sorted(cn[tdir.name], key=lambda s: (s.split("_")[0].zfill(2), int(s.split("_")[1]))):
            for i, s in enumerate(samples):
                a = sum(props[c][i] * cn[tdir.name][seg].get(c, (1, 1))[0] for c in props)
                b = sum(props[c][i] * cn[tdir.name][seg].get(c, (1, 1))[1] for c in props)
                a, b = (round(max(v + rng.uniform(-0.05, 0.05), 0.0), 4) for v in (a, b))
                out.append(f"{tdir.name},{seg},{s},{a},{b}")
        (tdir / "ALPACA_input_table.csv").write_text("\n".join(out) + "\n")
        n += 1
    print(f"wrote ALPACA_input_table.csv for {n} WES mock tumours under {ALPACA / 'input'}")


if __name__ == "__main__":
    main()
