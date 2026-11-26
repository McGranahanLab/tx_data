"""Common helpers for build scripts."""

from __future__ import annotations

from pathlib import Path

from tx_data.paths import REPO_ROOT

DATA_DIR = REPO_ROOT / "data"


def canonical_output_path(name: str) -> Path:
    """Where a canonical Parquet for a table lives."""
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    return DATA_DIR / f"{name}.parquet"


def log(msg: str) -> None:
    """Tiny logger — one function so build output is easy to grep."""
    print(f"[tx_data.build] {msg}", flush=True)


# ---- per-tumour build issues (Piotr, 2026-11-26) ------------------------------------------------
# Some tumours have missing or unusable inputs. Per-tumour builders don't stop on them: they leave the
# tumour out (severity `skipped`), keep it with something missing (`partial`) or note an oddity (`warning`),
# and record why. Each build of a table replaces data/_issues/<table>.tsv; build_all combines them into
# data/build_issues.tsv and the build report.
ISSUES_DIR = DATA_DIR / "_issues"
ISSUE_COLS = ("table", "tumour_id", "severity", "reason", "path")
SEVERITIES = ("skipped", "partial", "warning")


class Issues:
    """Collects the per-tumour issues of one table's build; `write()` at the end of the build."""

    def __init__(self, table: str) -> None:
        self.table, self.rows = table, []

    def add(self, tumour_id: str, severity: str, reason: str, path: object = "") -> None:
        assert severity in SEVERITIES, severity
        self.rows.append({"table": self.table, "tumour_id": tumour_id, "severity": severity,
                          "reason": str(reason).replace("\t", " ").replace("\n", " "), "path": str(path)})

    def skipped(self) -> set[str]:
        return {r["tumour_id"] for r in self.rows if r["severity"] == "skipped"}

    def write(self) -> Path:
        import polars as pl

        ISSUES_DIR.mkdir(parents=True, exist_ok=True)
        p = ISSUES_DIR / f"{self.table}.tsv"
        pl.DataFrame(self.rows, schema={c: pl.String for c in ISSUE_COLS}).write_csv(p, separator="\t")
        if self.rows:
            counts = {s: sum(r["severity"] == s for r in self.rows) for s in SEVERITIES}
            log(f"{self.table}: tumour issues " + ", ".join(f"{s}={n}" for s, n in counts.items() if n)
                + f" (see {p})")
            for r in self.rows:
                log(f"  {r['severity']:<8} {r['tumour_id']}: {r['reason']}")
        return p
