# Synthetic CONIPHER cases

Fully synthetic tumours simulated from a known clone tree, for testing
[CONIPHER](https://github.com/McGranahanLab/CONIPHER). Not derived from any patient data.
Unlike `nemo_mock/`, this does not mirror a `/Nemo` file; it is a benchmark set.

Regenerate (byte-identical, fixed seed):

```bash
pixi run python scripts/generate_conipher_synthetic.py
```

| Case | Tree | Regions | Mutations | Notes |
|---|---|---|---|---|
| SYN001 | linear 1→2→3 | 2 | 70 | diploid |
| SYN002 | 1→{2,3}, 2→4 | 3 | 85 | clone4 private to R1; clone2 absent from R3 |
| SYN003 | 1→{2,3} | 2 | 80 | clonal gain on chr5/chr8 (cn 2+1); purity 0.35 in R1 |

Per case:

- `input.tsv`: CONIPHER input table, same columns as CONIPHER's `inst/extdata/input_table.tsv`.
  `CLUSTER` is the true clone index; `CCF_OBS`, `CCF_PHYLO`, `MUT_COPY` are estimated from the
  simulated counts. CONIPHER's clustering ignores these four columns, so one file serves both
  end-to-end and tree-building-only runs.
- `truth_tree.tsv`: true clone → parent edges, mutation counts and per-region CCF.
- `truth_mutations.tsv`: true clone and multiplicity per mutation.

Model: `VAF = purity · CCF · multiplicity / (purity · (cn_a + cn_b) + 2 · (1 − purity))`,
`var_count ~ Binomial(depth, VAF)`, `depth ~ N(150, 37.5)` (min 30). Every chromosome 1–22
carries at least one mutation: CONIPHER's plotting (`fun.chrom.length`) fails on gaps.

These files are copied into CONIPHER's `tests/testthat/` fixtures; CONIPHER does not read this repo.
