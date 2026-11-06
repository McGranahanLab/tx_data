#!/usr/bin/env Rscript
# Generate synthetic ParallelGDDetect (PGDD) outputs for the WES mock cohort.
#
# Real location (old run, pre-new-ALPACA; to be regenerated with new hg19 ALPACA, then hg38):
#   /nemo/project/proj-tracerx-lung/tctProjects/lungTx/Tx842/release/tracerx_842/WGD/release_tx842/pgdd
# Mock location (mirrors /nemo):
#   nemo_mock/project/proj-tracerx-lung/tctProjects/lungTx/Tx842/release/tracerx_842/WGD/release_tx842/pgdd
#
# Unlike the other mocks, the output tables are NOT hand-written: this script
# builds a synthetic per-mutation PGDD *input* table from hand-designed
# scenarios (below), then runs the real `detect_par_gd()` with
# bootstrap_discovery = TRUE, so column names, types and the stable/unstable
# split come straight from PGDD. The input is kept at
# `synthetic/pgdd/pgdd_input_mock.tsv`; the scenario truth at
# `synthetic/pgdd/pgdd_mock_scenarios.json`.
#
# PGDD runs on clustered mutations BEFORE tree building, so it can call WGDs
# on clusters that CONIPHER later drops from the tree. The scenarios include
# such clusters (`in_tree = FALSE`) so downstream pruning can be tested.
#
# Requires: R + data.table + jsonlite, and a checkout of amf71/ParallelGDDetect
# (set PGDD_SRC, default ~/scratch/PGDD_amf71). Tested at master cb0437b.
#
# Run:
#   PGDD_SRC=/path/to/ParallelGDDetect Rscript scripts/generate_mock_pgdd.R
#
# Deterministic: fixed seed (input) + PGDD's own `seed` argument (bootstrap).

suppressPackageStartupMessages({
  library(data.table)
})

.script_dir <- (function() {
  all_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", all_args[grepl("^--file=", all_args)])
  if (length(file_arg) == 1 && nzchar(file_arg)) {
    return(normalizePath(dirname(file_arg)))
  }
  getwd()
})()
repo_root <- normalizePath(file.path(.script_dir, ".."))

pgdd_src <- Sys.getenv("PGDD_SRC", file.path(Sys.getenv("HOME"), "scratch", "PGDD_amf71"))
source(file.path(pgdd_src, "R", "parrallel_gd_caller.R"))

wes_root <- file.path(repo_root, "nemo_mock", "project", "proj-tracerx-lung", "tctProjects",
                      "lungTx", "Tx842", "release", "tracerx_842")
out_dir <- file.path(wes_root, "WGD", "release_tx842", "pgdd")
syn_dir <- file.path(repo_root, "synthetic", "pgdd")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(syn_dir, recursive = TRUE, showWarnings = FALSE)

SEED <- 20261106L
set.seed(SEED)

# ---- sample ids + tree clones: from the cohort manifest (generate_mocks.py) ----
# Same per-region hash as the muttable / purity_ploidy / cp_table / CONIPHER RDS,
# in region order (R1, R2, ...), so `gds` below lines up with the region.
manifest <- jsonlite::fromJSON(file.path(repo_root, "nemo_mock", "_cohort_manifest.json"),
                               simplifyVector = FALSE)
if (!identical(as.integer(manifest$schema_version), 2L)) {
  stop("cohort manifest schema_version must be 2; re-run `python scripts/generate_mocks.py`.")
}
manifest_tumours <- setNames(manifest$tumours,
                             vapply(manifest$tumours, function(t) t$tumour_id_muttable, ""))
samples_by_tumour <- lapply(manifest_tumours, function(t)
  vapply(t$samples, function(s) s$sample_name_hash, ""))
# Cross-check against the muttable the manifest was written alongside.
mt <- unique(fread(cmd = paste("zcat", shQuote(file.path(wes_root, "tx842_mutation_table.tsv.gz"))),
                   select = c("patient_tumour", "sample_name_hash")))
for (t in names(samples_by_tumour)) {
  stopifnot(setequal(samples_by_tumour[[t]], mt[patient_tumour == t, sample_name_hash]))
}

# ---- scenarios ----
# gds: num_gds per sample, in region order (R1, R2, ...). Consistent with mock
#      purity_ploidy (ploidy >= 3 -> 1 GD) where that table has the tumour.
# clusters: id, clonal?, in final CONIPHER tree?, number of mutations,
#           per-sample presence (TRUE/FALSE) and per-sample doubled fraction
#           (fraction of present mutations at mut_cpn ~2).
cl <- function(id, n, present, dbl, clonal = FALSE, in_tree = TRUE) {
  list(id = as.character(id), n = n, present = present, dbl = dbl,
       clonal = clonal, in_tree = in_tree)
}
LOW <- 0.04  # background rate of mut_cpn ~2 (noise, post-WGD gains)

scenarios <- list(
  LTX0001_tumour1 = list(
    note = "no WGD",
    gds = c(0, 0, 0),
    clusters = list(cl(1, 60, c(T, T, T), c(0, 0, 0), clonal = TRUE),
                    cl(2, 25, c(T, T, F), c(0, 0, 0)),
                    cl(3, 25, c(T, F, F), c(0, 0, 0)))),
  LTX0002_tumour1 = list(
    note = "clonal WGD",
    gds = c(1, 1, 1),
    clusters = list(cl(1, 60, c(T, T, T), c(.6, .6, .6), clonal = TRUE),
                    cl(2, 25, c(T, T, T), rep(LOW, 3)),
                    cl(3, 25, c(T, T, F), rep(LOW, 3)),
                    cl(4, 25, c(T, F, F), rep(LOW, 3)))),
  LTX0003_tumour1 = list(
    note = "clonal WGD",
    gds = c(1, 1),
    clusters = list(cl(1, 60, c(T, T), c(.55, .55), clonal = TRUE),
                    cl(2, 25, c(T, T), rep(LOW, 2)),
                    cl(3, 25, c(F, T), rep(LOW, 2)))),
  LTX0004_tumour1 = list(
    note = "single subclonal WGD (R1,R2,R4) supported by tree cluster 3; borderline cluster 5 in R4 (expected unstable)",
    gds = c(1, 1, 0, 1),
    clusters = list(cl(1, 60, c(T, T, T, T), rep(LOW, 4), clonal = TRUE),
                    cl(2, 25, c(T, T, T, T), rep(LOW, 4)),
                    cl(3, 25, c(T, T, F, T), c(.5, .5, 0, .5)),
                    cl(4, 25, c(F, F, T, F), rep(LOW, 4)),
                    cl(5, 12, c(F, F, F, T), c(0, 0, 0, .5)))),
  LTX0005_tumour1 = list(
    note = "parallel subclonal WGDs: R1 supported by clusters 3 + 6 ('3,6'; 6 not in tree), R3+R4 by cluster 4",
    gds = c(1, 0, 1, 1),
    clusters = list(cl(1, 60, c(T, T, T, T), rep(LOW, 4), clonal = TRUE),
                    cl(2, 25, c(T, T, T, T), rep(LOW, 4)),
                    cl(3, 25, c(T, F, F, F), c(.5, 0, 0, 0)),
                    cl(4, 25, c(F, F, T, T), c(0, 0, .5, .5)),
                    cl(5, 25, c(F, T, F, F), rep(LOW, 4)),
                    cl(6, 20, c(T, F, F, F), c(.5, 0, 0, 0), in_tree = FALSE))),
  LTX0006_tumour1 = list(
    note = "subclonal WGD (R1,R2) supported ONLY by cluster 7, which is not in the tree -> ploidy-only after pruning",
    gds = c(1, 1, 0),
    clusters = list(cl(1, 60, c(T, T, T), rep(LOW, 3), clonal = TRUE),
                    cl(2, 25, c(T, T, F), rep(LOW, 3)),
                    cl(3, 25, c(F, F, T), rep(LOW, 3)),
                    cl(4, 25, c(T, F, F), rep(LOW, 3)),
                    cl(7, 22, c(T, T, F), c(.5, .5, 0), in_tree = FALSE))),
  LTX0007_tumour1 = list(
    note = "parallel subclonal WGDs: R2 by cluster 2, R3 by cluster 12 (not in tree). Cluster ids 2/12 exercise PGDD's grepl substring match (R/parrallel_gd_caller.R:412)",
    gds = c(0, 1, 1),
    clusters = list(cl(1, 60, c(T, T, T), rep(LOW, 3), clonal = TRUE),
                    cl(2, 30, c(F, T, F), c(0, .65, 0)),
                    cl(3, 25, c(T, F, F), rep(LOW, 3)),
                    cl(4, 25, c(F, T, T), rep(LOW, 3)),
                    cl(5, 25, c(F, F, T), rep(LOW, 3)),
                    cl(12, 22, c(F, F, T), c(0, 0, .5), in_tree = FALSE))),
  LTX0008_tumour1 = list(
    note = "clonal WGD",
    gds = c(1, 1),
    clusters = list(cl(1, 60, c(T, T), c(.6, .6), clonal = TRUE),
                    cl(2, 25, c(T, F), rep(LOW, 2)),
                    cl(3, 25, c(F, T), rep(LOW, 2)))),
  LTX0008_tumour2 = list(
    note = "no WGD",
    gds = c(0, 0),
    clusters = list(cl(1, 60, c(T, T), c(0, 0), clonal = TRUE),
                    cl(2, 25, c(T, F), c(0, 0)))),
  LTX0009_tumour1 = list(
    note = "subclonal WGD (R1, cluster 2); tumour has no CONIPHER tree (unresolved) -> drop",
    gds = c(1, 0),
    clusters = list(cl(1, 60, c(T, T), rep(LOW, 2), clonal = TRUE, in_tree = FALSE),
                    cl(2, 30, c(T, F), c(.65, 0), in_tree = FALSE))),
  LTX0010_tumour1 = list(
    note = "clonal first WGD + subclonal second WGD in R1 (cluster 3); tumour has no CONIPHER tree -> drop",
    gds = c(2, 1),
    clusters = list(cl(1, 60, c(T, T), c(.5, .5), clonal = TRUE, in_tree = FALSE),
                    cl(2, 25, c(T, T), rep(LOW, 2), in_tree = FALSE),
                    cl(3, 25, c(T, F), c(.5, 0), in_tree = FALSE)))
)

stopifnot(setequal(names(scenarios), names(samples_by_tumour)))

# Scenario clusters must match the tree: in_tree clusters are exactly the
# manifest tree clones (= ALPACA / CONIPHER / clone_info), in_tree = FALSE are
# the manifest non_tree_clusters (or everything, for tumours with no tree).
for (tumour in names(scenarios)) {
  mt_t <- manifest_tumours[[tumour]]
  ids <- vapply(scenarios[[tumour]]$clusters, function(c) c$id, "")
  in_tree <- vapply(scenarios[[tumour]]$clusters, function(c) c$in_tree, TRUE)
  if (isTRUE(mt_t$has_tree)) {
    tree <- sub("^clone", "", vapply(mt_t$clones, function(c) c$name, ""))
    off <- as.character(unlist(mt_t$non_tree_clusters))
  } else {
    tree <- character(0)
    off <- ids
  }
  if (!setequal(ids[in_tree], tree) || !setequal(ids[!in_tree], off)) {
    stop(sprintf("%s: scenario clusters (tree: %s; off: %s) disagree with manifest (tree: %s; off: %s)",
                 tumour, toString(ids[in_tree]), toString(ids[!in_tree]), toString(tree), toString(off)))
  }
}

# Scenario keys follow the muttable form (LTX0001_tumour1, used to look up samples);
# PGDD is run on the canonical form used by the lab's PGDD runs (LTX0001-Tumour1).
canonical_tumour_id <- function(x) sub("_tumour", "-Tumour", x, fixed = TRUE)

# ---- build per-mutation input ----
bases <- c("A", "C", "G", "T")
rows <- list()
for (tumour in names(scenarios)) {
  sc <- scenarios[[tumour]]
  smp <- samples_by_tumour[[tumour]]
  stopifnot(length(smp) == length(sc$gds))
  for (c in sc$clusters) {
    chrom <- sample(as.character(1:22), c$n, replace = TRUE)
    pos <- sample.int(2e8L, c$n)
    ref <- sample(bases, c$n, replace = TRUE)
    alt <- vapply(ref, function(r) sample(setdiff(bases, r), 1), "")
    for (i in seq_along(smp)) {
      g <- sc$gds[i]
      # most mutations sit at the expected major CN for this sample's GD state
      majcn <- ifelse(runif(c$n) < 0.85, 2^g, pmax(1, 2^g + sample(c(-1, 1), c$n, replace = TRUE)))
      if (c$present[i]) {
        doubled <- runif(c$n) < c$dbl[i]
        mut_cpn <- ifelse(doubled, rnorm(c$n, 2, 0.15), rnorm(c$n, 1, 0.12))
      } else {
        mut_cpn <- abs(rnorm(c$n, 0, 0.03))
      }
      rows[[length(rows) + 1]] <- data.table(
        tumour_id = canonical_tumour_id(tumour), chromosome = chrom, position = pos, ref = ref, alt = alt,
        cluster_id = c$id, is_clonal_cluster = c$clonal, sample_id = smp[i],
        mut_cpn = round(pmax(mut_cpn, 0), 4), MajCN = majcn, num_gds = g)
    }
  }
}
input <- rbindlist(rows)
fwrite(input, file.path(syn_dir, "pgdd_input_mock.tsv"), sep = "\t")

# ---- scenario truth for tests ----
truth <- lapply(names(scenarios), function(t) {
  sc <- scenarios[[t]]
  list(tumour_id = canonical_tumour_id(t), note = sc$note,
       samples = samples_by_tumour[[t]], num_gds = sc$gds,
       clusters = lapply(sc$clusters, function(c) list(cluster_id = c$id, clonal = c$clonal,
                                                       in_tree = c$in_tree, n = c$n)))
})
jsonlite::write_json(list(seed = SEED, pgdd_commit = "cb0437b", tumours = truth),
                     file.path(syn_dir, "pgdd_mock_scenarios.json"),
                     auto_unbox = TRUE, pretty = TRUE)

# ---- run PGDD (bootstrap on, package defaults otherwise) ----
out <- detect_par_gd(copy(input), bootstrap_discovery = TRUE, n_boot = 500, seed = 1)

# File names as written by the lab's run wrapper. The unstable_calls_* files are
# for manual review of rejected calls only; nothing downstream ingests them. `*_mut_counts_not_called` and
# `pgdd_cohort_gd_summary.tsv` are not mocked (not needed for now, per Piotr 2026-11-06).
files <- c(GDs_per_tumour = "pgdd_GDs_per_tumour.tsv",
           GDs_per_region = "pgdd_GDs_per_region.tsv",
           GDs_events = "pgdd_GD_events.tsv",
           mut_counts = "pgdd_mut_counts.tsv",
           unstable_calls_GDs_per_tumour = "unstable_calls_pgdd_GDs_per_tumour.tsv",
           unstable_calls_GDs_per_region = "unstable_calls_pgdd_GDs_per_region.tsv",
           unstable_calls_GDs_events = "unstable_calls_pgdd_GD_events.tsv",
           unstable_calls_mut_counts = "unstable_calls_pgdd_mut_counts.tsv",
           bootstrap_discovery = "pgdd_bootstrap_discovery_summary.tsv",
           bootstrap_discovery_settings = "pgdd_bootstrap_discovery_settings.tsv")
for (el in names(files)) {
  fwrite(out[[el]], file.path(out_dir, files[[el]]), sep = "\t")
}
message("Wrote ", length(files), " PGDD mock files to ", out_dir)
