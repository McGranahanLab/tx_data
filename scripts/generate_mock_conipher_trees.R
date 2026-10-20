#!/usr/bin/env Rscript
# Generate a synthetic `tx842_trees.RDS` under `nemo_mock/`, matching the
# per-patient bundle shape observed in the real Tx842 release
# (`/nemo/project/proj-tracerx-lung/tctProjects/lungTx/Tx842/release/tracerx_842/tx842_trees.RDS`).
#
# Shape (per real example at `temp_data/tree_example.RDS`):
#   - top level: UNNAMED list; one element per patient-tumour
#   - each element: NAMED list, 15 slots; the two the ALPACA figure
#     pipeline consumes are:
#       $nested_pyclone$ccf_cluster_table      (matrix clones x regions, int 0..100+)
#       $clonality_out$clonality_table_corrected  (data.frame clones x regions, values {clonal, subclonal, absent})
#     with row names like "1","2","3",... (NOT "clone1") and column
#     names like "LTX0001_SU_T1.R1..6ecf51f11b1b".
#   - tumour identity recoverable from $parameters$sampleID (e.g. "LTX0001_1").
#
# This script reads `nemo_mock/_cohort_manifest.json` (produced by
# `generate_mocks.py`) so the mock cohort shape stays in lockstep with the
# other mock tables. Only `resolved` tumours are written; unresolved and
# needs_follow_up tumours are skipped (CONIPHER would not have produced an
# output for them in practice).
#
# Run:
#   Rscript scripts/generate_mock_conipher_trees.R
#
# Deterministic: fixed seed = manifest$random_seed.

suppressPackageStartupMessages({
  library(jsonlite)
})

# Resolve script dir from the `--file=` arg Rscript injects into commandArgs().
.script_dir <- (function() {
  all_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", all_args[grepl("^--file=", all_args)])
  if (length(file_arg) == 1 && nzchar(file_arg)) {
    return(normalizePath(dirname(file_arg)))
  }
  getwd()
})()

repo_root <- normalizePath(file.path(.script_dir, ".."))
mock_root <- file.path(repo_root, "nemo_mock")
manifest_path <- file.path(mock_root, "_cohort_manifest.json")
out_path <- file.path(mock_root, "tx842_trees.RDS")

if (!file.exists(manifest_path)) {
  stop(sprintf(
    "cohort manifest not found: %s\nRun `python scripts/generate_mocks.py` first.",
    manifest_path
  ))
}

manifest <- fromJSON(manifest_path, simplifyVector = FALSE)
set.seed(manifest$random_seed)

# ---- helpers ----------------------------------------------------------------

rand_hash <- function(nchar = 12) {
  paste0(sample(c(0:9, letters[1:6]), nchar, replace = TRUE), collapse = "")
}

region_colnames <- function(patient_id, tumour_ordinal, n_samples) {
  # Mirrors the real convention: "<patient>_SU_T<n>.R<r>..<hash>"
  sapply(seq_len(n_samples), function(r) {
    sprintf("%s_SU_T%d.R%d..%s", patient_id, tumour_ordinal, r, rand_hash())
  })
}

clonality_from_ccf <- function(ccf_val) {
  # Match the empirical banding seen in the real example:
  #   CCF >= 90        -> "clonal"
  #   10 <= CCF < 90   -> "subclonal"
  #   else             -> "absent"
  if (is.na(ccf_val)) return(NA_character_)
  if (ccf_val >= 90) return("clonal")
  if (ccf_val >= 10) return("subclonal")
  "absent"
}

# Build a CCF matrix for a single tumour: clones x regions, integer-valued
# "percentage CCF" (0..100ish). Clone1 (trunk) is clonal everywhere; descendants
# scatter: present in a subset of regions, with mixed clonality.
build_ccf_matrix <- function(n_clones, n_samples) {
  m <- matrix(0L, nrow = n_clones, ncol = n_samples)
  # Trunk: present everywhere.
  m[1, ] <- as.integer(sample(90:105, n_samples, replace = TRUE))
  if (n_clones < 2) return(m)
  for (ci in seq.int(2, n_clones)) {
    # Each descendant present in a random subset of regions (1..n_samples).
    k <- sample(seq_len(n_samples), 1)
    regions <- sample(seq_len(n_samples), k)
    for (r in regions) {
      # ~40% clonal (>=90), ~45% subclonal (10..89), ~15% barely detectable (<10).
      bucket <- sample(c("clonal", "subclonal", "trace"), 1,
                       prob = c(0.40, 0.45, 0.15))
      m[ci, r] <- as.integer(switch(
        bucket,
        clonal    = sample(90:105, 1),
        subclonal = sample(10:85, 1),
        trace     = sample(0:9, 1)
      ))
    }
  }
  m
}

build_clonality_df <- function(ccf_mat) {
  as.data.frame(
    apply(ccf_mat, c(1, 2), clonality_from_ccf),
    stringsAsFactors = FALSE
  )
}

# A per-patient bundle. Only the slots the figure pipeline reads carry real
# content; the rest are structurally present but minimally stubbed.
build_bundle <- function(spec) {
  patient <- spec$patient_id
  t_ord   <- spec$tumour_ordinal
  n_s     <- spec$n_samples
  n_c     <- length(spec$clones)

  sample_cols <- region_colnames(patient, t_ord, n_s)
  clone_rows  <- as.character(seq_len(n_c))

  ccf_mat <- build_ccf_matrix(n_c, n_s)
  rownames(ccf_mat) <- clone_rows
  colnames(ccf_mat) <- sample_cols

  clonality_df <- build_clonality_df(ccf_mat)
  rownames(clonality_df) <- clone_rows
  colnames(clonality_df) <- sample_cols

  # Short tumour ID in LTX _N form (parameters$sampleID flavour).
  short_id <- sprintf("%s_%d", patient, t_ord)
  long_id  <- sprintf("%s0_%d", patient, t_ord)  # extra zero mirrors real "LTX00001_1"

  # Minimal but shape-faithful stubs for the non-critical slots.
  empty_ccf_table <- data.frame(matrix(nrow = 0, ncol = 0))

  # graph_pyclone$edgelength: a 1-D `table` object keyed by clone id, carrying
  # the SNV count on the branch into each clone. Must be a `table` (not a
  # plain named vector) because analysis_functions.R::get_clone_genomic_distances
  # calls `as.data.table(edgelength)` and expects a 2-column result (V1=clone,
  # N=count) — `as.data.table` collapses a named vector to one column but
  # unpacks a `table` to two. Trunk gets ~500-700 muts, descendants 5-60.
  edgelength <- as.table(stats::setNames(
    as.integer(c(sample(500:700, 1), sample(5:60, n_c - 1, replace = TRUE))),
    clone_rows
  ))

  list(
    ccf_table_pyclone              = empty_ccf_table,
    ccf_table_absolute             = empty_ccf_table,
    ccf_table_pyclone_clean        = empty_ccf_table,
    ccf_table_absolute_clean       = empty_ccf_table,
    merged_clusters                = FALSE,
    noisy_clusters_adjusted        = FALSE,
    cpn_removed_clusters           = character(0),
    tree_removed_clusters          = character(0),
    graph_pyclone                  = list(
      sampleID       = short_id,
      long_sampleID  = long_id,
      edgelength     = edgelength
    ),
    parameters                     = list(
      sampleID = short_id,
      prefix   = sub("[0-9].*", "", patient)   # "LTX" / "CRUK"
    ),
    nested_pyclone                 = list(
      ccf_cluster_table = ccf_mat
    ),
    clonality_table                = clonality_df,         # original (same as corrected for the mock)
    clonality_out                  = list(
      clonality_table_corrected = clonality_df,
      clonality_table_original  = clonality_df
    ),
    clone_proportion_out           = list(),
    subclonal_expansion_score_out  = list()
  )
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- build the cohort list --------------------------------------------------

resolved <- Filter(function(x) identical(x$status, "resolved"), manifest$tumours)
cat(sprintf(
  "generating conipher_trees mock: %d resolved tumours (of %d total)\n",
  length(resolved), length(manifest$tumours)
))

trees <- lapply(resolved, build_bundle)   # unnamed list — matches real file

if (!dir.exists(mock_root)) dir.create(mock_root, recursive = TRUE)
saveRDS(trees, file = out_path)

cat(sprintf("wrote %s (%d bundles)\n", out_path, length(trees)))
cat("done.\n")
