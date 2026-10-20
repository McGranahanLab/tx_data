#!/usr/bin/env Rscript
# Dump the aggregated `conipher_trees` RDS to two long TSVs so the Python
# builder can canonicalise them into Parquet without an R dependency on the
# tx_data env itself. Called by `tx_data.builds.conipher_trees.build()`.
#
# Flags:
#   --input=<path>       aggregated tx842_trees.RDS
#   --output-dir=<path>  directory to write ccf.tsv + clonality.tsv into
#
# Output schema (both TSVs):
#   source_tumour_id   TAB   source_clone   TAB   region   TAB   value
#     - source_tumour_id : raw $parameters$sampleID (e.g. "LTX0001_1")
#     - source_clone     : raw rowname of the matrix (e.g. "1")
#     - region           : raw colname (e.g. "LTX0001_SU_T1.R1..6ecf51f11b1b")
#     - value            : integer CCF (0..~105) OR clonality string
# The Python builder owns canonicalisation of tumour_id / clone and any
# region-id reshaping decisions. We stay faithful to source here.

parse_flags <- function() {
  raw <- commandArgs(trailingOnly = TRUE)
  flags <- list()
  for (a in raw) {
    if (grepl("^--input=", a))           flags$input      <- sub("^--input=", "", a)
    else if (grepl("^--output-dir=", a)) flags$output_dir <- sub("^--output-dir=", "", a)
    else stop(sprintf("Unknown flag: %s", a))
  }
  for (k in c("input", "output_dir")) {
    if (is.null(flags[[k]])) stop(sprintf("Missing required flag: --%s=<path>", gsub("_", "-", k)))
  }
  flags
}

long_from_matrix_like <- function(tbl, source_tumour_id) {
  # Works for matrix OR data.frame; both have row/col names.
  mat <- as.matrix(tbl)
  rn  <- rownames(mat)
  cn  <- colnames(mat)
  if (is.null(rn) || is.null(cn)) {
    stop(sprintf("[dump_conipher] %s: table missing row or col names", source_tumour_id))
  }
  # Produce long form via `expand.grid` + `as.vector(t(mat))` (row-major).
  grid <- expand.grid(source_clone = rn, region = cn,
                      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  grid$value <- as.vector(mat)
  grid$source_tumour_id <- source_tumour_id
  grid[, c("source_tumour_id", "source_clone", "region", "value")]
}

main <- function() {
  flags <- parse_flags()
  trees <- readRDS(flags$input)
  if (!is.list(trees)) stop("aggregated RDS is not a list")

  dir.create(flags$output_dir, showWarnings = FALSE, recursive = TRUE)
  ccf_parts <- list()
  cln_parts <- list()
  n_bundles <- length(trees)
  n_dumped  <- 0L

  for (i in seq_along(trees)) {
    b <- trees[[i]]
    tum <- b$parameters$sampleID
    if (is.null(tum) || !nzchar(tum)) {
      message(sprintf("[dump_conipher] WARN bundle %d: missing parameters$sampleID — skipping", i))
      next
    }
    ccf <- b$nested_pyclone$ccf_cluster_table
    cln <- b$clonality_out$clonality_table_corrected
    if (is.null(ccf) || is.null(cln)) {
      message(sprintf(
        "[dump_conipher] WARN bundle %d (%s): missing ccf_cluster_table or clonality_table_corrected — skipping",
        i, tum))
      next
    }
    ccf_parts[[length(ccf_parts) + 1]] <- long_from_matrix_like(ccf, tum)
    cln_parts[[length(cln_parts) + 1]] <- long_from_matrix_like(cln, tum)
    n_dumped <- n_dumped + 1L
  }

  ccf_long <- do.call(rbind, ccf_parts)
  cln_long <- do.call(rbind, cln_parts)

  # Numeric cast for CCF; keep clonality as character.
  ccf_long$value <- as.integer(ccf_long$value)
  cln_long$value <- as.character(cln_long$value)

  ccf_out <- file.path(flags$output_dir, "conipher_ccf.tsv")
  cln_out <- file.path(flags$output_dir, "conipher_clonality.tsv")
  write.table(ccf_long, file = ccf_out, sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(cln_long, file = cln_out, sep = "\t", row.names = FALSE, quote = FALSE)

  cat(sprintf("[dump_conipher] dumped %d/%d bundles; ccf=%d rows -> %s; clonality=%d rows -> %s\n",
              n_dumped, n_bundles, nrow(ccf_long), ccf_out, nrow(cln_long), cln_out))
}

invisible(main())
