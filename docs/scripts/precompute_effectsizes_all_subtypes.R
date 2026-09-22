#!/usr/bin/env Rscript
# ============================================================
# Offline: one TSV per feasible OncotreeSubtype (same MIN_GE_CELL_LINES as the app).
# Run from project root:
#   Rscript docs/scripts/precompute_effectsizes_all_subtypes.R
#
# Writes: data/precomputed_effectsizes/<CRISPR|RNAi>_<sanitized_subtype>.rds
#         and optional per-dataset logs: precompute_log_<dataset>.csv (for your records only;
#         the Shiny app does not read these logs — it indexes precomputed_effectsizes/).
#
# Needs:  data/*gene effect*.rds, data/cancer_model_data.rds
#         Same R packages as the Shiny app (cdsrmodels, dplyr, limma, …).
# See docs/DATA_PROVENANCE.md.
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(magrittr)  # cdsrmodels uses %<>%
  library(tibble)    # cdsrmodels uses rownames_to_column
  library(cdsrmodels)
})

data_dir <- "data"
source(file.path("R", "functions.R"), local = FALSE)
# sanitize_subtype() comes from R/functions.R (cpt_sanitize_subtype_key alias)

if (!dir.exists("R") || !dir.exists(data_dir)) {
  stop("Run from CanProTarget project root", call. = FALSE)
}

out_dir <- file.path(data_dir, "precomputed_effectsizes")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

datasets <- c(
  "CRISPR (23Q4)" = "CRISPRGeneEffect_23Q4_clean.rds",
  "RNAi"          = "d2_gene_effect_headers_refined.rds"
)

meta_path_rds <- file.path(data_dir, "cancer_model_data.rds")
if (!file.exists(meta_path_rds)) {
  stop(
    "Missing ", meta_path_rds,
    ". Build it with docs/scripts/rebuild_depmap_23q4.R (needs Model_v2.csv). ",
    "See docs/DATA_PROVENANCE.md.",
    call. = FALSE
  )
}
meta <- readRDS(meta_path_rds)

for (ds_label in names(datasets)) {
  rds_fn <- datasets[[ds_label]]
  rds_path <- file.path(data_dir, rds_fn)
  if (!file.exists(rds_path)) {
    cat("SKIP (no file):", rds_fn, "\n")
    next
  }

  ds_short <- if (grepl("CRISPR", ds_label)) "CRISPR" else "RNAi"

  cat(sprintf("\n=== %s ===\n", ds_label))

  mat <- readRDS(rds_path)
  ids <- rownames(mat)
  subs <- feasible_oncotree_subtypes(meta, ids, min_lines = MIN_GE_CELL_LINES)
  subs <- sort(unique(subs))
  cat("Subtypes:", length(subs), "\n")

  log_rows <- list()

  for (st in subs) {
    fname <- paste0(ds_short, "_", sanitize_subtype(st), ".rds")
    fpath <- file.path(out_dir, fname)

    ok <- tryCatch(
      {
        pval_res <- lin_ass_pval(mat, meta, st)
        lm_res   <- run_lm_ge(mat, meta, st)
        out <- ge_analysis(
          pval_results              = pval_res,
          lm_results                = lm_res,
          apply_pvalue_filter       = FALSE,
          min_effect_size           = 0,
          exclude_common_essentials = FALSE,
          filter_pval_vs_noncancer  = FALSE,
          effect_size_matrix        = mat,
          cancer_model_data         = meta,
          selected_subtype          = st
        )
        saveRDS(out$all_gene_ge_df, fpath)
        cat("  OK:", fname, "\n")
        TRUE
      },
      error = function(e) {
        cat("  FAIL:", st, "—", conditionMessage(e), "\n")
        FALSE
      }
    )

    n_in <- sum(
      meta$OncotreeSubtype %in% st &
        meta$ModelID %in% ids &
        !is.na(meta$OncotreeSubtype)
    )
    log_rows[[length(log_rows) + 1]] <- data.frame(
      dataset = ds_short,
      subtype = st,
      n_lines = as.integer(n_in),
      n_genes = if (ok && file.exists(fpath)) {
        length(readr::read_lines(fpath)) - 1L
      } else {
        NA_integer_
      },
      file    = fname,
      status  = if (ok) "SUCCESS" else "FAIL",
      stringsAsFactors = FALSE
    )
  }

  rm(mat)
  gc()

  if (length(log_rows)) {
    log_df <- do.call(rbind, log_rows)
    utils::write.csv(
      log_df,
      file.path(out_dir, paste0("precompute_log_", ds_short, ".csv")),
      row.names = FALSE
    )
  }
}

cat("\nDone. Output:", normalizePath(out_dir), "\n")
