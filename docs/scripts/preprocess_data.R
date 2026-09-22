#!/usr/bin/env Rscript
# ============================================================
# Offline reference pipeline — NOT used by the Shiny app.
# Run from project root (directory containing data/, R/, docs/):
#   Rscript docs/scripts/preprocess_data.R
#
# See docs/DATA_PROVENANCE.md for what each output is and which raw files it came from.
#
# Produces / refreshes under data/:
#   - Gene effect *.rds (from DepMap-style CSV matrices)
#   - swissadme_preprocessed.rds (from swissadme.csv + optional Backus link)
#   - cancer_model_data.rds (from DepMap 23Q4 Model_v2.csv)
#   - cancer_subtypes_*.txt
#   - protein_binding_lookup_preprocessed.rds (via docs/scripts/preprocess_protein_binding.R)
# ============================================================

library(dplyr)
library(readxl)
library(tidyr)

if (!dir.exists("data") || !dir.exists("R")) {
  stop(
    "Run from CanProTarget project root (folder with data/ and R/), e.g.\n",
    "  cd cell_cpt && Rscript docs/scripts/preprocess_data.R",
    call. = FALSE
  )
}

data_dir <- "data"
source(file.path("R", "app_config.R"), local = FALSE)
source(file.path("R", "functions.R"), local = FALSE)
source(file.path("docs", "scripts", "reference_swissadme_from_csv.R"), local = FALSE)
source(file.path("docs", "scripts", "rebuild_depmap_23q4.R"), local = FALSE)

needs_rebuild <- function(out, srcs) {
  srcs <- srcs[file.exists(srcs)]
  if (!length(srcs)) {
    return(FALSE)
  }
  if (!file.exists(out)) {
    return(TRUE)
  }
  max(file.mtime(srcs), na.rm = TRUE) > file.mtime(out)
}

# --- 1. Gene effect matrices (CSV -> xz-compressed RDS) -----------------
cat("\n=== 1. Gene effect matrices ===\n")

# CRISPR 23Q4: official CSV, Entrez stripped, Excel truncation refused.
crispr_csv <- tryCatch(find_crispr_csv(data_dir), error = function(e) "")
crispr_rds <- file.path(data_dir, "CRISPRGeneEffect_23Q4_clean.rds")
crispr_mid <- file.path(data_dir, "CRISPRGeneEffect_23Q4_clean_modelids.rds")
if (!nzchar(crispr_csv)) {
  cat("  SKIP: CRISPRGeneEffect_23Q4.csv not found\n")
} else if (!needs_rebuild(crispr_rds, crispr_csv)) {
  cat("  OK (up to date): CRISPRGeneEffect_23Q4_clean.rds\n")
} else {
  cat("  Building: CRISPR (23Q4) from", basename(crispr_csv), "...\n")
  t0 <- Sys.time()
  mat <- read_crispr_23q4(crispr_csv)
  saveRDS(mat, crispr_rds, compress = "xz")
  saveRDS(rownames(mat), crispr_mid, compress = "xz")
  rm(mat)
  gc()
  cat("     done in", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "s\n")
}

gene_map <- list(
  list(
    label = "RNAi (D2 refined)",
    csv = "d2_gene_effect_headers_refined.csv",
    rds = "d2_gene_effect_headers_refined.rds"
  )
)

for (g in gene_map) {
  csv_path <- file.path(data_dir, g$csv)
  rds_path <- file.path(data_dir, g$rds)
  if (!file.exists(csv_path)) {
    cat("  SKIP:", g$csv, "not found\n")
    next
  }
  mid_path <- paste0(tools::file_path_sans_ext(rds_path), "_modelids.rds")
  if (!needs_rebuild(rds_path, csv_path)) {
    if (!file.exists(mid_path) && file.exists(rds_path)) {
      cat("  Writing model-ID sidecar:", basename(mid_path), "(one-time)\n")
      tmp <- readRDS(rds_path)
      saveRDS(rownames(tmp), mid_path, compress = "xz")
      rm(tmp)
      gc()
    }
    cat("  OK (up to date):", g$rds, "\n")
    next
  }
  cat("  Building:", g$label, "...\n")
  t0 <- Sys.time()
  mat <- utils::read.csv(csv_path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  cat("     ", nrow(mat), "x", ncol(mat), "— saving RDS (xz)...\n")
  saveRDS(mat, rds_path, compress = "xz")
  saveRDS(rownames(mat), mid_path, compress = "xz")
  rm(mat)
  gc()
  cat(
    "     done in",
    round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "s\n"
  )
}

# --- 2. SwissADME ------------------------------------------------------
cat("\n=== 2. SwissADME (swissadme_preprocessed.rds) ===\n")

backus_xlsx <- file.path(data_dir, "backus_swiss_link.xlsx")
backus_rds <- file.path(data_dir, "backus_swiss_link.rds")
if (file.exists(backus_xlsx)) {
  if (!file.exists(backus_rds) || file.mtime(backus_xlsx) > file.mtime(backus_rds)) {
    cat("  Refreshing backus_swiss_link.rds from xlsx...\n")
    saveRDS(readxl::read_excel(backus_xlsx), backus_rds, compress = "xz")
  }
}

swiss_out <- file.path(data_dir, "swissadme_preprocessed.rds")
swiss_src <- c(
  file.path(data_dir, "swissadme.csv"),
  file.path(data_dir, "backus_swiss_link.xlsx"),
  file.path(data_dir, "backus_swiss_link.rds")
)

if (!file.exists(file.path(data_dir, "swissadme.csv"))) {
  cat("  SKIP: swissadme.csv missing\n")
} else if (!needs_rebuild(swiss_out, swiss_src[file.exists(swiss_src)])) {
  cat("  OK (up to date): swissadme_preprocessed.rds\n")
} else {
  cat("  Building from CSV + optional link file...\n")
  adme <- read_swissadme_from_raw_files(data_dir)
  saveRDS(adme, swiss_out, compress = "xz")
  cat("  Rows:", nrow(adme), " cols:", ncol(adme), "\n")
}

# --- 3. Cancer model metadata ------------------------------------------
cat("\n=== 3. Cancer model (cancer_model_data.rds) ===\n")

model_csv <- tryCatch(find_model_csv(data_dir), error = function(e) "")
model_out <- file.path(data_dir, "cancer_model_data.rds")

if (!nzchar(model_csv)) {
  cat("  SKIP: Model_v2.csv missing (do not use retracted Portal Model.csv)\n")
} else if (!needs_rebuild(model_out, model_csv)) {
  cat("  OK (up to date): cancer_model_data.rds\n")
} else {
  cat("  Reading", basename(model_csv), "...\n")
  meta <- utils::read.csv(model_csv, stringsAsFactors = FALSE, check.names = FALSE)
  saveRDS(meta, model_out, compress = "xz")
  cat("  Rows:", nrow(meta), "\n")
}

# --- 3b. Cancer subtype picker lists (plain text per gene-effect dataset) -----
cat("\n=== 3b. Cancer subtype lists (cancer_subtypes_*.txt) ===\n")

if (!file.exists(model_out)) {
  cat("  SKIP: cancer_model_data.rds missing\n")
} else {
  meta <- readRDS(model_out)
  for (ds_label in names(app_config$file_map)) {
    rds_name <- app_config$file_map[[ds_label]]
    txt_name <- app_config$cancer_subtype_list_files[[ds_label]]
    if (is.null(txt_name) || !nzchar(txt_name)) {
      next
    }
    rds_path <- file.path(data_dir, rds_name)
    txt_path <- file.path(data_dir, txt_name)
    if (!file.exists(rds_path)) {
      cat("  SKIP:", ds_label, "-", basename(rds_path), "not found\n")
      next
    }
    mid_path <- paste0(tools::file_path_sans_ext(rds_path), "_modelids.rds")
    ids <- if (file.exists(mid_path)) {
      readRDS(mid_path)
    } else {
      rownames(readRDS(rds_path))
    }
    short <- dataset_short_name(ds_label)
    ch <- cancer_subtype_choices_build(meta, ids, data_dir, short)
    hdr <- c(
      "# OncotreeSubtype labels for the Cancer Dependencies picker.",
      "# Regenerated by docs/scripts/preprocess_data.R — edit only if you know the exact DepMap strings."
    )
    writeLines(c(hdr, ch), txt_path, useBytes = TRUE)
    cat("  Wrote", length(ch), "lines ->", basename(txt_path), "\n")
  }
}

# --- 4. Protein binding (lookup RDS) -----------------------------------
cat("\n=== 4. Protein binding (preprocess_protein_binding.R) ===\n")

pb_script <- file.path("docs", "scripts", "preprocess_protein_binding.R")
if (!file.exists(pb_script)) {
  stop("Missing ", pb_script, call. = FALSE)
}

pb_out <- file.path(data_dir, "protein_binding_lookup_preprocessed.rds")
pb_required <- c(
  file.path(data_dir, "table-s2.xlsx"),
  file.path(data_dir, "id_mapping.tsv")
)
pb_swiss <- file.path(data_dir, "swissadme_preprocessed.rds")
pb_src <- c(pb_required, pb_swiss)
pb_src <- pb_src[file.exists(pb_src)]

if (!all(file.exists(pb_required))) {
  cat("  SKIP: need both table-s2.xlsx and id_mapping.tsv\n")
} else if (!needs_rebuild(pb_out, pb_src)) {
  cat("  OK (up to date): protein_binding_lookup_preprocessed.rds\n")
} else {
  cat("  Running protein binding pipeline (may take a few minutes)...\n")
  source(pb_script, local = FALSE)
}

cat("\n=== Preprocessing finished ===\n")
