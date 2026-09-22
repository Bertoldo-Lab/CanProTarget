#!/usr/bin/env Rscript
# Rebuild DepMap 23Q4 CRISPR + Model_v2 products for CanProTarget.
#
# Run from the project root:
#   Rscript docs/scripts/rebuild_depmap_23q4.R
#
# Expected inputs (gitignored; do not commit the CSVs):
#   data/CRISPRGeneEffect_23Q4.csv   DepMap Public 23Q4 gene effect
#   data/Model_v2.csv                DepMap 23Q4 Model_v2 (Portal Model.csv is retracted)
#
# Writes:
#   data/CRISPRGeneEffect_23Q4_clean.rds
#   data/CRISPRGeneEffect_23Q4_clean_modelids.rds
#   data/cancer_model_data.rds
#   data/cancer_subtypes_CRISPR.txt
#   data/cancer_subtypes_RNAi.txt
#
# Never open the CRISPR CSV in Excel. Excel's 16,384-column limit truncates
# the matrix at TNFRSF10C and silently drops TP53 and ~2,000 other genes.

CRISPR_EXPECTED_N_MODELS <- 1100L
CRISPR_EXPECTED_N_GENES  <- 18443L
CRISPR_EXPECTED_LAST     <- "ZZZ3"

strip_depmap_entrez <- function(names) {
  sub(" \\([0-9]+\\)$", "", names)
}

find_crispr_csv <- function(data_dir) {
  candidates <- c(
    file.path(data_dir, "CRISPRGeneEffect_23Q4.csv"),
    file.path(data_dir, "CRISPRGeneEffect_23Q4_clean.csv")
  )
  found <- candidates[file.exists(candidates)]
  if (!length(found)) {
    stop(
      "Missing CRISPR CSV. Place DepMap Public 23Q4 CRISPRGeneEffect.csv at\n",
      "  ", candidates[[1]], "\n",
      "Do not open or save that file in Excel.",
      call. = FALSE
    )
  }
  found[[1]]
}

find_model_csv <- function(data_dir) {
  v2 <- file.path(data_dir, "Model_v2.csv")
  if (file.exists(v2)) {
    return(v2)
  }
  stop(
    "Missing data/Model_v2.csv.\n",
    "Use the 23Q4 Model_v2 download. Do not use the current Portal Model.csv;\n",
    "DepMap retracted it because it includes models that do not exist.",
    call. = FALSE
  )
}

assert_crispr_23q4 <- function(mat) {
  last_gene <- colnames(mat)[ncol(mat)]
  if (identical(last_gene, "TNFRSF10C")) {
    stop(
      "CRISPR matrix ends at TNFRSF10C (", ncol(mat), " genes). ",
      "That is the Excel 16,384-column truncation. Re-download 23Q4 and ",
      "read it in R without opening it in Excel.",
      call. = FALSE
    )
  }
  if (nrow(mat) != CRISPR_EXPECTED_N_MODELS || ncol(mat) != CRISPR_EXPECTED_N_GENES) {
    stop(
      "Unexpected CRISPR 23Q4 shape: ", nrow(mat), " x ", ncol(mat),
      " (expected ", CRISPR_EXPECTED_N_MODELS, " x ", CRISPR_EXPECTED_N_GENES, ").",
      call. = FALSE
    )
  }
  if (!"TP53" %in% colnames(mat)) {
    stop("TP53 is missing after Entrez stripping; gene names look wrong.", call. = FALSE)
  }
  if (!identical(last_gene, CRISPR_EXPECTED_LAST)) {
    stop(
      "Last CRISPR gene is ", last_gene, " (expected ", CRISPR_EXPECTED_LAST, ").",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

read_crispr_23q4 <- function(csv_path) {
  cat("Reading ", csv_path, " ...\n", sep = "")
  mat <- utils::read.csv(
    csv_path,
    row.names = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  colnames(mat) <- strip_depmap_entrez(colnames(mat))
  if (anyDuplicated(colnames(mat))) {
    stop("Duplicate gene symbols after stripping Entrez IDs.", call. = FALSE)
  }
  assert_crispr_23q4(mat)
  cat(
    "  ", nrow(mat), " models x ", ncol(mat), " genes; last=",
    colnames(mat)[ncol(mat)], "; TP53=TRUE; ACH-002509=",
    "ACH-002509" %in% rownames(mat), "\n",
    sep = ""
  )
  mat
}

rebuild_depmap_23q4 <- function(data_dir = "data") {
  if (!dir.exists("R") || !dir.exists(data_dir)) {
    stop("Run from the CanProTarget project root (folder with data/ and R/).", call. = FALSE)
  }
  source(file.path("R", "app_config.R"), local = FALSE)
  source(file.path("R", "functions.R"), local = FALSE)

  csv_path <- find_crispr_csv(data_dir)
  rds_path <- file.path(data_dir, "CRISPRGeneEffect_23Q4_clean.rds")
  mid_path <- file.path(data_dir, "CRISPRGeneEffect_23Q4_clean_modelids.rds")
  mat <- read_crispr_23q4(csv_path)
  cat("Writing ", rds_path, " (xz) ...\n", sep = "")
  # gzip, not xz. These matrices are read at runtime whenever the Target tab
  # draws its per-cell-line boxes, and xz costs 4 s for the RNAi matrix and
  # 8.5 s for CRISPR against 0.4 s and 0.6 s for gzip. The saving is 3 MB and
  # 6 MB respectively, which is not worth several seconds on every first use.
  saveRDS(mat, rds_path, compress = "gzip")
  saveRDS(rownames(mat), mid_path, compress = "gzip")

  model_path <- find_model_csv(data_dir)
  cat("Reading ", model_path, " ...\n", sep = "")
  meta <- utils::read.csv(model_path, stringsAsFactors = FALSE, check.names = FALSE)
  model_out <- file.path(data_dir, "cancer_model_data.rds")
  saveRDS(meta, model_out, compress = "xz")
  cat("  ", nrow(meta), " models x ", ncol(meta), " columns\n", sep = "")

  n_overlap <- sum(rownames(mat) %in% meta$ModelID)
  cat("  CRISPR ModelIDs in Model_v2: ", n_overlap, "/", nrow(mat), "\n", sep = "")
  am <- meta$ModelID[meta$OncotreeSubtype == "Acral Melanoma" & !is.na(meta$OncotreeSubtype)]
  cat(
    "  Acral Melanoma in metadata: ", length(am),
    "; in CRISPR matrix: ", sum(am %in% rownames(mat)), "\n",
    sep = ""
  )

  rm(mat)
  invisible(gc())

  cat("Writing subtype lists ...\n")
  for (ds_label in names(app_config$file_map)) {
    rds_name <- app_config$file_map[[ds_label]]
    txt_name <- app_config$cancer_subtype_list_files[[ds_label]]
    if (is.null(txt_name) || !nzchar(txt_name)) {
      next
    }
    ge_path <- file.path(data_dir, rds_name)
    txt_path <- file.path(data_dir, txt_name)
    if (!file.exists(ge_path)) {
      cat("  SKIP ", ds_label, ": ", rds_name, " not found\n", sep = "")
      next
    }
    ge_mid <- paste0(tools::file_path_sans_ext(ge_path), "_modelids.rds")
    ids <- if (file.exists(ge_mid)) {
      readRDS(ge_mid)
    } else {
      rownames(readRDS(ge_path))
    }
    short <- dataset_short_name(ds_label)
    ch <- cancer_subtype_choices_build(meta, ids, data_dir, short)
    hdr <- c(
      "# OncotreeSubtype labels for the Cancer Dependencies picker.",
      "# Regenerated from Model_v2 intersected with this gene-effect matrix."
    )
    writeLines(c(hdr, ch), txt_path, useBytes = TRUE)
    cat("  ", ds_label, ": ", length(ch), " -> ", basename(txt_path), "\n", sep = "")
  }

  invisible(TRUE)
}

if (sys.nframe() == 0L) {
  rebuild_depmap_23q4("data")
  cat("Done.\n")
}
