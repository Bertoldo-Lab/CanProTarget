#!/usr/bin/env Rscript
# Recover the deployable CanProTarget data products from the original local
# ChemProTarget data directory. Raw DepMap and chemoproteomics files remain
# outside the application repository; only compact, app-ready RDS files and
# subtype lists are written under data/.

suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
})

args <- commandArgs(trailingOnly = TRUE)
source_dir <- if (length(args) >= 1L) args[[1L]] else Sys.getenv("CPT_ORIGINAL_DATA_DIR")
output_dir <- if (length(args) >= 2L) args[[2L]] else "data"

if (!nzchar(source_dir) || !dir.exists(source_dir)) {
  stop(
    "Provide the original ChemProTarget data directory as the first argument, e.g.\n",
    "  Rscript docs/scripts/recover_original_data.R /path/to/ChemProTarget/data",
    call. = FALSE
  )
}
if (!dir.exists("R")) {
  stop("Run this script from the CanProTarget project root.", call. = FALSE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
source(file.path("R", "app_config.R"), local = FALSE)
source(file.path("R", "functions.R"), local = FALSE)
source(file.path("docs", "scripts", "reference_swissadme_from_csv.R"), local = FALSE)

need <- function(name) {
  path <- file.path(source_dir, name)
  if (!file.exists(path)) stop("Missing original source file: ", path, call. = FALSE)
  path
}

save_matrix <- function(source_name, output_name) {
  source_path <- need(source_name)
  output_path <- file.path(output_dir, output_name)
  sidecar_path <- paste0(tools::file_path_sans_ext(output_path), "_modelids.rds")
  cat("Reading ", source_name, " ...\n", sep = "")
  matrix <- utils::read.csv(
    source_path,
    row.names = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  colnames(matrix) <- sub(" \\([0-9]+\\)$", "", colnames(matrix))
  cat("  ", nrow(matrix), " models x ", ncol(matrix), " genes\n", sep = "")
  if (identical(output_name, "CRISPRGeneEffect_23Q4_clean.rds")) {
    last_gene <- colnames(matrix)[ncol(matrix)]
    if (identical(last_gene, "TNFRSF10C") || ncol(matrix) != 18443L || nrow(matrix) != 1100L) {
      stop(
        "CRISPR 23Q4 matrix looks truncated or the wrong file: ",
        nrow(matrix), " x ", ncol(matrix), ", last gene ", last_gene,
        call. = FALSE
      )
    }
  }
  # gzip, not xz. These matrices are read at runtime whenever the Target tab
  # draws its per-cell-line boxes, and xz costs 4 s for the RNAi matrix and
  # 8.5 s for CRISPR against 0.4 s and 0.6 s for gzip. The saving is 3 MB and
  # 6 MB respectively, which is not worth several seconds on every first use.
  saveRDS(matrix, output_path, compress = "gzip")
  saveRDS(rownames(matrix), sidecar_path, compress = "gzip")
  rm(matrix)
  invisible(gc())
}

save_matrix("CRISPRGeneEffect_23Q4.csv", "CRISPRGeneEffect_23Q4_clean.rds")
save_matrix("d2_gene_effect_headers_refined.csv", "d2_gene_effect_headers_refined.rds")

cat("Recovering DepMap model metadata and subtype lists ...\n")
if (!file.exists(file.path(source_dir, "Model_v2.csv"))) {
  stop(
    "Need Model_v2.csv in the source directory.\n",
    "Do not use the current Portal Model.csv; DepMap retracted it because it ",
    "includes models that do not exist.",
    call. = FALSE
  )
}
model <- utils::read.csv(need("Model_v2.csv"), stringsAsFactors = FALSE, check.names = FALSE)
saveRDS(model, file.path(output_dir, "cancer_model_data.rds"), compress = "xz")

for (dataset_label in names(app_config$file_map)) {
  rds_path <- file.path(output_dir, app_config$file_map[[dataset_label]])
  model_ids <- readRDS(paste0(tools::file_path_sans_ext(rds_path), "_modelids.rds"))
  choices <- cancer_subtype_choices_build(
    model,
    model_ids,
    output_dir,
    dataset_short_name(dataset_label)
  )
  list_path <- file.path(output_dir, app_config$cancer_subtype_list_files[[dataset_label]])
  writeLines(
    c(
      "# OncotreeSubtype labels available in this gene-effect matrix.",
      "# Recovered from the original CanProTarget data sources.",
      choices
    ),
    list_path,
    useBytes = TRUE
  )
  cat("  ", dataset_label, ": ", length(choices), " selectable subtypes\n", sep = "")
}

cat("Recovering the original chemoproteomics binding lookup ...\n")
binding <- readRDS(need("protein_binding_preprocessed.rds")) %>%
  dplyr::mutate(
    probe_name = toupper(as.character(.data$probe_name)),
    CR = suppressWarnings(as.numeric(.data$CR)),
    ligandable = tolower(as.character(.data$ligandable))
  )

target_counts <- binding %>%
  dplyr::filter(!is.na(.data$CR), .data$CR >= 4) %>%
  dplyr::group_by(.data$probe_name) %>%
  dplyr::summarise(n_targets = dplyr::n_distinct(.data$proteinid), .groups = "drop")

# Restore compound metadata using the same legacy probe renumbering that was
# used when the binding RDS was produced.
normalise_compound_name <- function(name) {
  name <- toupper(as.character(name))
  number <- suppressWarnings(as.integer(sub("^(ACRYL|CL)_", "", name)))
  is_acryl <- grepl("^ACRYL_[0-9]+$", name)
  is_cl <- grepl("^CL_[0-9]+$", name)
  out <- name
  out[is_acryl & number >= 1 & number <= 3] <- paste0("ACRYL_", number[is_acryl & number >= 1 & number <= 3] - 1)
  out[is_acryl & number >= 5 & number <= 11] <- paste0("ACRYL_", number[is_acryl & number >= 5 & number <= 11] - 2)
  out[is_acryl & number == 12] <- "OTHER_6"
  out[is_acryl & number == 13] <- "OTHER_7"
  out[is_acryl & number >= 14] <- paste0("ACRYL_", number[is_acryl & number >= 14] - 4)
  out[is_cl & number >= 20] <- paste0("CL_", number[is_cl & number >= 20] - 1)
  out
}

compound <- readxl::read_excel(need("table-s2.xlsx"), sheet = "Compound Keys") %>%
  dplyr::transmute(
    probe_name = normalise_compound_name(.data$Compound_Name),
    Dataset = as.character(.data$Dataset),
    Cell_Line = as.character(.data$Cell_Line),
    SMILES = as.character(.data$SMILES)
  ) %>%
  dplyr::group_by(.data$probe_name) %>%
  dplyr::summarise(
    Dataset = dplyr::first(stats::na.omit(.data$Dataset), default = NA_character_),
    Cell_Line = dplyr::first(stats::na.omit(.data$Cell_Line), default = NA_character_),
    SMILES = dplyr::first(stats::na.omit(.data$SMILES), default = NA_character_),
    .groups = "drop"
  )

binding <- binding %>%
  dplyr::left_join(target_counts, by = "probe_name") %>%
  dplyr::left_join(compound, by = "probe_name") %>%
  dplyr::mutate(
    n_targets = dplyr::coalesce(.data$n_targets, 0L),
    targets_total = .data$n_targets,
    gene_name_key = tolower(as.character(.data$gene_name))
  ) %>%
  dplyr::select(
    probe_name, proteinid, gene_name, CR,
    n_targets, targets_total, cysteineid,
    ligandable, Dataset, Cell_Line, SMILES,
    gene_name_key
  )

saveRDS(
  binding,
  file.path(output_dir, "protein_binding_lookup_preprocessed.rds"),
  compress = "xz"
)
cat(
  "  ", nrow(binding), " binding rows; ",
  dplyr::n_distinct(binding$gene_name), " genes; ",
  dplyr::n_distinct(binding$cysteineid), " cysteine sites\n",
  sep = ""
)

cat("Recovering SwissADME reference data ...\n")
swiss <- read_swissadme_from_raw_files(source_dir)
saveRDS(swiss, file.path(output_dir, "swissadme_preprocessed.rds"), compress = "xz")
cat("  ", nrow(swiss), " SwissADME rows\n", sep = "")

cat("Recovery complete. App-ready files are under ", output_dir, ".\n", sep = "")
