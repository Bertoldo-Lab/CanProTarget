#!/usr/bin/env Rscript

# Build the compact CanProTarget functional-cysteine atlas from a local clone of
# https://github.com/cravattlab/Cys_editing.
#
# Usage (from the CanProTarget project root):
#   Rscript docs/scripts/import_cys_editing_atlas.R ../Cys_editing_reference

args <- commandArgs(trailingOnly = TRUE)
source_root <- if (length(args)) args[[1]] else "../Cys_editing_reference"
output_path <- if (length(args) >= 2L) args[[2]] else file.path("data", "cys_editing_atlas.rds")

source("R/cys_editing_functions.R")

load_object <- function(path, expected_name) {
  if (!file.exists(path)) {
    stop("Missing Cys_editing source file: ", path, call. = FALSE)
  }
  environment <- new.env(parent = emptyenv())
  loaded_names <- load(path, envir = environment)
  if (!expected_name %in% loaded_names) {
    stop(
      basename(path), " does not contain expected object ", expected_name, ".",
      call. = FALSE
    )
  }
  environment[[expected_name]]
}

rdat_dir <- file.path(source_root, "Part5_global_analysis", "Rdat")
cys_dropout <- load_object(file.path(rdat_dir, "Cys_dropout.Rdat"), "Cys_dropout")
dep_score <- load_object(file.path(rdat_dir, "ceres_used.Rdat"), "ceres_used")
af_rsa <- load_object(file.path(rdat_dir, "AF_RSA.Rdat"), "AF_RSA")

# Conservation, ortholog, and ClinVar data (optional but recommended)
conserve_score <- tryCatch(
  load_object(file.path(rdat_dir, "Cys_conserve_score.Rdat"), "Cys_conserve_score"),
  error = function(e) {
    message("Note: Cys_conserve_score.Rdat not found, skipping conservation data.")
    NULL
  }
)
ortho_count <- tryCatch(
  load_object(file.path(rdat_dir, "ortho_count_df_20230419.Rdat"), "ortho_count_df"),
  error = function(e) {
    message("Note: ortho_count_df_20230419.Rdat not found, skipping ortholog data.")
    NULL
  }
)
clinvar <- tryCatch(
  load_object(file.path(rdat_dir, "clinvar_data.Rdat"), "clinvar_data"),
  error = function(e) {
    message("Note: clinvar_data.Rdat not found, skipping ClinVar data.")
    NULL
  }
)

source_commit <- tryCatch(
  system2("git", c("-C", shQuote(source_root), "rev-parse", "HEAD"), stdout = TRUE),
  warning = function(e) NA_character_,
  error = function(e) NA_character_
)
source_commit <- if (length(source_commit)) source_commit[[1]] else NA_character_

atlas <- cpt_cys_build_atlas(
  cys_dropout = cys_dropout,
  dep_score = dep_score,
  af_rsa = af_rsa,
  conserve_score = conserve_score,
  ortho_count = ortho_count,
  clinvar = clinvar,
  source_commit = source_commit
)
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(atlas, output_path, compress = "xz")

cat(
  "Wrote ", output_path, "\n",
  "  sites: ", nrow(atlas), "\n",
  "  genes: ", length(unique(atlas$gene_symbol)), "\n",
  "  functional sites: ", sum(atlas$functional), "\n",
  "  functional + ligandable sites: ", sum(atlas$functional_ligandable), "\n",
  "  with conservation data: ", sum(!is.na(atlas$conservation_score)), "\n",
  "  with ortholog data: ", sum(!is.na(atlas$ortholog_cys_count)), "\n",
  "  with ClinVar pathogenic: ", sum(atlas$clinvar_pathogenic), "\n",
  sep = ""
)
