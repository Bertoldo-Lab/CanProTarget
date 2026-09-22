#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(rvest)
})

args <- commandArgs(trailingOnly = TRUE)
benchmark_dir <- if (length(args)) {
  normalizePath(args[[1]], mustWork = TRUE)
} else {
  normalizePath(file.path(getwd(), "manuscript_revision_2026-08-03", "software_benchmark"),
                mustWork = TRUE)
}

input_dir <- file.path(benchmark_dir, "data_external")
page_files <- file.path(input_dir, sprintf("covpdb_cysteine_page_%d.html", 1:6))
if (!all(file.exists(page_files))) {
  stop("Missing one or more CovPDB cysteine index pages.", call. = FALSE)
}

read_page <- function(path, page_number) {
  page <- read_html(path)
  table_node <- html_element(page, "table.result_table")
  if (is.na(table_node)) {
    stop("Could not find CovPDB result table in ", path, call. = FALSE)
  }
  tab <- as.data.frame(html_table(table_node, fill = TRUE), stringsAsFactors = FALSE)
  if (ncol(tab) != 5L) {
    stop("Unexpected CovPDB table structure in ", path, call. = FALSE)
  }
  names(tab) <- c("covpdb_index", "protein_name", "organism", "uniprot_label",
                  "n_protein_ligand_complexes")
  tab$source_page <- page_number
  tab
}

all_rows <- do.call(rbind, Map(read_page, page_files, seq_along(page_files)))
all_rows$uniprot_accession <- sub("\\s.*$", "", trimws(all_rows$uniprot_label))
all_rows$uniprot_accession[all_rows$uniprot_accession %in% c("", "None", "NA")] <- NA_character_
all_rows$protein_name <- trimws(all_rows$protein_name)
all_rows$organism <- trimws(all_rows$organism)
all_rows$n_protein_ligand_complexes <- suppressWarnings(
  as.integer(all_rows$n_protein_ligand_complexes)
)

all_rows <- all_rows[order(all_rows$covpdb_index), c(
  "covpdb_index", "protein_name", "organism", "uniprot_accession",
  "n_protein_ligand_complexes", "source_page"
)]
human_rows <- all_rows[grepl("^Homo sapiens", all_rows$organism) &
                         !is.na(all_rows$uniprot_accession), , drop = FALSE]

if (nrow(all_rows) != 291L) {
  warning("Expected 291 CovPDB cysteine proteins but parsed ", nrow(all_rows), ".")
}
if (anyDuplicated(all_rows$covpdb_index)) {
  stop("Duplicate CovPDB row indices after parsing.", call. = FALSE)
}

write.csv(all_rows, file.path(input_dir, "covpdb_cysteine_proteins_all.csv"),
          row.names = FALSE, na = "")
write.csv(human_rows, file.path(input_dir, "covpdb_cysteine_targets_human.csv"),
          row.names = FALSE, na = "")

cat(sprintf(
  "Parsed %d cysteine-bearing CovPDB proteins; %d are human, spanning %d complexes.\n",
  nrow(all_rows), nrow(human_rows), sum(human_rows$n_protein_ligand_complexes, na.rm = TRUE)
))
