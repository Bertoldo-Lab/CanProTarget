# Offline-only helpers to build swissadme_preprocessed.rds from raw exports.
# Sourced by docs/scripts/preprocess_data.R (not loaded by the Shiny app).

read_swissadme_from_raw_files <- function(data_dir) {
  csv_path <- file.path(data_dir, "swissadme.csv")
  if (!file.exists(csv_path)) {
    stop("swissadme.csv not found in ", data_dir, call. = FALSE)
  }

  adme <- utils::read.csv(csv_path, stringsAsFactors = FALSE, check.names = TRUE)

  link_rds <- file.path(data_dir, "backus_swiss_link.rds")
  link_xlsx <- file.path(data_dir, "backus_swiss_link.xlsx")
  link <- if (file.exists(link_rds)) {
    readRDS(link_rds)
  } else if (file.exists(link_xlsx)) {
    readxl::read_excel(link_xlsx)
  } else {
    NULL
  }

  if (!is.null(link)) {
    link_long <- link %>%
      tidyr::separate_rows(updated_compound_names, sep = ",\\s*") %>%
      dplyr::rename(probe_name = updated_compound_names)

    adme <- adme %>%
      dplyr::left_join(link_long, by = "Molecule")
  }

  if (is.null(link)) {
    if ("Molecule" %in% colnames(adme)) {
      adme <- dplyr::rename(adme, probe_name = Molecule)
    }
  }

  if (!"probe_name" %in% colnames(adme)) {
    adme$probe_name <- NA_character_
  }

  adme
}
