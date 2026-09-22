#!/usr/bin/env Rscript
# Build data/backus_swiss_link.rds from data/backus_swiss_link.xlsx.
# Run from project root:  Rscript docs/scripts/backus_swiss_link_xlsx_to_rds.R

args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[[1]] else "data"

xlsx <- file.path(data_dir, "backus_swiss_link.xlsx")
rds <- file.path(data_dir, "backus_swiss_link.rds")

if (!file.exists(xlsx)) {
  stop("Missing ", xlsx, call. = FALSE)
}

if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("Install readxl: install.packages(\"readxl\")", call. = FALSE)
}

cat("Reading ", xlsx, " ...\n", sep = "")
tbl <- readxl::read_excel(xlsx)
saveRDS(tbl, rds, compress = "xz")
cat("Wrote ", rds, " (", nrow(tbl), " rows x ", ncol(tbl), " cols)\n", sep = "")
