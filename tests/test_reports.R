#!/usr/bin/env Rscript --vanilla
# ============================================================
# Test:     test_reports.R
# Purpose:  Verify report generation works end-to-end.
#           Run from project root: Rscript --vanilla tests/test_reports.R
# ============================================================

cat("=== CanProTarget Report Tests ===\n\n")

# Setup
project_root <- getwd()
source("R/app_helpers.R")
source("R/api_functions.R")
source("R/canprotarget_score.R")
source("R/pancancer_profile.R")
source("R/theme_canprotarget.R")
source("R/report_generator.R")

pass <- 0L
fail <- 0L
errors <- character(0)

test <- function(desc, expr) {
  result <- tryCatch(
    {
      val <- eval(expr)
      if (isTRUE(val)) {
        cat("  PASS:", desc, "\n")
        pass <<- pass + 1L
      } else {
        cat("  FAIL:", desc, "(returned", deparse(val), ")\n")
        fail <<- fail + 1L
        errors <<- c(errors, desc)
      }
    },
    error = function(e) {
      cat("  FAIL:", desc, "(error:", conditionMessage(e), ")\n")
      fail <<- fail + 1L
      errors <<- c(errors, paste0(desc, ": ", conditionMessage(e)))
    }
  )
  invisible(NULL)
}

output_dir <- file.path(tempdir(), "cpt_test_reports")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 1. Template discovery ----
cat("[Template Discovery]\n")

test("cpt_template_dir finds templates directory", {
  d <- cpt_template_dir(project_root)
  dir.exists(d)
})

test("cpt_available_reports lists gene_dependency", {
  "gene_dependency" %in% cpt_available_reports(project_root)
})

test("cpt_available_reports lists cysteine_target", {
  "cysteine_target" %in% cpt_available_reports(project_root)
})

test("cpt_render_report rejects unknown type", {
  tryCatch(
    { cpt_render_report("fake_type", project_root = project_root); FALSE },
    error = function(e) grepl("Unknown report type", conditionMessage(e))
  )
})

# ---- 2. Data helpers ----
cat("\n[Data Helpers]\n")

test("data_provenance_string works with NULL input", {
  result <- data_provenance_string(NULL)
  is.character(result) && nzchar(result)
})

test("data_provenance_string works with missing datasets", {
  fake <- list(datasets = list(list(name = "Test", version = "1.0", status = "missing")))
  result <- data_provenance_string(fake)
  is.character(result) && nzchar(result)
})

test("data_provenance_string formats available datasets", {
  fake <- list(datasets = list(
    list(name = "DepMap CRISPR", version = "23Q4", status = "ok"),
    list(name = "Bad", version = "x", status = "missing")
  ))
  result <- data_provenance_string(fake)
  grepl("DepMap CRISPR 23Q4", result) && !grepl("Bad", result)
})

# ---- 3. Cysteine report rendering ----
cat("\n[Cysteine Target Report]\n")

test("cysteine_target report renders for EGFR", {
  out <- cpt_report_cysteine_target(
    gene = "EGFR",
    project_root = project_root,
    output_dir = output_dir
  )
  file.exists(out) && file.size(out) > 1000
})

test("cysteine_target report contains branding", {
  files <- list.files(output_dir, pattern = "cysteine_target_EGFR", full.names = TRUE)
  if (!length(files)) return(FALSE)
  content <- readLines(files[1], warn = FALSE)
  text <- paste(content, collapse = " ")
  grepl("CanProTarget", text) && grepl("Bertoldo Lab", text) && grepl("github.com/Bertoldo-Lab", text)
})

test("cysteine_target report contains gene data", {
  files <- list.files(output_dir, pattern = "cysteine_target_EGFR", full.names = TRUE)
  if (!length(files)) return(FALSE)
  content <- paste(readLines(files[1], warn = FALSE), collapse = " ")
  grepl("EGFR", content) && grepl("Functional", content)
})

test("cysteine_target report with specific site_id", {
  out <- cpt_report_cysteine_target(
    gene = "EGFR",
    site_id = "EGFR_797",
    project_root = project_root,
    output_dir = output_dir
  )
  file.exists(out) && file.size(out) > 1000
})

# ---- 4. Gene dependency report rendering ----
cat("\n[Gene Dependency Report]\n")

test("gene_dependency report renders for KRAS/Pancreatic", {
  out <- cpt_report_gene_dependency(
    gene = "KRAS",
    subtype = "Pancreatic Adenocarcinoma",
    dataset = "CRISPR",
    project_root = project_root,
    output_dir = output_dir
  )
  file.exists(out) && file.size(out) > 1000
})

test("gene_dependency report contains branding", {
  files <- list.files(output_dir, pattern = "gene_dependency_KRAS", full.names = TRUE)
  if (!length(files)) return(FALSE)
  content <- paste(readLines(files[1], warn = FALSE), collapse = " ")
  grepl("CanProTarget", content) && grepl("Bertoldo Lab", content) && grepl("github.com/Bertoldo-Lab", content)
})

test("gene_dependency report contains analysis results", {
  files <- list.files(output_dir, pattern = "gene_dependency_KRAS", full.names = TRUE)
  if (!length(files)) return(FALSE)
  content <- paste(readLines(files[1], warn = FALSE), collapse = " ")
  grepl("KRAS", content) && grepl("Pancreatic", content) && grepl("CPT Score", content)
})

test("gene_dependency report contains pan-cancer plot", {
  files <- list.files(output_dir, pattern = "gene_dependency_KRAS", full.names = TRUE)
  if (!length(files)) return(FALSE)
  content <- paste(readLines(files[1], warn = FALSE), collapse = " ")
  # SVG or PNG image data present (base64 encoded images in self-contained HTML)
  grepl("data:image/png;base64", content)
})

# ---- 5. MCP API wrapper ----
cat("\n[MCP API Wrapper]\n")

test("api_generate_report validates report_type", {
  data_env <- list(data_dir = file.path(project_root, "data"))
  tryCatch(
    { api_generate_report("nonexistent", list(gene = "X"), data_env); FALSE },
    error = function(e) grepl("Unknown report_type", conditionMessage(e))
  )
})

test("api_generate_report validates required params for gene_dependency", {
  data_env <- list(data_dir = file.path(project_root, "data"))
  tryCatch(
    { api_generate_report("gene_dependency", list(), data_env); FALSE },
    error = function(e) grepl("gene", conditionMessage(e))
  )
})

test("api_generate_report produces report and returns metadata", {
  data_env <- list(data_dir = file.path(project_root, "data"))
  result <- api_generate_report(
    "cysteine_target",
    list(gene = "TP53"),
    data_env
  )
  result$status == "success" &&
    file.exists(result$file_path) &&
    !is.null(result$citation) &&
    !is.null(result$source) &&
    grepl("github.com/Bertoldo-Lab", result$source)
})

# ---- Summary ----
cat("\n=== Results ===\n")
cat(sprintf("  %d passed, %d failed\n", pass, fail))
if (fail > 0L) {
  cat("\nFailed tests:\n")
  for (e in errors) cat("  -", e, "\n")
  quit(status = 1L)
} else {
  cat("All tests passed.\n")
}

# Cleanup
unlink(output_dir, recursive = TRUE)
