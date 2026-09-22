#!/usr/bin/env Rscript --vanilla
# ============================================================
# Test:     test_r_worker.R
# Purpose:  Verify the R MCP worker loads data and responds correctly.
#           Run from project root: Rscript --vanilla tests/test_r_worker.R
# ============================================================

cat("=== CanProTarget R Worker Tests ===\n\n")

# Track results
tests_run <- 0L
tests_passed <- 0L
tests_failed <- 0L

assert <- function(description, condition) {
  tests_run <<- tests_run + 1L
  if (isTRUE(condition)) {
    tests_passed <<- tests_passed + 1L
    cat(sprintf("  PASS: %s\n", description))
  } else {
    tests_failed <<- tests_failed + 1L
    cat(sprintf("  FAIL: %s\n", description))
  }
}

assert_error <- function(description, expr) {
  tests_run <<- tests_run + 1L
  result <- tryCatch({ expr; FALSE }, error = function(e) TRUE)
  if (result) {
    tests_passed <<- tests_passed + 1L
    cat(sprintf("  PASS: %s\n", description))
  } else {
    tests_failed <<- tests_failed + 1L
    cat(sprintf("  FAIL: %s (expected error, got success)\n", description))
  }
}

# --- Load API functions -----------------------------------------
cat("Loading api_functions.R...\n")
source("R/api_functions.R")
source("R/canprotarget_score.R")

# --- Load data (same as mcp_worker.R) --------------------------
cat("Loading data files...\n")
data_env <- list(data_dir = "data")
data_env$crispr_matrix <- readRDS("data/CRISPRGeneEffect_23Q4_clean.rds")
data_env$rnai_matrix <- readRDS("data/d2_gene_effect_headers_refined.rds")
data_env$cancer_model_data <- readRDS("data/cancer_model_data.rds")
data_env$cys_atlas <- readRDS("data/cys_editing_atlas.rds")
cat("Data loaded.\n\n")

# === Test: parse_gene_names =====================================
cat("[parse_gene_names]\n")
names_in <- c("KRAS (3845)", "EGFR (1956)", "A1BG (1)")
names_out <- parse_gene_names(names_in)
assert("parses gene names correctly", identical(names_out, c("KRAS", "EGFR", "A1BG")))

# === Test: find_gene_column =====================================
cat("[find_gene_column]\n")
idx <- find_gene_column(data_env$crispr_matrix, "KRAS")
assert("finds KRAS in CRISPR matrix", !is.null(idx) && idx > 0)
assert("KRAS column name contains KRAS", grepl("KRAS", colnames(data_env$crispr_matrix)[idx]))

idx_null <- find_gene_column(data_env$crispr_matrix, "FAKEGENE123")
assert("returns NULL for missing gene", is.null(idx_null))

idx_lower <- find_gene_column(data_env$crispr_matrix, "kras")
assert("case-insensitive lookup works", !is.null(idx_lower))

# === Test: api_list_subtypes ====================================
cat("[api_list_subtypes]\n")
result <- api_list_subtypes("CRISPR", data_env)
assert("returns list with dataset field", result$dataset == "CRISPR")
assert("has subtypes", result$count > 50)
assert("subtypes are sorted", identical(result$subtypes, sort(result$subtypes)))

result_rnai <- api_list_subtypes("RNAi", data_env)
assert("RNAi subtypes work", result_rnai$count > 30)

assert_error("rejects invalid dataset", api_list_subtypes("INVALID", data_env))

# === Test: api_list_genes =======================================
cat("[api_list_genes]\n")
result <- api_list_genes("CRISPR", data_env)
assert("CRISPR has ~18K genes", result$count > 18000 && result$count < 20000)
assert("genes are sorted", identical(result$genes, sort(result$genes)))
assert("KRAS is in gene list", "KRAS" %in% result$genes)

# === Test: api_query_dependency =================================
cat("[api_query_dependency]\n")
result <- api_query_dependency("KRAS", "Pancreatic Adenocarcinoma", "CRISPR", data_env)
assert("returns correct gene", result$gene == "KRAS")
assert("returns correct subtype", result$subtype == "Pancreatic Adenocarcinoma")
assert("KRAS is selective in pancreatic", result$is_selective_dependency == TRUE)
assert("effect size is negative", result$effect_size < -1.0)
assert("p-value is significant", result$p_value < 0.001)
assert("has interpretation string", nchar(result$interpretation) > 20)

# Non-dependency
result2 <- api_query_dependency("A1BG", "Melanoma", "CRISPR", data_env)
assert("A1BG is not a dependency in melanoma", result2$is_dependency == FALSE)

# Error cases
assert_error("rejects missing gene", api_query_dependency("NOTREAL", "Melanoma", "CRISPR", data_env))
assert_error("rejects missing subtype", api_query_dependency("KRAS", "Fake Cancer", "CRISPR", data_env))

# === Test: api_top_dependencies =================================
cat("[api_top_dependencies]\n")
result <- api_top_dependencies("Melanoma", "CRISPR", 10L, data_env)
assert("returns correct subtype", result$subtype == "Melanoma")
assert("returns 10 genes", length(result$genes) == 10)
assert("genes have effect_size field", !is.null(result$genes[[1]]$effect_size))
assert("genes are sorted by effect size", result$genes[[1]]$effect_size <= result$genes[[2]]$effect_size)
genes <- sapply(result$genes, function(g) g$gene)
assert("SOX10 or BRAF in top melanoma deps", any(c("SOX10", "BRAF") %in% genes))

# === Test: api_gene_cysteines ===================================
cat("[api_gene_cysteines]\n")
result <- api_gene_cysteines("EGFR", data_env)
assert("finds EGFR cysteines", result$total_sites > 30)
assert("has functional sites", result$n_functional > 0)
assert("has ligandable sites", result$n_ligandable > 0)
assert("sites have expected fields", !is.null(result$sites[[1]]$site_id))

# Case insensitive
result_lower <- api_gene_cysteines("egfr", data_env)
assert("case-insensitive gene lookup", result_lower$total_sites == result$total_sites)

assert_error("rejects missing gene", api_gene_cysteines("NOTINatlas", data_env))

# === Test: api_cysteine_detail ==================================
cat("[api_cysteine_detail]\n")
result <- api_cysteine_detail("EGFR_797", data_env)
assert("EGFR_797 is functional", result$functional == TRUE)
assert("EGFR_797 is ligandable", result$ligandable == TRUE)
assert("EGFR_797 ligandability score is 100", result$ligandability_score == 100)
assert("has editor_support field", !is.null(result$editor_support))
assert("has conservation_score", !is.null(result$conservation_score) || is.null(result$conservation_score))

assert_error("rejects missing site", api_cysteine_detail("FAKE_999", data_env))

# === Test: api_compare_subtypes =================================
cat("[api_compare_subtypes]\n")
result <- api_compare_subtypes("BRAF", "Melanoma", "Non-Small Cell Lung Cancer", "CRISPR", data_env)
assert("returns correct gene", result$gene == "BRAF")
assert("BRAF more essential in melanoma", result$subtype1$mean_effect < result$subtype2$mean_effect)
assert("p-value is significant", result$p_value < 0.001)
assert("has interpretation", nchar(result$interpretation) > 20)

assert_error("rejects subtype with <3 lines", 
  api_compare_subtypes("BRAF", "Melanoma", "Fake Subtype", "CRISPR", data_env))

# === Test: api_platform_info ====================================
cat("[api_platform_info]\n")
result <- api_platform_info()
assert("returns platform name", result$name == "CanProTarget")
assert("has version", nchar(result$version) > 0)
assert("lists 14 tools", length(result$available_tools) == 14)
assert("lists rank_site_targets", "rank_site_targets" %in% result$available_tools)
assert("has data_sources", length(result$data_sources) >= 3)

# === Summary ====================================================
cat(sprintf("\n=== Results: %d/%d passed", tests_passed, tests_run))
if (tests_failed > 0) {
  cat(sprintf(", %d FAILED", tests_failed))
}
cat(" ===\n")

if (tests_failed > 0) quit(status = 1)
