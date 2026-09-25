#!/usr/bin/env Rscript --vanilla
# ============================================================
# Test:     test_cpt_agent_workflows.R
# Purpose:  CPT scoring correctness (background ranks), assess_target,
#           rank_targets, and light validation fixtures.
#           Run from project root:
#             Rscript --vanilla tests/test_cpt_agent_workflows.R
#
# Note:     Fixture expectations cover score computation and DepMap-style
#           dependency patterns. They are regression anchors, not a
#           biological validation of the weights.
# ============================================================

cat("=== CanProTarget CPT / agent workflow tests ===\n\n")

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

cat("Loading sources and data...\n")
source("R/api_functions.R")
source("R/canprotarget_score.R")
source("R/pancancer_profile.R")

data_env <- list(data_dir = "data")
data_env$crispr_matrix <- readRDS("data/CRISPRGeneEffect_23Q4_clean.rds")
data_env$rnai_matrix <- readRDS("data/d2_gene_effect_headers_refined.rds")
data_env$cancer_model_data <- readRDS("data/cancer_model_data.rds")
data_env$cys_atlas <- readRDS("data/cys_editing_atlas.rds")
cat("Data loaded.\n\n")

# --- Single-gene score must not collapse dep/sel to 50 ----------
cat("[cpt_score_single background ranking]\n")
kras_pdac <- api_canprotarget_score(
  "KRAS", "Pancreatic Adenocarcinoma", "CRISPR", data_env
)
assert("returns gene KRAS", identical(kras_pdac$gene, "KRAS"))
assert("returns a numeric CPT score", is.numeric(kras_pdac$cpt_score) && !is.na(kras_pdac$cpt_score))
assert("dependency dimension is active",
       "dependency_strength" %in% kras_pdac$dimensions_active)
assert("selectivity dimension is active",
       "cancer_selectivity" %in% kras_pdac$dimensions_active)
# Regression: old single-row percentile always returned 50
dep_sc <- kras_pdac$dimensions$dependency_strength$score
sel_sc <- kras_pdac$dimensions$cancer_selectivity$score
assert("KRAS dependency percentile is not the dummy median 50",
       !is.null(dep_sc) && !is.na(dep_sc) && abs(dep_sc - 50) > 5)
assert("KRAS selectivity percentile is not the dummy median 50",
       !is.null(sel_sc) && !is.na(sel_sc) && abs(sel_sc - 50) > 5)
# Strong PDAC selective dep should sit high on both axes
assert("KRAS ranks high on dependency in PDAC", dep_sc > 80)
assert("KRAS ranks high on selectivity in PDAC", sel_sc > 80)
assert("reports active/missing dimension lists",
       length(kras_pdac$dimensions_active) >= 2 &&
         length(kras_pdac$dimensions_missing) >= 1)
assert("ADME is listed as missing (expected until data lands)",
       "adme_druggability" %in% kras_pdac$dimensions_missing)
assert("includes caveats", length(kras_pdac$caveats) >= 1)
assert("n_dimensions_used matches active count",
       kras_pdac$n_dimensions_used == length(kras_pdac$dimensions_active))

# Non-essential gene: raw mean near zero; percentile near mid-pack is normal
# (most genes sit around 0 effect). Must still score well below KRAS/PDAC.
a1bg <- api_canprotarget_score("A1BG", "Melanoma", "CRISPR", data_env)
a1bg_dep <- a1bg$dimensions$dependency_strength$score
a1bg_raw <- a1bg$dimensions$dependency_strength$raw_value
assert("A1BG is not essential by raw mean in melanoma",
       !is.null(a1bg_raw) && a1bg_raw > -0.5)
assert("A1BG dependency percentile well below KRAS/PDAC",
       !is.null(a1bg_dep) && a1bg_dep < dep_sc - 20)
assert("A1BG CPT score below KRAS/PDAC",
       is.numeric(a1bg$cpt_score) && a1bg$cpt_score < kras_pdac$cpt_score)

# --- rank_targets -----------------------------------------------
cat("[api_rank_targets]\n")
ranked <- api_rank_targets(
  subtype = "Pancreatic Adenocarcinoma",
  dataset = "CRISPR",
  n = 15L,
  pool = 80L,
  require_ligandable = FALSE,
  data_env = data_env
)
assert("returns requested subtype", ranked$subtype == "Pancreatic Adenocarcinoma")
assert("returns up to n genes", length(ranked$genes) > 0 && length(ranked$genes) <= 15)
assert("genes sorted by cpt_score descending", {
  scores <- sapply(ranked$genes, function(g) g$cpt_score)
  all(diff(scores) <= 0)
})
ranked_genes <- sapply(ranked$genes, function(g) g$gene)
assert("KRAS appears in top PDAC CPT ranks (pool of selective deps)",
       "KRAS" %in% ranked_genes)

# --- assess_target ----------------------------------------------
cat("[api_assess_target]\n")
assessed <- api_assess_target(
  "KRAS", "Pancreatic Adenocarcinoma", "CRISPR", data_env
)
assert("summary is non-empty", nchar(assessed$summary) > 40)
assert("includes dependency block", isTRUE(assessed$dependency$is_selective_dependency))
assert("includes cpt_score block", !is.na(assessed$cpt_score$cpt_score))
assert("includes next_actions", length(assessed$next_actions) >= 2)
assert("next_actions name real tools", {
  tools <- sapply(assessed$next_actions, function(a) a$tool)
  all(tools %in% c(
    "cysteine_detail", "gene_cysteines", "pancancer_profile",
    "generate_report", "canprotarget_score", "compare_subtypes"
  ))
})
assert("has citation string", grepl("CanProTarget", assessed$citation))
assert("assess includes external_resources", !is.null(assessed$external_resources$depmap))
assert("assess includes provenance", !is.null(assessed$provenance$version))

egfr_nsclc <- api_assess_target(
  "EGFR", "Non-Small Cell Lung Cancer", "CRISPR", data_env
)
assert("EGFR assess returns gene", identical(egfr_nsclc$gene, "EGFR"))
assert("EGFR has cysteine atlas coverage",
       isTRUE(egfr_nsclc$cysteines$available) &&
         egfr_nsclc$cysteines$total_sites > 0)

# --- platform_info tool list ------------------------------------
cat("[api_platform_info tools]\n")
info <- api_platform_info()
assert("lists assess_target", "assess_target" %in% info$available_tools)
assert("lists rank_targets", "rank_targets" %in% info$available_tools)
assert("lists rank_site_targets", "rank_site_targets" %in% info$available_tools)
assert("tool count is 14", length(info$available_tools) == 14L)
assert("credits authors present", grepl("Ong JP", info$credits$authors))

# Cysteine ordering + links
cat("[cysteine quality]\n")
egfr_cys <- api_gene_cysteines("EGFR", data_env)
assert("gene_cysteines has external_resources", !is.null(egfr_cys$external_resources$alphafold) ||
         !is.null(egfr_cys$external_resources$uniprot_search))
assert("first site is functional+ligandable when any exist", {
  if (egfr_cys$n_functional_ligandable > 0) {
    isTRUE(egfr_cys$sites[[1]]$functional_ligandable)
  } else {
    TRUE
  }
})
detail <- api_cysteine_detail("EGFR_797", data_env)
assert("cysteine_detail has alphafold or uniprot link",
       !is.null(detail$external_resources$alphafold) ||
         !is.null(detail$external_resources$uniprot))

# --- Known-target recovery fixtures (score computation + DepMap patterns)
# These check that well-established targets are still recovered in their
# expected context. They are not a validation of the CPT weights themselves.
cat("[validation fixtures]\n")

fixtures <- list(
  list(
    id = "KRAS_PDAC",
    gene = "KRAS",
    subtype = "Pancreatic Adenocarcinoma",
    expect_selective_dependency = TRUE,
    expect_min_dep_percentile = 80,
    expect_min_sel_percentile = 80,
    note = "Well-known selective dependency in PDAC CRISPR screens"
  ),
  list(
    id = "BRAF_Melanoma",
    gene = "BRAF",
    subtype = "Melanoma",
    expect_selective_dependency = TRUE,
    expect_min_dep_percentile = 70,
    expect_min_sel_percentile = 70,
    note = "Well-known melanoma dependency; used elsewhere in test suite"
  ),
  list(
    id = "EGFR_C797_exists",
    gene = "EGFR",
    site_id = "EGFR_797",
    expect_site_functional = TRUE,
    expect_site_ligandable = TRUE,
    note = "Atlas site used as cysteine report example; not a CPT score claim"
  )
)

for (fx in fixtures) {
  if (!is.null(fx$site_id)) {
    site <- api_cysteine_detail(fx$site_id, data_env)
    assert(paste0(fx$id, ": site found"), identical(site$site_id, fx$site_id) ||
             grepl(fx$gene, site$gene_symbol, ignore.case = TRUE))
    if (isTRUE(fx$expect_site_functional)) {
      assert(paste0(fx$id, ": functional flag"), isTRUE(site$functional))
    }
    if (isTRUE(fx$expect_site_ligandable)) {
      assert(paste0(fx$id, ": ligandable flag"), isTRUE(site$ligandable))
    }
  } else {
    dep <- api_query_dependency(fx$gene, fx$subtype, "CRISPR", data_env)
    sc <- api_canprotarget_score(fx$gene, fx$subtype, "CRISPR", data_env)
    if (isTRUE(fx$expect_selective_dependency)) {
      assert(paste0(fx$id, ": selective dependency"),
             isTRUE(dep$is_selective_dependency))
    }
    if (!is.null(fx$expect_min_dep_percentile)) {
      assert(paste0(fx$id, ": dependency percentile threshold"),
             sc$dimensions$dependency_strength$score >= fx$expect_min_dep_percentile)
    }
    if (!is.null(fx$expect_min_sel_percentile)) {
      assert(paste0(fx$id, ": selectivity percentile threshold"),
             sc$dimensions$cancer_selectivity$score >= fx$expect_min_sel_percentile)
    }
  }
}

# --- Summary ----------------------------------------------------
cat(sprintf("\n=== Results: %d/%d passed", tests_passed, tests_run))
if (tests_failed > 0) {
  cat(sprintf(", %d FAILED", tests_failed))
}
cat(" ===\n")
if (tests_failed > 0) quit(status = 1)
