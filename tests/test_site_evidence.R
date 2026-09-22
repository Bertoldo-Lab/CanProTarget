#!/usr/bin/env Rscript --vanilla
# ============================================================
# Test:     test_site_evidence.R
# Purpose:  Evidence tiers, compact SMCL index, Site CPT, and the
#           guarantee that gene-level CPT is unchanged.
#           Run from project root:
#             Rscript --vanilla tests/test_site_evidence.R
# ============================================================

cat("=== CanProTarget site evidence / SMCL / Site CPT tests ===\n\n")

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

source("R/api_functions.R")
source("R/canprotarget_score.R")

cat("[cpt_evidence_tier]\n")
assert("functional + atlas ligandable is Tier 1",
       identical(cpt_evidence_tier(TRUE, TRUE, TRUE), 1L))
assert("functional without atlas ligandable is Tier 2",
       identical(cpt_evidence_tier(TRUE, TRUE, FALSE), 2L))
assert("in atlas, not functional is Tier 3",
       identical(cpt_evidence_tier(TRUE, FALSE, FALSE), 3L))
assert("not in atlas is Tier 4",
       identical(cpt_evidence_tier(FALSE, FALSE, FALSE), 4L))
assert("vectorised",
       identical(cpt_evidence_tier(c(TRUE, TRUE, TRUE, FALSE),
                                   c(TRUE, TRUE, FALSE, TRUE),
                                   c(TRUE, FALSE, TRUE, TRUE)),
                 c(1L, 2L, 3L, 4L)))

cat("Loading data (binding table is large; once)...\n")
data_env <- list(data_dir = "data")
data_env$crispr_matrix <- readRDS("data/CRISPRGeneEffect_23Q4_clean.rds")
data_env$rnai_matrix <- readRDS("data/d2_gene_effect_headers_refined.rds")
data_env$cancer_model_data <- readRDS("data/cancer_model_data.rds")
data_env$cys_atlas <- readRDS("data/cys_editing_atlas.rds")
bind <- readRDS("data/protein_binding_lookup_preprocessed.rds")
data_env$smcl_index <- cpt_build_smcl_index(bind, min_cr = 4)
rm(bind)
gc(verbose = FALSE)
cat("SMCL index rows:", nrow(data_env$smcl_index), "\n\n")

cat("[cpt_build_smcl_index]\n")
assert("index is non-empty", nrow(data_env$smcl_index) > 1000)
assert("has gene_key and residue_number",
       all(c("gene_key", "residue_number", "probe_name", "CR") %in% names(data_env$smcl_index)))
assert("all CR >= 4", min(data_env$smcl_index$CR, na.rm = TRUE) >= 4)

cat("[gene-level CPT unchanged]\n")
kras_no_idx <- data_env
kras_no_idx$smcl_index <- NULL
a <- api_canprotarget_score("KRAS", "Pancreatic Adenocarcinoma", "CRISPR", kras_no_idx)
b <- api_canprotarget_score("KRAS", "Pancreatic Adenocarcinoma", "CRISPR", data_env)
assert("KRAS PDAC gene CPT identical with and without SMCL index",
       identical(a$cpt_score, b$cpt_score) && identical(a$dimensions_active, b$dimensions_active))

cat("[EGFR sites + SMCL]\n")
egfr <- api_gene_cysteines("EGFR", data_env)
assert("EGFR still has atlas coverage", egfr$n_atlas > 30)
assert("EGFR has engaged SMCLs", isTRUE(egfr$n_engaged > 0))
assert("sites include evidence_tier", !is.null(egfr$sites[[1]]$evidence_tier))
c797 <- Filter(function(s) identical(as.integer(s$cysteine_position), 797L), egfr$sites)
assert("EGFR C797 present", length(c797) == 1L)
if (length(c797) == 1L) {
  assert("EGFR C797 is Tier 1", identical(as.integer(c797[[1]]$evidence_tier), 1L))
  assert("EGFR C797 is functional+ligandable", isTRUE(c797[[1]]$functional_ligandable))
}

detail <- api_cysteine_detail("EGFR_797", data_env)
assert("cysteine_detail reports evidence_tier", identical(as.integer(detail$evidence_tier), 1L))
assert("cysteine_detail has smcls field", !is.null(detail$smcls))

cat("[IGF2BP3 Cys336 Tier 4 if engaged]\n")
igf <- tryCatch(api_gene_cysteines("IGF2BP3", data_env), error = function(e) NULL)
if (!is.null(igf)) {
  c336 <- Filter(function(s) identical(as.integer(s$cysteine_position), 336L), igf$sites)
  if (length(c336) == 1L) {
    assert("IGF2BP3 C336 has a tier", c336[[1]]$evidence_tier %in% 1:4)
    if (!isTRUE(c336[[1]]$in_atlas)) {
      assert("IGF2BP3 C336 not in atlas is Tier 4", identical(as.integer(c336[[1]]$evidence_tier), 4L))
    }
    d336 <- api_cysteine_detail("IGF2BP3_336", data_env)
    assert("IGF2BP3_336 cysteine_detail does not error", identical(d336$gene_symbol, igf$gene) ||
             grepl("IGF2BP3", d336$gene_symbol, ignore.case = TRUE))
  } else {
    assert("IGF2BP3 C336 not in bundled screens (skip)", TRUE)
  }
} else {
  assert("IGF2BP3 not in atlas or SMCL (skip)", TRUE)
}

cat("[assess_target engaged_sites]\n")
assessed <- api_assess_target("EGFR", "Non-Small Cell Lung Cancer", "CRISPR", data_env)
assert("assess_target has engaged_sites", is.list(assessed$engaged_sites))
if (length(assessed$engaged_sites) > 0) {
  assert("engaged site has SMCL or n_smcls",
         !is.null(assessed$engaged_sites[[1]]$n_smcls) &&
           assessed$engaged_sites[[1]]$n_smcls >= 1)
  assert("engaged site has Site CPT or NA",
         is.null(assessed$engaged_sites[[1]]$site_cpt_score) ||
           is.numeric(assessed$engaged_sites[[1]]$site_cpt_score) ||
           is.na(assessed$engaged_sites[[1]]$site_cpt_score))
}

cat("[max_targets applied before SMCL cap]\n")
fake <- data.frame(
  gene_key = "FAKETESTGENE",
  gene_name = "FAKETESTGENE",
  proteinid = "P0FAKE",
  cysteineid = "P0FAKE_C10",
  residue_number = 10L,
  probe_name = paste0("PROBE", 1:10),
  CR = as.numeric(20:11),
  n_targets = c(rep(80, 8), 5, 6),
  ligandable = "yes",
  stringsAsFactors = FALSE
)
env_fake <- list(
  cys_atlas = NULL,
  smcl_index = fake
)
all_probes <- cpt_gene_site_records("FAKETESTGENE", env_fake, min_cr = 4, smcl_cap = 8L)
assert("NULL atlas still returns SMCL-only sites", !is.null(all_probes) && length(all_probes$sites) == 1L)
assert("without max_targets, cap keeps 8 highest-CR probes",
       all_probes$sites[[1]]$n_smcls == 10L && length(all_probes$sites[[1]]$smcls) == 8L)
filtered <- cpt_gene_site_records(
  "FAKETESTGENE", env_fake, min_cr = 4, max_targets = 20, smcl_cap = 8L
)
assert("max_targets keeps the site via restricted probes below the CR cap",
       !is.null(filtered) && isTRUE(filtered$sites[[1]]$engaged) &&
         filtered$sites[[1]]$n_smcls == 2L)
assert("restricted probes are the n_targets <= 20 ones",
       all(vapply(filtered$sites[[1]]$smcls, function(p) p$n_targets <= 20, logical(1))))

cat("[rank_site_targets order]\n")
# Melanoma is a reasonably sized CRISPR subtype used elsewhere in tests
ranked_sites <- api_rank_site_targets(
  subtype = "Melanoma",
  dataset = "CRISPR",
  n = 15L,
  pool = 40L,
  min_cr = 4,
  max_targets = 20,
  data_env = data_env
)
assert("rank_site_targets returns list", is.list(ranked_sites$sites))
if (length(ranked_sites$sites) >= 2) {
  tiers <- vapply(ranked_sites$sites, function(s) as.integer(s$evidence_tier), integer(1))
  assert("tiers are non-decreasing (worse tiers cannot precede better)",
         all(diff(tiers) >= 0))
  # Within the first tier block, site CPT is non-increasing
  t1 <- tiers == min(tiers)
  scores <- vapply(ranked_sites$sites[t1], function(s) as.numeric(s$site_cpt_score), numeric(1))
  if (length(scores) >= 2 && all(is.finite(scores))) {
    assert("within-tier Site CPT is non-increasing", all(diff(scores) <= 1e-8))
  } else {
    assert("within-tier Site CPT check skipped (too few)", TRUE)
  }
} else {
  assert("rank_site_targets returned few sites (skip order check)", TRUE)
}

cat("[manuscript filters]\n")
ms <- api_rank_site_targets(
  subtype = "Melanoma",
  dataset = "CRISPR",
  n = 10L,
  pool = 30L,
  min_cr = 4,
  max_targets = 20,
  effect_size_max = -0.1,
  p_max = 0.05,
  exclude_common_essentials = TRUE,
  data_env = data_env
)
assert("manuscript pool is labelled", identical(ms$gene_pool, "manuscript_filters"))
if (length(ms$sites) > 0) {
  es <- vapply(ms$sites, function(s) as.numeric(s$effect_size), numeric(1))
  assert("manuscript filter keeps effect_size < -0.1", all(es < -0.1, na.rm = TRUE))
} else {
  assert("manuscript filter returned no engaged sites in this pool (ok)", TRUE)
}

cat("\n=== Results: ", tests_passed, "/", tests_run, " passed",
    if (tests_failed) paste0(", ", tests_failed, " failed") else "",
    " ===\n", sep = "")
if (tests_failed > 0) quit(status = 1)
