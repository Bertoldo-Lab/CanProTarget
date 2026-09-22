#!/usr/bin/env Rscript
# ============================================================
# Tests for the ADME / developability CPT dimension.
#
# The contract being protected:
#   1. ADME is COMPUTED for genes with covalent probe coverage.
#   2. Its default weight is 0, so composite scores are IDENTICAL to the
#      pre-implementation baseline. This is the guarantee that lets the
#      dimension ship before the formula is signed off.
#   3. Raising the weight does change scores (the slider is real).
#   4. The score is decoupled from probe count, i.e. it measures developability
#      and not how well-screened a target happens to be.
#
# Run from project root:  Rscript --vanilla tests/test_adme_dimension.R
# ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(magrittr); library(tibble)
})

pass <- 0L
fail <- 0L
failed <- character(0)
assert <- function(label, cond) {
  if (isTRUE(cond)) {
    pass <<- pass + 1L
    cat("  PASS:", label, "\n")
  } else {
    fail <<- fail + 1L
    failed <<- c(failed, label)
    cat("  FAIL:", label, "\n")
  }
}

source("R/functions.R", local = TRUE)
source("R/canprotarget_score.R", local = TRUE)

# ---- 1. defaults ------------------------------------------------
cat("\n[defaults]\n")
assert("ADME default weight is 0 (composite must not move)",
       identical(as.numeric(CPT_DEFAULT_WEIGHTS$adme_druggability), 0))
assert("proposed weight constant is documented",
       exists("CPT_ADME_PROPOSED_WEIGHT") &&
         is.numeric(CPT_ADME_PROPOSED_WEIGHT))

# ---- 2. window helper -------------------------------------------
cat("\n[cpt_window_score]\n")
assert("inside the ideal window scores 100",
       identical(as.numeric(cpt_window_score(250, 0, 300, -Inf, 500)), 100))
assert("at the zero point scores 0",
       identical(as.numeric(cpt_window_score(500, 0, 300, -Inf, 500)), 0))
assert("beyond the zero point clamps at 0",
       identical(as.numeric(cpt_window_score(900, 0, 300, -Inf, 500)), 0))
assert("midway down the shoulder is between 0 and 100", {
  v <- as.numeric(cpt_window_score(400, 0, 300, -Inf, 500))
  v > 40 && v < 60
})
assert("low shoulder penalises too-polar logP", {
  v <- as.numeric(cpt_window_score(-1.5, 0, 3, -3, 6))
  v > 0 && v < 100
})
assert("NA in, NA out", is.na(cpt_window_score(NA, 0, 300, -Inf, 500)))

# ---- 3. graceful degradation ------------------------------------
cat("\n[score_adme_druggability fallbacks]\n")
assert("NULL lookup yields all NA",
       all(is.na(score_adme_druggability(c("EGFR", "KRAS"), NULL))))
assert("wrong-shaped frame yields all NA",
       all(is.na(score_adme_druggability("EGFR", data.frame(x = 1)))))
assert("returns one value per gene asked",
       length(score_adme_druggability(c("A", "B", "C"), NULL)) == 3L)

# ---- 4. live data ----------------------------------------------
bind_p <- file.path("data", "protein_binding_lookup_preprocessed.rds")
adme_p <- file.path("data", "swissadme_preprocessed.rds")
atlas_p <- file.path("data", "cys_editing_atlas.rds")

if (!file.exists(bind_p) || !file.exists(adme_p)) {
  cat("\n[live data] SKIP - chemoproteomics RDS not present\n")
} else {
  cat("\n[live lookup build]\n")
  bind <- readRDS(bind_p)
  adme <- readRDS(adme_p)
  lk <- cpt_build_adme_gene_scores(bind, adme)

  assert("lookup is non-empty", nrow(lk) > 1000)
  assert("has the expected columns",
         all(c("gene_key", "adme_score", "n_probes", "best_probe") %in% names(lk)))
  assert("gene keys are upper case", all(lk$gene_key == toupper(lk$gene_key)))
  assert("scores are within 0-100",
         all(lk$adme_score >= 0 & lk$adme_score <= 100, na.rm = TRUE))
  assert("no duplicate genes", !any(duplicated(lk$gene_key)))

  # The bug this guards: taking max() across probes made the score a proxy for
  # how many probes hit the gene (Spearman +0.52). The mean decouples it.
  rho <- suppressWarnings(cor(lk$adme_score, lk$n_probes,
                              method = "spearman", use = "complete.obs"))
  assert("score is decoupled from probe count (|rho| < 0.2)", abs(rho) < 0.2)
  cat("       (spearman vs n_probes =", round(rho, 3), ")\n")

  assert("score actually discriminates (>100 distinct values)",
         length(unique(round(lk$adme_score, 1))) > 100)

  cat("\n[gene-level expectations]\n")
  assert("EGFR has a score (probe-covered)",
         !is.na(score_adme_druggability("EGFR", lk)))
  assert("KRAS has no score (no probe at CR >= 4)",
         is.na(score_adme_druggability("KRAS", lk)))
  assert("Entrez-suffixed names still resolve",
         !is.na(score_adme_druggability("EGFR (1956)", lk)))
  assert("lower case resolves", !is.na(score_adme_druggability("egfr", lk)))

  if (file.exists(atlas_p)) {
    cat("\n[composite invariance at weight 0]\n")
    atlas <- readRDS(atlas_p)
    stats <- data.frame(
      gene_name   = c("KRAS", "EGFR", "BRAF", "VCP", "WEE1", "MYC"),
      cancer_mean = c(-2.0312, -0.2723, -1.0150, -2.6028, -2.4001, -1.1000),
      effect_size = c(-1.3865,  0.0753, -0.9202, -0.1067,  0.0175, -0.3000),
      stringsAsFactors = FALSE
    )
    without <- cpt_score(stats, cys_atlas = atlas, adme_data = NULL)
    with_it <- cpt_score(stats, cys_atlas = atlas, adme_data = lk)

    assert("composite identical with vs without ADME data",
           isTRUE(all.equal(without$cpt_score, with_it$cpt_score)))
    assert("ADME values are still reported",
           any(!is.na(with_it$dim_adme)))
    assert("dimension count rises where ADME exists",
           any(with_it$n_dimensions > without$n_dimensions))

    cat("\n[weight actually does something]\n")
    w <- CPT_DEFAULT_WEIGHTS
    w$adme_druggability <- CPT_ADME_PROPOSED_WEIGHT
    raised <- cpt_score(stats, cys_atlas = atlas, adme_data = lk, weights = w)
    assert("raising the weight changes scores",
           !isTRUE(all.equal(with_it$cpt_score, raised$cpt_score)))
    assert("genes without ADME are unaffected by the weight", {
      i <- which(is.na(with_it$dim_adme))
      length(i) == 0 || isTRUE(all.equal(with_it$cpt_score[i], raised$cpt_score[i]))
    })
  }
}

cat(sprintf("\n=== Results: %d/%d passed ===\n", pass, pass + fail))
if (fail > 0) {
  cat("Failed tests:\n")
  for (f in failed) cat("  -", f, "\n")
  quit(status = 1)
}
