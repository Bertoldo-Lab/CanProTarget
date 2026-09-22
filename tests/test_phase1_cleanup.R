#!/usr/bin/env Rscript --vanilla
# ============================================================
# Test:     test_phase1_cleanup.R
# Purpose:  Phase 1.1-1.4 cleanup regressions — gene suggestions,
#           group-comparison data isolation, BOILED-Egg axes,
#           missing-data UI helpers, shared_data key convention,
#           and static wiring checks in modules.
# Run from project root:
#   Rscript --vanilla tests/test_phase1_cleanup.R
# ============================================================

cat("=== Phase 1 cleanup tests ===\n\n")

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

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

source("R/app_helpers.R", local = TRUE)
source("R/functions.R", local = TRUE)

# ---------------------------------------------------------------------------
cat("[cpt_suggest_symbols]\n")
# ---------------------------------------------------------------------------
cand <- c("KRAS", "HRAS", "NRAS", "EGFR", "BRAF", "TP53", "CDK2")

sug <- cpt_suggest_symbols(c("KRAS", "kras", "EGF", "KRASX", "ZZZNOPE"), cand, n = 5)
assert("exact match returns that gene", identical(sug[["KRAS"]], "KRAS"))
assert("case-insensitive exact match", identical(sug[["kras"]], "KRAS"))
assert("prefix EGF suggests EGFR", "EGFR" %in% sug[["EGF"]])
assert("typo KRASX still suggests something RAS-like",
       length(sug[["KRASX"]]) >= 1 && any(grepl("RAS", sug[["KRASX"]], ignore.case = TRUE)))
assert("unknown gene can return empty suggestions",
       is.character(sug[["ZZZNOPE"]]))

empty_q <- cpt_suggest_symbols(character(0), cand)
assert("empty query returns empty list", length(empty_q) == 0)

empty_c <- cpt_suggest_symbols("KRAS", character(0))
assert("empty candidates returns empty for query",
       length(empty_c[["KRAS"]]) == 0)

capped <- cpt_suggest_symbols("R", c("RA", "RB", "RC", "RD", "RE", "RF", "RG"), n = 3)
assert("respects n cap", length(capped[["R"]]) <= 3)

# ---------------------------------------------------------------------------
cat("\n[cpt_data_unavailable_card]\n")
# ---------------------------------------------------------------------------
card <- cpt_data_unavailable_card("Missing probes", "Need the RDS file.", "CPT-1007")
assert("card is a shiny tag", inherits(card, "shiny.tag"))
assert("card class is alert-warning",
       grepl("alert-warning", paste(card$attribs$class, collapse = " ")))
card_html <- as.character(card)
assert("card embeds CPT code", grepl("CPT-1007", card_html, fixed = TRUE))
assert("card mentions DATA_PROVENANCE.md", grepl("DATA_PROVENANCE.md", card_html, fixed = TRUE))
assert("card embeds title", grepl("Missing probes", card_html, fixed = TRUE))

card_no_code <- cpt_data_unavailable_card("T", "D", NULL)
assert("card works without code", inherits(card_no_code, "shiny.tag"))

# ---------------------------------------------------------------------------
cat("\n[sanitize_subtype]\n")
# ---------------------------------------------------------------------------
assert("sanitize_subtype alias exists", exists("sanitize_subtype", mode = "function"))
assert("sanitize strips punctuation",
       identical(sanitize_subtype("Non-Small Cell Lung Cancer"),
                 "Non_Small_Cell_Lung_Cancer"))
assert("cpt_sanitize_subtype_key same as alias",
       identical(sanitize_subtype("Melanoma"), cpt_sanitize_subtype_key("Melanoma")))

# ---------------------------------------------------------------------------
cat("\n[cpt_group_comparison_df]\n")
# ---------------------------------------------------------------------------
# Tiny synthetic matrix: 6 lines x 2 genes
mat <- matrix(
  c(
    -1.2, -0.1,   # CL1 cancer of interest, KRAS dep
    -1.0,  0.0,   # CL2 cancer of interest
     0.1, -0.2,   # CL3 other cancer
     0.0,  0.1,   # CL4 other cancer
     0.2,  0.0,   # CL5 non-cancer
     0.1,  0.1    # CL6 non-cancer
  ),
  nrow = 6, byrow = TRUE,
  dimnames = list(
    paste0("ACH-", 1:6),
    c("KRAS", "A1BG")
  )
)
meta <- data.frame(
  ModelID = paste0("ACH-", 1:6),
  OncotreeSubtype = c(
    "Pancreatic Adenocarcinoma", "Pancreatic Adenocarcinoma",
    "Melanoma", "Melanoma",
    "Normal", "Normal"
  ),
  OncotreePrimaryDisease = c(
    "Pancreatic Cancer", "Pancreatic Cancer",
    "Melanoma", "Melanoma",
    "Non-Cancerous", "Non-Cancerous"
  ),
  StrippedCellLineName = paste0("LINE", 1:6),
  stringsAsFactors = FALSE
)

g_kras <- cpt_group_comparison_df(
  mat, meta, "KRAS", "Pancreatic Adenocarcinoma"
)
assert("returns data.frame for valid gene", is.data.frame(g_kras) && nrow(g_kras) == 6)
assert("has required columns",
       all(c("score", "group", "cell_line") %in% colnames(g_kras)))
assert("cancer of interest has 2 lines",
       sum(g_kras$group == "Cancer of Interest") == 2)
assert("non-cancer has 2 lines",
       sum(g_kras$group == "Non-Cancer") == 2)
assert("other cancers has 2 lines",
       sum(g_kras$group == "Other Cancers") == 2)
assert("cancer scores match matrix for KRAS",
       all.equal(
         sort(g_kras$score[g_kras$group == "Cancer of Interest"]),
         sort(c(-1.2, -1.0)),
         tolerance = 1e-9
       ))

g_a1bg <- cpt_group_comparison_df(mat, meta, "A1BG", "Pancreatic Adenocarcinoma")
assert("different gene yields different scores (isolation)",
       !isTRUE(all.equal(
         g_kras$score[g_kras$group == "Cancer of Interest"],
         g_a1bg$score[g_a1bg$group == "Cancer of Interest"]
       )))
assert("A1BG cancer scores from correct column",
       all.equal(
         sort(g_a1bg$score[g_a1bg$group == "Cancer of Interest"]),
         sort(c(-0.1, 0.0)),
         tolerance = 1e-9
       ))

assert("NULL gene returns NULL", is.null(cpt_group_comparison_df(mat, meta, NULL, "Melanoma")))
assert("empty gene returns NULL", is.null(cpt_group_comparison_df(mat, meta, "", "Melanoma")))
assert("unknown gene returns NULL", is.null(cpt_group_comparison_df(mat, meta, "NOTREAL", "Melanoma")))
assert("NULL mat returns NULL", is.null(cpt_group_comparison_df(NULL, meta, "KRAS", "Melanoma")))
assert("NULL meta returns NULL", is.null(cpt_group_comparison_df(mat, NULL, "KRAS", "Melanoma")))

# Regression for the old shared-reactive bug: calling with gene A then gene B
# must not contaminate results (function is pure / argument-driven).
g1 <- cpt_group_comparison_df(mat, meta, "KRAS", "Pancreatic Adenocarcinoma")
g2 <- cpt_group_comparison_df(mat, meta, "A1BG", "Pancreatic Adenocarcinoma")
g1b <- cpt_group_comparison_df(mat, meta, "KRAS", "Pancreatic Adenocarcinoma")
assert("repeated KRAS call identical (pure function)",
       identical(g1$score, g1b$score) && identical(g1$group, g1b$group))
assert("KRAS and A1BG remain distinct after interleaved calls",
       !identical(g1$score, g2$score))

# ---------------------------------------------------------------------------
cat("\n[cpt_boiled_egg_gg axes]\n")
# ---------------------------------------------------------------------------
adme_toy <- data.frame(
  probe_name = c("AC1", "AC2", "AC3"),
  WLOGP = c(1.0, 2.5, 4.0),
  TPSA = c(40, 80, 120),
  stringsAsFactors = FALSE
)

p_all <- cpt_boiled_egg_gg(adme_toy, selected_probe = "AC2", show_mode = "all")
assert("returns ggplot", inherits(p_all, "ggplot"))
assert("x axis label is WLOGP", identical(p_all$labels$x, "WLOGP"))
assert("y axis label is TPSA", grepl("TPSA", p_all$labels$y, fixed = TRUE))
assert("title is BOILED-Egg", identical(p_all$labels$title, "BOILED-Egg"))
assert("subtitle explains HIA white / BBB yolk",
       !is.null(p_all$labels$subtitle) &&
         grepl("HIA", p_all$labels$subtitle, fixed = TRUE) &&
         grepl("BBB", p_all$labels$subtitle, fixed = TRUE))

# Coordinate limits: WLOGP on x (-3..7), TPSA on y (0..180)
xlim <- NULL
ylim <- NULL
for (s in p_all$scales$scales) {
  if (inherits(s, "ScaleContinuousPosition")) {
    if ("x" %in% s$aesthetics) xlim <- s$limits
    if ("y" %in% s$aesthetics) ylim <- s$limits
  }
}
# coord_cartesian stores limits on coordinates, not always scales
if (!is.null(p_all$coordinates$limits)) {
  if (!is.null(p_all$coordinates$limits$x)) xlim <- p_all$coordinates$limits$x
  if (!is.null(p_all$coordinates$limits$y)) ylim <- p_all$coordinates$limits$y
}
assert("x limits look like WLOGP range (not 0-180)",
       !is.null(xlim) && xlim[1] < 0 && xlim[2] < 50)
assert("y limits look like TPSA range (0-180)",
       !is.null(ylim) && ylim[1] <= 0 && ylim[2] >= 150)

# Build data layers: point layer aes must map x=WLOGP y=TPSA
point_layers <- Filter(function(ly) inherits(ly$geom, "GeomPoint"), p_all$layers)
assert("has point layer(s)", length(point_layers) >= 1)
# Extract mapping from first point layer
pm <- point_layers[[1]]$mapping
# aes may be quosures; as.character of rlang labels
x_map <- if (!is.null(pm$x)) rlang::as_label(pm$x) else NA_character_
y_map <- if (!is.null(pm$y)) rlang::as_label(pm$y) else NA_character_
# Fallback without rlang
if (is.na(x_map) || x_map == "NA") {
  x_map <- as.character(pm$x)[1]
  y_map <- as.character(pm$y)[1]
}
assert("point layer x aesthetic is WLOGP", grepl("WLOGP", x_map, fixed = TRUE))
assert("point layer y aesthetic is TPSA", grepl("TPSA", y_map, fixed = TRUE))

# Polygon egg layers also WLOGP x TPSA y
poly_layers <- Filter(function(ly) inherits(ly$geom, "GeomPolygon"), p_all$layers)
assert("has egg polygon layers", length(poly_layers) >= 2)
poly_x <- rlang::as_label(poly_layers[[1]]$mapping$x)
poly_y <- rlang::as_label(poly_layers[[1]]$mapping$y)
assert("egg polygon x is WLOGP", grepl("WLOGP", poly_x, fixed = TRUE))
assert("egg polygon y is TPSA", grepl("TPSA", poly_y, fixed = TRUE))

# selected_only mode still builds
p_sel <- cpt_boiled_egg_gg(adme_toy, selected_probe = "AC2", show_mode = "selected_only")
assert("selected_only returns ggplot", inherits(p_sel, "ggplot"))
assert("selected_only keeps WLOGP x", identical(p_sel$labels$x, "WLOGP"))

# missing columns
p_bad <- cpt_boiled_egg_gg(
  data.frame(probe_name = "X", foo = 1),
  selected_probe = "X"
)
assert("missing WLOGP/TPSA still returns ggplot", inherits(p_bad, "ggplot"))

# alias columns
adme_alias <- data.frame(
  probe_name = "P1",
  wlogp = 1.5,
  tpsa = 50,
  stringsAsFactors = FALSE
)
p_alias <- cpt_boiled_egg_gg(adme_alias, "P1")
assert("accepts lowercase wlogp/tpsa aliases", inherits(p_alias, "ggplot"))
assert("alias plot still labels x as WLOGP", identical(p_alias$labels$x, "WLOGP"))

# ---------------------------------------------------------------------------
cat("\n[process_swissadme_data]\n")
# ---------------------------------------------------------------------------
empty_swiss <- process_swissadme_data(tempdir())
assert("missing swissadme rds returns empty tibble/df",
       is.data.frame(empty_swiss) && nrow(empty_swiss) == 0)

# ---------------------------------------------------------------------------
cat("\n[source wiring / static regressions]\n")
# ---------------------------------------------------------------------------
# Modules must call pure helpers and tab-specific gene inputs (not shared reactive)
deps_src <- paste(readLines("R/dependencies_module.R", warn = FALSE), collapse = "\n")
swiss_src <- paste(readLines("R/swissadme_module.R", warn = FALSE), collapse = "\n")
app_src <- paste(readLines("app.R", warn = FALSE), collapse = "\n")

assert("deps uses cpt_group_comparison_df",
       grepl("cpt_group_comparison_df", deps_src, fixed = TRUE))
assert("cancer plot uses selected_gene_cancer",
       grepl("make_group_gg\\(input\\$selected_gene_cancer\\)", deps_src))
assert("all-genes plot uses selected_gene_all",
       grepl("make_group_gg\\(input\\$selected_gene_all\\)", deps_src))
# Old bug: one reactive preferred cancer picker then fell back to all.
old_shared_picker <- grepl(
  "if\\s*\\(!is\\.null\\(input\\$selected_gene_cancer\\)",
  deps_src
) && grepl(
  "nzchar\\(input\\$selected_gene_cancer\\)",
  deps_src
)
assert("no shared prefer-cancer-then-all gene picker", !isTRUE(old_shared_picker))

count_re <- function(pattern, text) {
  m <- gregexpr(pattern, text)[[1]]
  ml <- attr(m, "match.length")
  if (length(ml) == 1L && ml[1] == -1L) 0L else length(m)
}
n_cancer_calls <- count_re("make_group_gg\\(input\\$selected_gene_cancer\\)", deps_src)
n_all_calls <- count_re("make_group_gg\\(input\\$selected_gene_all\\)", deps_src)
assert("make_group_gg called for cancer gene at least twice (plot+download)",
       n_cancer_calls >= 2)
assert("make_group_gg called for all gene at least twice (plot+download)",
       n_all_calls >= 2)

assert("HTML report gene follows the open Dependencies subtab",
       grepl("resolve_report_gene", deps_src, fixed = TRUE) &&
         grepl("input\\$dep_subtabs", deps_src) &&
         grepl("input\\$selected_gene_all", deps_src) &&
         grepl("input\\$selected_gene_cancer", deps_src))

assert("swissadme uses cpt_boiled_egg_gg",
       grepl("cpt_boiled_egg_gg", swiss_src, fixed = TRUE))
# Ensure we did not reintroduce swapped aes in the module itself
assert("swissadme module has no aes(x = TPSA, y = WLOGP)",
       !grepl("aes\\(x\\s*=\\s*TPSA\\s*,\\s*y\\s*=\\s*WLOGP", swiss_src))

assert("app.R defines protein_binding_lookup key",
       grepl("protein_binding_lookup", app_src, fixed = TRUE))
assert("app.R keeps proteinbindinglookup alias",
       grepl("proteinbindinglookup", app_src, fixed = TRUE))

# precompute script must not redefine sanitize_subtype locally
pre_src <- paste(readLines("docs/scripts/precompute_effectsizes_all_subtypes.R", warn = FALSE),
                 collapse = "\n")
assert("precompute sources functions.R",
       grepl('source\\(file\\.path\\("R",\\s*"functions\\.R"', pre_src) ||
         grepl("functions.R", pre_src, fixed = TRUE))
assert("precompute does not redefine sanitize_subtype",
       !grepl("sanitize_subtype\\s*<-\\s*function", pre_src))

# ---------------------------------------------------------------------------
cat("\n[live data smoke: group comparison on real matrix]\n")
# ---------------------------------------------------------------------------
if (file.exists("data/CRISPRGeneEffect_23Q4_clean.rds") &&
    file.exists("data/cancer_model_data.rds")) {
  # Small slice only — full matrix is large; take a few columns after load head
  # Load full is heavy but already used by other tests; use on-the-fly subset.
  mat_live <- readRDS("data/CRISPRGeneEffect_23Q4_clean.rds")
  meta_live <- readRDS("data/cancer_model_data.rds")
  # Column names are often "GENE (ENTREZ)"
  gene_cols <- colnames(mat_live)
  kras_col <- gene_cols[grepl("^KRAS\\b", gene_cols)][1]
  braf_col <- gene_cols[grepl("^BRAF\\b", gene_cols)][1]
  if (!is.na(kras_col) && !is.na(braf_col)) {
    g_k <- cpt_group_comparison_df(
      mat_live, meta_live, kras_col, "Pancreatic Adenocarcinoma"
    )
    g_b <- cpt_group_comparison_df(
      mat_live, meta_live, braf_col, "Melanoma"
    )
    assert("live KRAS/PDAC group df non-empty",
           is.data.frame(g_k) && nrow(g_k) > 10)
    assert("live BRAF/Melanoma group df non-empty",
           is.data.frame(g_b) && nrow(g_b) > 10)
    assert("live KRAS cancer group has rows",
           sum(g_k$group == "Cancer of Interest") >= 1)
    # Different genes => different mean cancer scores (almost always)
    mk <- mean(g_k$score[g_k$group == "Cancer of Interest"], na.rm = TRUE)
    mb <- mean(g_b$score[g_b$group == "Cancer of Interest"], na.rm = TRUE)
    assert("live KRAS and BRAF cancer means differ",
           is.finite(mk) && is.finite(mb) && abs(mk - mb) > 1e-6)
  } else {
    cat("  SKIP: KRAS/BRAF columns not found in matrix\n")
  }
} else {
  cat("  SKIP: live DepMap RDS not present\n")
}

# ---------------------------------------------------------------------------
cat(sprintf("\n=== Results: %d/%d passed", tests_passed, tests_run))
if (tests_failed > 0) {
  cat(sprintf(" (%d FAILED) ===\n", tests_failed))
  quit(status = 1)
}
cat(" ===\n")
