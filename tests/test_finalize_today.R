#!/usr/bin/env Rscript --vanilla
# ============================================================
# Test: finalize-today features (1.5 data versions, CPT table
# annotate, radar, weight sensitivity, compare wiring).
# Run: Rscript --vanilla tests/test_finalize_today.R
# ============================================================

cat("=== Finalize-today tests ===\n\n")

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
source("R/canprotarget_score.R", local = TRUE)
source("R/functions.R", local = TRUE)
source("R/api_functions.R", local = TRUE)

# ---- 1.5 data versions ----------------------------------------
cat("[data_versions]\n")
dv <- read_data_versions("data")
assert("reads data_versions.yaml", !is.null(dv) && !is.null(dv$datasets))
tbl <- data_versions_table(dv)
assert("versions table has rows", is.data.frame(tbl) && nrow(tbl) >= 4)
assert("CRISPR dataset listed", any(grepl("CRISPR", tbl$name)))
assert("every dataset carries a status", all(nzchar(tbl$status)))
assert("chemoproteomics + SwissADME now available",
       all(tbl$status[grepl("Chemoproteomics|SwissADME", tbl$name)] == "available"))
prov <- data_provenance_string(dv)
assert("provenance string non-empty", nzchar(prov) && !grepl("^Data versions not", prov))
assert("provenance includes chemoproteomics now that it is available",
       grepl("Chemoproteomics", prov))
ui <- cpt_data_versions_ui(dv)
assert("About UI is shiny tag / tagList",
       inherits(ui, "shiny.tag") || inherits(ui, "shiny.tag.list") ||
         inherits(ui, "list"))

# ---- CPT annotate table ---------------------------------------
cat("\n[cpt_annotate_gene_table]\n")
bg <- data.frame(
  gene_name = c("A", "B", "C", "D"),
  Cancer_Avg = c(-1.2, -0.8, -0.1, 0.2),
  EffectSize = c(-0.9, -0.5, -0.05, 0.1),
  stringsAsFactors = FALSE
)
filt <- bg[1:2, , drop = FALSE]
ann <- cpt_annotate_gene_table(filt, bg, cys_atlas = NULL)
assert("annotates CPT_Score", "CPT_Score" %in% colnames(ann) && all(!is.na(ann$CPT_Score)))
assert("keeps both genes", nrow(ann) == 2)
assert("stronger dep gene ranks better or equal",
       ann$CPT_Score[ann$gene_name == "A"] >= ann$CPT_Score[ann$gene_name == "B"] - 1e-6)
assert("dep percentiles not all 50",
       !all(abs(ann$dim_dependency - 50) < 1e-6))

# ---- radar + sensitivity --------------------------------------
cat("\n[radar + weight sensitivity]\n")
dim_scores <- c(
  dependency_strength = 90,
  cancer_selectivity = 80,
  cysteine_ligandability = 60,
  conservation = NA_real_,
  clinical_evidence = 40,
  adme_druggability = NA_real_
)
radar <- cpt_dimension_radar_gg(dim_scores, title = "Test")
assert("radar is ggplot", inherits(radar, "ggplot"))
sens <- cpt_weight_sensitivity(dim_scores)
assert("sensitivity has baseline row", any(sens$dimension == "(baseline)"))
assert("sensitivity has multiple rows", nrow(sens) >= 4)
assert("scores stay in 0-100-ish",
       all(sens$cpt_score >= 0 & sens$cpt_score <= 100, na.rm = TRUE))
# Perturbing dependency weight should move score vs baseline
base <- sens$cpt_score[sens$dimension == "(baseline)"][1]
dep_hi <- sens$cpt_score[sens$dimension == "dependency_strength" & sens$factor == 2]
assert("x2 dependency weight changes composite",
       length(dep_hi) == 1 && abs(dep_hi - base) > 1e-6)

# Radar omits NA dims (subtitle should mention missing when any NA)
assert("radar subtitle mentions missing dims",
       grepl("missing", radar$labels$subtitle, ignore.case = TRUE))

# ---- priority guards + dataset normalize + atlas fairness ----
cat("\n[priority guards + dataset + ligandability NA]\n")
assert("high demoted when not dep",
       identical(cpt_apply_priority_guards("high", FALSE, TRUE), "moderate"))
assert("moderate demoted when not dep and not selective",
       identical(cpt_apply_priority_guards("moderate", FALSE, FALSE), "low"))
assert("high stays when dependency",
       identical(cpt_apply_priority_guards("high", TRUE, TRUE), "high"))
assert("normalize CRISPR ok", identical(cpt_normalize_dataset("CRISPR"), "CRISPR"))
assert("normalize RNAi ok", identical(cpt_normalize_dataset("RNAi"), "RNAi"))
assert("normalize CRISPR (23Q4) ok",
       identical(cpt_normalize_dataset("CRISPR (23Q4)"), "CRISPR"))
assert("clamp negative n to 1",
       identical(cpt_clamp_int(-5, 1L, 100L, 20L), 1L))
assert("clamp huge n to max",
       identical(cpt_clamp_int(500, 1L, 100L, 20L), 100L))
bad_ds <- tryCatch(cpt_normalize_dataset("FOO"), error = function(e) e$message)
assert("invalid dataset errors", is.character(bad_ds) && grepl("CRISPR|RNAi", bad_ds))

# Atlas: gene not present → NA ligandability (not 0)
fake_atlas <- data.frame(
  gene_symbol = c("GENE_IN", "GENE_IN"),
  ligandable = c(FALSE, FALSE),
  clinvar_pathogenic = c(FALSE, FALSE),
  stringsAsFactors = FALSE
)
lig_missing <- score_cysteine_ligandability(c("NOT_IN_ATLAS", "GENE_IN"), fake_atlas)
assert("not-in-atlas ligandability is NA", is.na(lig_missing[1]))
assert("explicit non-ligandable is 0", identical(as.numeric(lig_missing[2]), 0))
cli_missing <- score_clinical_evidence(c("NOT_IN_ATLAS", "GENE_IN"), fake_atlas)
assert("not-in-atlas clinical is NA", is.na(cli_missing[1]))
assert("in-atlas no pathogenic is 0", identical(as.numeric(cli_missing[2]), 0))

prio_disp <- cpt_priority_display(list(
  cpt_score = 85, priority = "moderate", is_dependency = FALSE,
  is_selective_by_effect_size = TRUE
))
assert("display uses guarded priority not raw score",
       grepl("Moderate", prio_disp, ignore.case = TRUE))

# ---- compare subtypes (live data) -----------------------------
cat("\n[compare + recovery fixtures]\n")
if (file.exists("data/CRISPRGeneEffect_23Q4_clean.rds") &&
    file.exists("data/cancer_model_data.rds")) {
  data_env <- list(data_dir = "data")
  data_env$crispr_matrix <- readRDS("data/CRISPRGeneEffect_23Q4_clean.rds")
  data_env$cancer_model_data <- readRDS("data/cancer_model_data.rds")
  data_env$cys_atlas <- if (file.exists("data/cys_editing_atlas.rds")) {
    readRDS("data/cys_editing_atlas.rds")
  } else NULL

  # Note: Oncotree "Colorectal Adenocarcinoma" can have <3 CRISPR lines;
  # use Colon Adenocarcinoma (n large enough) for recovery checks.
  cmp <- api_compare_subtypes(
    "BRAF", "Melanoma", "Colon Adenocarcinoma", "CRISPR", data_env
  )
  assert("compare returns gene BRAF", identical(cmp$gene, "BRAF"))
  assert("compare has two subtype blocks",
         !is.null(cmp$subtype1$mean_effect) && !is.null(cmp$subtype2$mean_effect))
  assert("melanoma mean more negative than colon adeno for BRAF (typical)",
         cmp$subtype1$mean_effect < cmp$subtype2$mean_effect)

  # Recovery: KRAS PDAC
  sc <- api_canprotarget_score(
    "KRAS", "Pancreatic Adenocarcinoma", "CRISPR", data_env
  )
  assert("KRAS PDAC has CPT score", is.numeric(sc$cpt_score) && !is.na(sc$cpt_score))
  assert("KRAS PDAC is dependency", isTRUE(sc$is_dependency))
  assert("KRAS PDAC CPT above floor", sc$cpt_score >= 40)

  # EGFR_797
  if (!is.null(data_env$cys_atlas)) {
    det <- api_cysteine_detail("EGFR_797", data_env)
    assert("EGFR_797 found", !is.null(det) && !is.null(det$site_id))
    assert("EGFR_797 functional", isTRUE(det$functional) || isTRUE(det$is_functional))
  }

  # annotate real cancer-like slice
  q_kras <- api_query_dependency("KRAS", "Pancreatic Adenocarcinoma", "CRISPR", data_env)
  q_a1bg <- api_query_dependency("A1BG", "Pancreatic Adenocarcinoma", "CRISPR", data_env)
  bg_live <- data.frame(
    gene_name = c("KRAS", "A1BG"),
    Cancer_Avg = c(q_kras$cancer_mean_effect, q_a1bg$cancer_mean_effect),
    EffectSize = c(q_kras$effect_size, q_a1bg$effect_size),
    stringsAsFactors = FALSE
  )
  # Expand background with a few more synthetic rows so percentiles differ
  bg_big <- rbind(
    bg_live,
    data.frame(gene_name = paste0("G", 1:20),
               Cancer_Avg = seq(-0.4, 0.3, length.out = 20),
               EffectSize = seq(-0.2, 0.2, length.out = 20),
               stringsAsFactors = FALSE)
  )
  ann_live <- cpt_annotate_gene_table(bg_live, bg_big, data_env$cys_atlas)
  assert("live annotate KRAS score > A1BG",
         ann_live$CPT_Score[ann_live$gene_name == "KRAS"] >
           ann_live$CPT_Score[ann_live$gene_name == "A1BG"])
} else {
  cat("  SKIP: DepMap RDS missing\n")
}

# ---- fixtures file present ------------------------------------
cat("\n[case study fixtures file]\n")
assert("recovery fixtures exist",
       file.exists("tests/fixtures/known_target_recovery.yaml"))
assert("Target tab exposes subtype comparison",
       grepl("Compare two subtypes", paste(readLines("R/gene_module.R"), collapse = "\n")))
assert("app sources api_functions",
       grepl("api_functions", paste(readLines("app.R"), collapse = "\n")))
assert("About uses about_data_versions",
       grepl("about_data_versions", paste(readLines("app.R"), collapse = "\n")))

cat(sprintf("\n=== Results: %d/%d passed", tests_passed, tests_run))
if (tests_failed > 0) {
  cat(sprintf(" (%d FAILED) ===\n", tests_failed))
  quit(status = 1)
}
cat(" ===\n")
