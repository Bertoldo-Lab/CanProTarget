project_root <- normalizePath(file.path(testthat::test_path(), "..", ".."))
suppressPackageStartupMessages(library(shiny))
source(file.path(project_root, "R", "functions.R"))
source(file.path(project_root, "R", "api_functions.R"))
source(file.path(project_root, "R", "canprotarget_score.R"))
source(file.path(project_root, "R", "report_generator.R"))
source(file.path(project_root, "R", "gene_module.R"))

data_dir <- file.path(project_root, "data")
mb_path <- file.path(data_dir, "precomputed_effectsizes", "RNAi_Medulloblastoma.rds")

test_that("dependency selection honours every Dependency-box setting", {
  df <- data.frame(
    gene_name = c("SEL", "ESS", "WEAK", "NOTSIG", "NCFAIL"),
    EffectSize = c(-0.3, -0.3, -0.05, -0.3, -0.3),
    p_value = c(1e-6, 1e-6, 1e-6, 0.2, 1e-6),
    Avg = c(-0.1, -0.8, -0.1, -0.1, -0.1),
    pval_vs_NonCancer = c(0.01, 0.01, 0.01, 0.01, 0.3),
    stringsAsFactors = FALSE
  )
  sel <- function(...) df$gene_name[cpt_dependency_status(df, -0.1, ...)$selected]
  expect_equal(sel(TRUE, TRUE, FALSE), c("SEL", "NCFAIL"))
  expect_equal(sel(TRUE, FALSE, FALSE), c("SEL", "ESS", "NCFAIL"))
  expect_equal(sel(FALSE, TRUE, FALSE), c("SEL", "NOTSIG", "NCFAIL"))
  expect_equal(sel(TRUE, TRUE, TRUE), "SEL")
  st <- cpt_dependency_status(df, -0.1, TRUE, TRUE, FALSE)
  expect_equal(st$status[st$gene_name == "ESS"], "Common essential")
})

test_that("medulloblastoma RNAi selection matches the Dependency analysis", {
  skip_if_not(file.exists(mb_path))
  df <- cpt_read_effectsizes(mb_path)
  on <- cpt_dependency_status(df, -0.1, TRUE, TRUE, FALSE)
  off <- cpt_dependency_status(df, -0.1, TRUE, FALSE, FALSE)
  expect_equal(sum(on$selected), 1154)
  expect_equal(sum(off$selected), 1199)
  expect_true(on$selected[on$gene_name == "SAFB2"])
})

test_that("Target subtype count is per screen", {
  idx_path <- file.path(data_dir, "gene_index.rds")
  skip_if_not(file.exists(idx_path))
  idx <- readRDS(idx_path)
  expect_equal(cpt_dep_n_subtypes(idx, "SAFB2", "RNAi"), 6L)
  expect_equal(cpt_dep_n_subtypes(idx, "SAFB2", "CRISPR"), 3L)
  expect_equal(cpt_dep_n_subtypes(idx, "NOT_A_GENE", "RNAi"), 0L)
})

test_that("reports render when the caller passes undeclared params", {
  skip_if_not(file.exists(file.path(data_dir, "cys_editing_atlas.rds")))
  out_dir <- tempfile("cpt_reports_")
  # The MCP server always sends report_type and dataset.
  out <- cpt_render_report(
    "cysteine_target",
    params = list(gene = "EGFR", dataset = "RNAi", report_type = "cysteine_target"),
    output_dir = out_dir, project_root = project_root
  )
  expect_true(file.exists(out))
})

test_that("SMCL index and platform info use canonical names and live counts", {
  binding <- data.frame(
    gene_name = "SAFB2", proteinid = "Q14151", cysteineid = "Q14151_C672",
    probe_name = c("CL_344", "ACRYL_216", "OTHER_6"), CR = c(4.77, 6.29, 5),
    n_targets = c(11, 11, 3), stringsAsFactors = FALSE
  )
  idx <- cpt_build_smcl_index(binding, min_cr = 4)
  expect_setequal(idx$probe_name, c("CL344", "AC216", "OTHER_6"))

  env <- list(rnai_matrix = matrix(0, nrow = 3, ncol = 5),
              crispr_matrix = matrix(0, nrow = 2, ncol = 4),
              cys_atlas = data.frame(x = 1:7))
  info <- api_platform_info(env)
  expect_equal(info$data_sources[[2]]$genes, 5L)
  expect_equal(info$data_sources[[2]]$cell_lines, 3L)
  expect_equal(info$data_sources[[1]]$genes, 4L)
  expect_equal(info$data_sources[[3]]$sites, 7L)
})

test_that("a zero-weight dimension is reported but not counted as active", {
  rep <- cpt_dimension_report(
    list(dependency_strength = 90, cancer_selectivity = 99,
         cysteine_ligandability = NA, conservation = NA,
         clinical_evidence = NA, adme_druggability = 98),
    weights = CPT_DEFAULT_WEIGHTS
  )
  expect_equal(rep$n_dimensions_used, 2L)
  expect_equal(rep$dimensions_unweighted, "adme_druggability")
  expect_false("adme_druggability" %in% rep$dimensions_missing)
  expect_equal(rep$dimensions$adme_druggability$score, 98)

  w <- CPT_DEFAULT_WEIGHTS; w$adme_druggability <- 0.5
  rep2 <- cpt_dimension_report(
    list(dependency_strength = 90, cancer_selectivity = 99, adme_druggability = 98),
    weights = w
  )
  expect_equal(rep2$n_dimensions_used, 3L)
})
