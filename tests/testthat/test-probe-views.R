project_root <- normalizePath(file.path(testthat::test_path(), "..", ".."))
suppressPackageStartupMessages(library(shiny))
source(file.path(project_root, "R", "canprotarget_score.R"))
source(file.path(project_root, "R", "targets_panel.R"))
source(file.path(project_root, "R", "gene_module.R"))

index_path <- file.path(project_root, "data", "gene_index.rds")

test_that("collapsing keeps a passing probe over a stronger promiscuous one", {
  es <- data.frame(
    gene_name = "SAFB2", cysteineid = "Q14151_C672",
    probe_name = c("AC34", "CL344"), CR = c(17.78, 4.77), n_targets = c(364, 11),
    stringsAsFactors = FALSE
  )
  for (grain in c("gene", "cysteine")) {
    out <- cpt_collapse_engaged(es, grain, cr = 4, mt = 20)
    expect_equal(nrow(out), 1)
    expect_true(out$passes)
    expect_equal(out$probe_name, "CL344")
  }
  # Nothing passes: still one row, the strongest record, flagged as failing.
  out <- cpt_collapse_engaged(es, "gene", cr = 4, mt = 5)
  expect_equal(out$probe_name, "AC34")
  expect_false(out$passes)
  # Probe grain keeps every record.
  expect_equal(nrow(cpt_collapse_engaged(es, "probe", 4, 20)), 2)
})

test_that("SAFB2 passes the ligandability browse at every grain on the real index", {
  skip_if_not(file.exists(index_path))
  es <- readRDS(index_path)$engaged_sites
  es <- es[es$gene_name == "SAFB2", , drop = FALSE]
  for (grain in c("gene", "cysteine")) {
    out <- cpt_collapse_engaged(es, grain, cr = 4, mt = 20)
    expect_true(all(out$passes[out$cysteineid == "Q14151_C672"]))
    expect_lte(out$n_targets[out$cysteineid == "Q14151_C672"], 20)
  }
})

test_that("Target probe table lists every engagement, not the top 10", {
  skip_if_not(file.exists(index_path))
  idx <- readRDS(index_path)
  cr4 <- readRDS(file.path(project_root, "data", "protein_binding_cr4.rds"))
  cr4$probe_name <- cpt_canonical_probe_name(cr4$probe_name)

  df <- cpt_gene_probe_rows(idx, "SAFB2", cr4)
  expect_true("CL344" %in% df$probe_name)
  cl344 <- df[df$probe_name == "CL344", , drop = FALSE]
  expect_equal(nrow(cl344), 1)
  expect_equal(cl344$cysteineid, "Q14151_C672")
  expect_equal(cl344$CR, 4.77)
  expect_false(is.na(cl344$SMILES))
  # All 15 CR >= 4 engagements, none duplicated, sorted by CR.
  eng <- df[df$CR >= 4, , drop = FALSE]
  expect_equal(nrow(eng), 15)
  expect_false(any(duplicated(paste(df$probe_name, df$cysteineid))))
  expect_false(is.unsorted(rev(df$CR)))
})
