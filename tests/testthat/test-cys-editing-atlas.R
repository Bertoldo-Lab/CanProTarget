project_root <- normalizePath(file.path(testthat::test_path(), "..", ".."))
source(file.path(project_root, "R", "cys_editing_functions.R"))

test_that("within-gene FDR is calculated separately for each protein", {
  result <- cpt_cys_within_gene_fdr(
    neg_log10_p = c(2, 1, 2),
    gene_symbol = c("A", "A", "B")
  )

  expected_a <- round(-log10(p.adjust(c(0.01, 0.1), method = "fdr")), 2)
  expected_b <- 2
  expect_equal(result, c(expected_a, expected_b))
})

test_that("bundled atlas preserves the publication-derived calls", {
  atlas_path <- file.path(project_root, "data", "cys_editing_atlas.rds")
  expect_true(file.exists(atlas_path))

  atlas <- readRDS(atlas_path)
  expect_equal(nrow(atlas), 13872)
  expect_equal(length(unique(atlas$gene_symbol)), 1778)
  expect_equal(sum(atlas$functional), 1718)
  expect_equal(sum(atlas$functional_ligandable), 160)
  expect_equal(attr(atlas, "source_commit"), "89bc6a268ea4385a687f2e193046b4486edd0185")

  egfr <- atlas[atlas$site_id == "EGFR_797", , drop = FALSE]
  expect_equal(nrow(egfr), 1)
  expect_true(egfr$functional)
  expect_true(egfr$ligandable)
  expect_equal(egfr$ligandability_score, 100)
  expect_true(egfr$ligandability_manual_annotation)
})

test_that("residue mapping status is explicit for every site", {
  atlas <- readRDS(file.path(project_root, "data", "cys_editing_atlas.rds"))
  expect_false(anyNA(atlas$residue_mapping_status))
  expect_setequal(
    unique(atlas$residue_mapping_status),
    c(
      "Source sequence matched",
      "Source sequence mismatch",
      "Unmapped to source AlphaFold table"
    )
  )
})

test_that("probe rows are joined to the exact functional cysteine site", {
  atlas <- readRDS(file.path(project_root, "data", "cys_editing_atlas.rds"))
  site <- atlas[atlas$functional & !is.na(atlas$uniprot_accession), , drop = FALSE][1, ]
  probes <- data.frame(
    gene_name = c(site$gene_symbol, site$gene_symbol, "NOT_A_GENE"),
    proteinid = c(site$uniprot_accession, site$uniprot_accession, "P00000"),
    cysteineid = c(
      paste0(site$uniprot_accession, "_C", site$cysteine_position),
      paste0(site$uniprot_accession, "_C999999"),
      "P00000_C1"
    ),
    stringsAsFactors = FALSE
  )

  joined <- cpt_cys_annotate_probe_table(probes, atlas)

  expect_identical(joined$cys_site_id[[1]], site$site_id)
  expect_true(joined$cys_site_in_atlas[[1]])
  expect_true(joined$cys_functional[[1]])
  expect_identical(joined$cys_join_status[[1]], "UniProt + cysteine position")
  expect_false(joined$cys_site_in_atlas[[2]])
  expect_false(joined$cys_site_in_atlas[[3]])
})

test_that("gene-position fallback is labelled and never becomes a gene-only join", {
  atlas <- readRDS(file.path(project_root, "data", "cys_editing_atlas.rds"))
  site <- atlas[atlas$functional, , drop = FALSE][1, ]
  probes <- data.frame(
    gene_name = c(site$gene_symbol, site$gene_symbol),
    cysteineid = c(
      paste0("UNMAPPED_C", site$cysteine_position),
      paste0("UNMAPPED_C", site$cysteine_position + 1L)
    ),
    stringsAsFactors = FALSE
  )

  joined <- cpt_cys_annotate_probe_table(probes, atlas)

  expect_identical(joined$cys_site_id[[1]], site$site_id)
  expect_identical(joined$cys_join_status[[1]], "Gene + cysteine position")
  expect_false(joined$cys_site_in_atlas[[2]])
})

test_that("functional-site linking distinguishes exact-site and same-gene evidence", {
  atlas <- readRDS(file.path(project_root, "data", "cys_editing_atlas.rds"))
  candidates <- atlas[atlas$functional & !duplicated(atlas$gene_symbol), , drop = FALSE]
  site <- candidates[1, ]
  probes <- data.frame(
    gene_name = c(site$gene_symbol, site$gene_symbol),
    cysteineid = c(
      paste0("UNIPROT_C", site$cysteine_position),
      paste0("UNIPROT_C", site$cysteine_position + 100000L)
    ),
    stringsAsFactors = FALSE
  )

  linked <- cpt_cys_link_functional_sites(probes, atlas)
  linked <- linked[linked$cys_site_id == site$site_id, , drop = FALSE]

  expect_equal(nrow(linked), 2)
  expect_setequal(linked$cys_probe_site_match, c(TRUE, FALSE))
  expect_setequal(
    linked$cys_relationship,
    c("Probe binds this functional site", "Functional site in same dependency gene")
  )
})
