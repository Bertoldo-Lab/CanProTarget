#!/usr/bin/env Rscript
# ============================================================
# Cross-source gene-name matching.
#
# DepMap gene-effect matrix columns carry an Entrez suffix ("KRAS (3845)"), which
# flows into every gene_name produced by ge_analysis(). The chemoproteomics
# binding table and the cysteine atlas use bare symbols ("KRAS"). Joining the raw
# strings matches NOTHING and fails silently: the Ligandability panel and the
# ligandable volcano simply render empty, with no error.
#
# That was a live bug (Acral Melanoma showed zero ligandable proteins locally
# while the deployed app showed several). These tests lock the fix in.
#
# Run from project root:  Rscript --vanilla tests/test_gene_key_matching.R
# ============================================================

suppressPackageStartupMessages({ library(dplyr) })

pass <- 0L; fail <- 0L; failed <- character(0)
assert <- function(label, cond) {
  if (isTRUE(cond)) { pass <<- pass + 1L; cat("  PASS:", label, "\n") }
  else { fail <<- fail + 1L; failed <<- c(failed, label); cat("  FAIL:", label, "\n") }
}

source("R/functions.R", local = TRUE)

cat("\n[cpt_gene_match_key]\n")
assert("helper exists", exists("cpt_gene_match_key", mode = "function"))
assert("strips the Entrez suffix",
       identical(cpt_gene_match_key("KRAS (3845)"), "KRAS"))
assert("strips multi-digit Entrez ids",
       identical(cpt_gene_match_key("AGAP3 (116988)"), "AGAP3"))
assert("leaves bare symbols untouched",
       identical(cpt_gene_match_key("KRAS"), "KRAS"))
assert("upper-cases for case-insensitive matching",
       identical(cpt_gene_match_key("kras"), "KRAS"))
assert("is vectorised",
       identical(cpt_gene_match_key(c("A1BG (1)", "egfr", "TP53")),
                 c("A1BG", "EGFR", "TP53")))
assert("keeps hyphenated symbols intact",
       identical(cpt_gene_match_key("HLA-A"), "HLA-A"))
assert("only strips a TRAILING parenthesised number",
       identical(cpt_gene_match_key("ABC (2) X"), "ABC (2) X"))
assert("tolerates NA", is.na(cpt_gene_match_key(NA)))

cat("\n[all five cross-source join sites use the helper]\n")
mod <- paste(readLines("R/dependencies_module.R", warn = FALSE), collapse = "\n")
cys <- paste(readLines("R/cys_editing_functions.R", warn = FALSE), collapse = "\n")
cmd <- paste(readLines("R/cys_editing_module.R", warn = FALSE), collapse = "\n")

# 1. probe table filtered to cancer dependency genes
# Matched on the helper being applied at the join, not on one spelling of it:
# the filter now compares against toupper(keep_keys).
assert("1. probe cancer-gene filter is normalised",
       grepl("cpt_gene_match_key\\(gene_name\\) %in% [a-z(]*keep_keys", mod))
assert("   no raw join against cg_df remains",
       !grepl("filter(gene_name %in% cg_df$gene_name)", mod, fixed = TRUE))
# 2. ligandable volcano (probe genes ∩ dependency volcano)
assert("2. ligandable volcano join is normalised",
       grepl("cpt_gene_match_key(gene_name) %in% lig_keys", mod, fixed = TRUE))
assert("   no raw join against lig_genes remains",
       !grepl("filter(gene_name %in% lig_genes)", mod, fixed = TRUE))
# 3. combined dependency + probe table (feeds the selectivity scatter)
assert("3. combined table joins on a normalised key",
       grepl(".cpt_gene_key = cpt_gene_match_key(gene_name)", mod, fixed = TRUE))
assert("   no raw left_join by gene_name remains",
       !grepl('left_join(selective_probes_per_protein, by = "gene_name")',
              mod, fixed = TRUE))
# 4. Functional Cysteines dependency intersection
assert("4. cys module normalises incoming dependency genes",
       grepl("cpt_gene_match_key(dependency_genes())", cmd, fixed = TRUE))
assert("   no bare toupper() of dependency genes remains",
       !grepl("unique(toupper(as.character(dependency_genes())))", cmd, fixed = TRUE))
# 5. atlas annotation + functional-site linking
assert("5a. annotate_probe_table normalises",
       grepl("cpt_gene_match_key(out$gene_name)", cys, fixed = TRUE))
assert("5b. link_functional_sites joins on a normalised key",
       grepl('by = ".cpt_gene_key"', cys, fixed = TRUE))
assert("    no raw inner_join by gene_name remains",
       !grepl('by = "gene_name",\n    relationship', cys, fixed = TRUE))
assert("cys_editing_functions.R is self-contained when sourced standalone",
       grepl('if (!exists("cpt_gene_match_key", mode = "function"))', cys, fixed = TRUE))

# ---- live data: the join must actually return rows ---------------
bind_p <- "data/protein_binding_lookup_preprocessed.rds"
mat_p  <- "data/CRISPRGeneEffect_23Q4_clean.rds"

if (!file.exists(bind_p) || !file.exists(mat_p)) {
  cat("\n[live data] SKIP - chemoproteomics or CRISPR matrix absent\n")
} else {
  cat("\n[live cross-source join]\n")
  pb <- readRDS(bind_p)
  mat <- readRDS(mat_p)

  dep_names <- colnames(mat)          # as ge_analysis() emits them
  probe_names <- unique(pb$gene_name) # bare symbols

  # The suffixed convention is what made the normaliser necessary, but the
  # matrices are not always built that way: the 23Q4 rebuild emits bare
  # symbols. Assert the property that must hold either way -- the normalised
  # join recovers the genes -- and only assert the degenerate raw join when
  # the labels actually differ between the two sources.
  suffixed <- any(grepl("\\(\\d+\\)$", dep_names))
  cat(sprintf("       matrix labels carry an Entrez suffix: %s\n", suffixed))
  assert("probe gene names really are bare",
         !any(grepl("\\(\\d+\\)$", probe_names)))

  raw <- sum(probe_names %in% dep_names)
  fixed <- length(intersect(cpt_gene_match_key(probe_names),
                            cpt_gene_match_key(dep_names)))
  cat(sprintf("       raw string overlap = %d, normalised overlap = %d\n", raw, fixed))
  if (suffixed) {
    assert("raw join is the degenerate case it was", raw == 0L)
  } else {
    assert("raw join already works when both sides are bare", raw > 5000L)
  }
  assert("normalised join recovers thousands of genes", fixed > 5000L)

  cat("\n[regression guard: a known ligandable dependency survives the join]\n")
  # CCND1 is the case that exposed the bug on Acral Melanoma.
  assert("CCND1 is present in the probe table",
         "CCND1" %in% cpt_gene_match_key(probe_names))
  assert("CCND1 is present in the CRISPR matrix",
         "CCND1" %in% cpt_gene_match_key(dep_names))
  assert("CCND1 survives the normalised join",
         "CCND1" %in% intersect(cpt_gene_match_key(probe_names),
                                cpt_gene_match_key(dep_names)))

  # ---- the atlas linkers must tolerate BOTH label conventions -------
  atlas_p <- "data/cys_editing_atlas.rds"
  if (file.exists(atlas_p)) {
    cat("\n[atlas linkers accept either label convention]\n")
    source("R/cys_editing_functions.R", local = TRUE)
    atlas <- readRDS(atlas_p)

    lig <- pb[toupper(as.character(pb$ligandable)) %in% c("YES", "TRUE", "1", "Y"), ,
              drop = FALSE]
    bare <- head(lig, 4000)
    sfx <- bare
    sfx$gene_name <- paste0(sfx$gene_name, " (999)")   # force the suffixed form

    ann_b <- cpt_cys_annotate_probe_table(bare, atlas)
    ann_s <- cpt_cys_annotate_probe_table(sfx, atlas)
    assert("annotate: bare labels map",
           sum(ann_b$cys_site_in_atlas, na.rm = TRUE) > 0)
    assert("annotate: suffixed labels map identically",
           identical(sum(ann_b$cys_site_in_atlas, na.rm = TRUE),
                     sum(ann_s$cys_site_in_atlas, na.rm = TRUE)))

    lk_b <- cpt_cys_link_functional_sites(bare, atlas)
    lk_s <- cpt_cys_link_functional_sites(sfx, atlas)
    assert("link: bare labels produce rows", nrow(lk_b) > 0)
    assert("link: suffixed labels produce the same row count",
           nrow(lk_b) == nrow(lk_s))
    assert("link: display gene_name is preserved, not clobbered",
           "gene_name" %in% names(lk_s) &&
             all(grepl("\\(999\\)$", lk_s$gene_name)))
    assert("link: helper key column is not leaked into output",
           !(".cpt_gene_key" %in% names(lk_s)))
  }
}

cat(sprintf("\n=== Results: %d/%d passed ===\n", pass, pass + fail))
if (fail > 0) {
  cat("Failed tests:\n"); for (f in failed) cat("  -", f, "\n")
  quit(status = 1)
}
