# ============================================================
# Script:   build_gene_index.R
# Purpose:  Build data/gene_index.rds — the per-gene index behind the Gene tab.
#
# Why this exists
#   Answering "what do we know about GENE X?" used to require running a
#   subtype analysis (140 MB matrix) and then loading the 10.6M-row binding
#   table (~20 s, ~0.9 GB) — for a question that has nothing to do with either.
#   This script pre-aggregates the same sources offline so the Gene tab reads
#   one small RDS and answers instantly.
#
# Inputs  (all already required by the app)
#   data/precomputed_effectsizes/{CRISPR,RNAi}_<subtype>.rds
#   data/protein_binding_lookup_preprocessed.rds
#   data/swissadme_preprocessed.rds
#   data/cys_editing_atlas.rds
#
# Output
#   data/gene_index.rds — list(dependency, dep_summary, probes, probe_summary,
#                              cys_sites_summary, adme, genes, built_at)
#
# Usage
#   Rscript docs/scripts/build_gene_index.R [data_dir]
# ============================================================

suppressMessages({
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args)) args[[1]] else "data"

source("R/app_helpers.R")
source("R/functions.R")
source("R/cys_editing_functions.R")
source("R/canprotarget_score.R")

msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n")

# Keep the index small: only rows that clear the app's own default filters,
# and at most TOP_N_SUBTYPES subtypes per gene (ranked by effect size).
EFFECT_MAX      <- -0.1
PVALUE_MAX      <- 0.05
TOP_N_SUBTYPES  <- 25
TOP_N_PROBES    <- 10

# ---- 1. Dependency: fold the 139 per-subtype TSVs into one long table -------
msg("Reading per-subtype effect sizes...")
tsv_dir <- file.path(data_dir, "precomputed_effectsizes")
tsv_files <- list.files(tsv_dir, pattern = "\\.(rds|tsv)$", full.names = TRUE)
stopifnot(length(tsv_files) > 0)

# Cache filenames are sanitize_subtype() output, which collapses "/" and other
# punctuation to "_". Reversing that with gsub("_", " ") produces a label that
# no longer matches OncotreeSubtype ("Acute Monoblastic/Monocytic Leukemia"),
# which silently breaks both the metadata join below and the Gene -> Discover
# hand-off. Map back through the canonical lists instead.
subtype_lookup <- list()
for (.ds in c("CRISPR", "RNAi")) {
  .f <- file.path(data_dir, paste0("cancer_subtypes_", .ds, ".txt"))
  .subs <- if (file.exists(.f)) read_cancer_subtype_list_file(.f) else character(0)
  subtype_lookup[[.ds]] <- stats::setNames(.subs, sanitize_subtype(.subs))
}

dep_rows <- lapply(tsv_files, function(f) {
  base <- sub("\\.(rds|tsv)$", "", basename(f))
  ds <- if (startsWith(base, "CRISPR_")) "CRISPR" else "RNAi"
  subtype <- sub("^(CRISPR|RNAi)_", "", base)
  df <- if (grepl("\\.rds$", f)) readRDS(f) else
    suppressWarnings(readr::read_tsv(f, show_col_types = FALSE, progress = FALSE))
  if (!nrow(df) || !all(c("gene_name", "EffectSize", "p_value") %in% names(df))) {
    return(NULL)
  }
  df %>%
    dplyr::filter(!is.na(.data$EffectSize), .data$EffectSize <= EFFECT_MAX,
                  !is.na(.data$p_value), .data$p_value < PVALUE_MAX) %>%
    dplyr::transmute(
      gene_key      = cpt_gene_match_key(.data$gene_name),
      dataset       = ds,
      subtype       = unname(subtype_lookup[[ds]][subtype]) %||% gsub("_", " ", subtype),
      effect_size   = .data$EffectSize,
      p_value       = .data$p_value,
      avg           = .data$Avg,
      cancer_avg    = if ("Cancer_Avg" %in% names(df)) .data$Cancer_Avg else NA_real_,
      noncancer_avg = if ("NonCancer_Avg" %in% names(df)) .data$NonCancer_Avg else NA_real_,
      pval_vs_noncancer = if ("pval_vs_NonCancer" %in% names(df)) .data$pval_vs_NonCancer else NA_real_
    )
})
dependency <- dplyr::bind_rows(dep_rows)
msg("  dependency rows:", format(nrow(dependency), big.mark = ","))

dep_summary <- dependency %>%
  dplyr::group_by(.data$gene_key, .data$dataset) %>%
  dplyr::summarise(
    n_subtypes     = dplyr::n(),
    best_subtype   = .data$subtype[which.min(.data$effect_size)],
    best_effect    = min(.data$effect_size, na.rm = TRUE),
    median_effect  = stats::median(.data$effect_size, na.rm = TRUE),
    .groups = "drop"
  )

dependency <- dependency %>%
  dplyr::group_by(.data$gene_key, .data$dataset) %>%
  dplyr::slice_min(.data$effect_size, n = TOP_N_SUBTYPES, with_ties = FALSE) %>%
  dplyr::ungroup()
msg("  dependency rows after top-N trim:", format(nrow(dependency), big.mark = ","))

# ---- 2. Cysteine evidence + tiers ------------------------------------------
# Tiers are the PAPER's definition, computed by the repo's own
# cpt_evidence_tier(): 1 functional with atlas ligandability, 2 functional
# without it, 3 tested but not functional, 4 absent from the atlas (untested).
# Note tier 4 is a property of an *engaged* residue that the atlas never
# tested, so it cannot arise from atlas rows alone — it is computed in the
# engagement section below.
msg("Summarising cysteine atlas...")
atlas <- readRDS(file.path(data_dir, "cys_editing_atlas.rds"))
atlas_t <- atlas %>%
  dplyr::mutate(
    gene_key = cpt_gene_match_key(.data$gene_symbol),
    is_func  = cpt_is_true(.data$functional),
    is_lig   = cpt_is_true(.data$ligandable),
    tier     = cpt_evidence_tier(TRUE, .data$is_func, .data$is_lig)
  )

cys_sites_summary <- atlas_t %>%
  dplyr::group_by(.data$gene_key) %>%
  dplyr::summarise(
    n_sites               = dplyr::n(),
    n_functional          = sum(.data$is_func, na.rm = TRUE),
    n_ligandable          = sum(.data$is_lig, na.rm = TRUE),
    n_functional_ligandable = sum(.data$is_func & .data$is_lig, na.rm = TRUE),
    n_clinvar             = sum(cpt_is_true(.data$clinvar_pathogenic), na.rm = TRUE),
    best_tier             = min(.data$tier, na.rm = TRUE),  # atlas sites only (1-3)
    .groups = "drop"
  )
msg("  genes with cysteine evidence:", nrow(cys_sites_summary))

# ---- 3. Probes / ligandability ---------------------------------------------
msg("Loading protein binding table (slow, once)...")
binding <- readRDS(file.path(data_dir, "protein_binding_lookup_preprocessed.rds")) %>%
  dplyr::mutate(probe_name = simplify_probe_name(.data$probe_name))

probe_summary <- binding %>%
  dplyr::mutate(gene_key = cpt_gene_match_key(.data$gene_name)) %>%
  dplyr::group_by(.data$gene_key) %>%
  dplyr::summarise(
    n_probe_rows = dplyr::n(),
    n_probes     = dplyr::n_distinct(.data$probe_name),
    n_cys_sites  = dplyr::n_distinct(.data$cysteineid),
    max_cr       = suppressWarnings(max(.data$CR, na.rm = TRUE)),
    n_ligandable_rows = sum(cpt_is_true(.data$ligandable), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(max_cr = ifelse(is.finite(.data$max_cr), .data$max_cr, NA_real_))
msg("  genes with probe evidence:", nrow(probe_summary))

# Engaged residues (CR >= 4, the paper's chemoproteomic record) carry the
# tiers the paper actually ranks. Reuses the app's own probe-to-atlas join so
# the index cannot drift from what the Ligandability tab and MCP server report.
msg("Tiering engaged cysteines (CR >= 4)...")
engaged <- binding %>%
  dplyr::filter(!is.na(.data$CR), .data$CR >= 4) %>%
  dplyr::select(dplyr::any_of(c("gene_name", "cysteineid", "probe_name", "CR", "n_targets")))
engaged_tiered <- tryCatch(
  cpt_annotate_engaged_tiers(as.data.frame(engaged), atlas),
  error = function(e) { msg("  engaged-tier annotation failed:", conditionMessage(e)); NULL }
)
# ClinVar is a per-site annotation the layered Discover view filters on; carry
# it (and the atlas site id) alongside the tier so the app need not re-join.
if (!is.null(engaged_tiered) && nrow(engaged_tiered) && "cys_site_id" %in% names(engaged_tiered)) {
  engaged_tiered$cys_clinvar_pathogenic <- cpt_is_true(
    atlas$clinvar_pathogenic[match(engaged_tiered$cys_site_id, atlas$site_id)]
  )
}

engaged_summary <- NULL
if (!is.null(engaged_tiered) && nrow(engaged_tiered)) {
  engaged_summary <- engaged_tiered %>%
    dplyr::mutate(gene_key = cpt_gene_match_key(.data$gene_name)) %>%
    dplyr::group_by(.data$gene_key) %>%
    dplyr::summarise(
      n_engaged_sites = dplyr::n_distinct(.data$cysteineid),
      engaged_tier1 = dplyr::n_distinct(.data$cysteineid[.data$evidence_tier == 1L]),
      engaged_tier2 = dplyr::n_distinct(.data$cysteineid[.data$evidence_tier == 2L]),
      engaged_tier3 = dplyr::n_distinct(.data$cysteineid[.data$evidence_tier == 3L]),
      engaged_tier4 = dplyr::n_distinct(.data$cysteineid[.data$evidence_tier == 4L]),
      best_engaged_tier = min(.data$evidence_tier, na.rm = TRUE),
      # Sites matched on gene symbol + position because the atlas row carries
      # no UniProt accession. Position numbering is only meaningful within a
      # sequence, so these are weaker than accession matches and are counted
      # separately rather than hidden.
      n_symbol_matched = dplyr::n_distinct(
        .data$cysteineid[.data$cys_join_status == "Gene + cysteine position"]
      ),
      .groups = "drop"
    )
  msg("  genes with engaged cysteines:", nrow(engaged_summary),
      "| tier1 sites:", sum(engaged_summary$engaged_tier1))
}

probes <- binding %>%
  dplyr::mutate(gene_key = cpt_gene_match_key(.data$gene_name)) %>%
  dplyr::filter(!is.na(.data$CR)) %>%
  dplyr::group_by(.data$gene_key) %>%
  dplyr::slice_max(.data$CR, n = TOP_N_PROBES, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::select(dplyr::any_of(c("gene_key", "gene_name", "probe_name", "CR",
                                "n_targets", "cysteineid", "ligandable",
                                "Dataset", "Cell_Line", "SMILES")))
msg("  top-probe rows:", format(nrow(probes), big.mark = ","))

# ---- 4. ADME (same builder the app uses) ------------------------------------
msg("Building ADME gene scores...")
swiss <- process_swissadme_data(data_dir)
adme <- tryCatch(
  cpt_build_adme_gene_scores(cpt_index_protein_binding(binding), swiss),
  error = function(e) { msg("  ADME build failed:", conditionMessage(e)); NULL }
)
if (!is.null(adme) && nrow(adme)) {
  names(adme)[names(adme) == "gene"] <- "gene_key_src"
  msg("  ADME rows:", nrow(adme))
}

# ---- 5. Gene roster ---------------------------------------------------------
gene_display <- dplyr::bind_rows(
  atlas_t %>% dplyr::distinct(gene_key, symbol = .data$gene_symbol),
  binding %>% dplyr::transmute(gene_key = cpt_gene_match_key(.data$gene_name),
                               symbol = .data$gene_name) %>% dplyr::distinct()
) %>%
  dplyr::filter(!is.na(.data$gene_key), nzchar(.data$gene_key)) %>%
  dplyr::distinct(.data$gene_key, .keep_all = TRUE)

genes <- dplyr::tibble(gene_key = sort(unique(c(
  dependency$gene_key, cys_sites_summary$gene_key, probe_summary$gene_key
)))) %>%
  dplyr::left_join(gene_display, by = "gene_key") %>%
  dplyr::mutate(symbol = ifelse(is.na(.data$symbol), .data$gene_key, .data$symbol)) %>%
  dplyr::left_join(dplyr::select(cys_sites_summary, gene_key, n_sites, n_functional,
                                 n_ligandable, n_functional_ligandable, n_clinvar,
                                 best_atlas_site_tier = best_tier),
                   by = "gene_key") %>%
  dplyr::left_join(dplyr::select(probe_summary, gene_key, n_probes, max_cr), by = "gene_key") %>%
  { if (is.null(engaged_summary)) . else dplyr::left_join(., engaged_summary, by = "gene_key") } %>%
  dplyr::left_join(
    dep_summary %>%
      dplyr::group_by(.data$gene_key) %>%
      dplyr::summarise(
        dep_n_subtypes = sum(.data$n_subtypes),
        dep_best_effect = min(.data$best_effect, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "gene_key"
  )
msg("  genes in roster:", format(nrow(genes), big.mark = ","))

# ---- 6. Effect-size distributions for box plots -----------------------------
# The group-comparison box plot is the one Gene-tab visual that would otherwise
# need a 140 MB matrix at runtime. Five-number summaries per gene x subtype are
# all a box plot consumes, so they are computed here instead. Groups match
# cpt_group_comparison_df(): the subtype, all other cancer lines, and the
# non-cancerous lines (the last is constant per gene, so stored once).
msg("Computing effect-size distributions for box plots...")
DIST_PROBS <- c(0, 0.25, 0.5, 0.75, 1)
meta <- readRDS(file.path(data_dir, "cancer_model_data.rds"))
dist_parts <- list()
point_parts <- list()

summarise_cols <- function(mat, rows, cols, gene_keys, label, dataset, subtype) {
  if (!length(rows) || !length(cols)) return(NULL)
  sub <- mat[rows, cols, drop = FALSE]
  q <- matrixStats::colQuantiles(sub, probs = DIST_PROBS, na.rm = TRUE)
  n <- matrixStats::colCounts(!is.na(sub), value = TRUE)
  # mean and sd travel with the quantiles so a Welch t-test between two
  # subtypes can be computed exactly from summaries, with no matrix at runtime.
  data.frame(
    gene_key = gene_keys, dataset = dataset, subtype = subtype, group = label,
    n = as.integer(n),
    ymin = q[, 1], lower = q[, 2], middle = q[, 3], upper = q[, 4], ymax = q[, 5],
    mean = matrixStats::colMeans2(sub, na.rm = TRUE),
    sd   = matrixStats::colSds(sub, na.rm = TRUE),
    stringsAsFactors = FALSE, row.names = NULL
  )
}

for (ds in c("CRISPR", "RNAi")) {
  mpath <- file.path(data_dir, if (identical(ds, "CRISPR")) {
    "CRISPRGeneEffect_23Q4_clean.rds"
  } else {
    "d2_gene_effect_headers_refined.rds"
  })
  if (!file.exists(mpath)) { msg("  skipping", ds, "- matrix not present"); next }
  mat <- readRDS(mpath)
  # The gene-effect objects are data.frames on disk; matrixStats needs a matrix.
  if (!is.matrix(mat)) mat <- as.matrix(mat)
  meta_a <- meta[match(rownames(mat), meta$ModelID), , drop = FALSE]
  col_key <- toupper(colnames(mat))
  noncancer_idx <- which(meta_a$OncotreePrimaryDisease == "Non-Cancerous")

  dep_ds <- dependency[dependency$dataset == ds, , drop = FALSE]
  if (!nrow(dep_ds)) next

  genes_ds <- unique(dep_ds$gene_key)
  gcols <- match(genes_ds, col_key)
  keep <- !is.na(gcols)
  dist_parts[[length(dist_parts) + 1L]] <- summarise_cols(
    mat, noncancer_idx, gcols[keep], genes_ds[keep], "Non-Cancer", ds, NA_character_
  )
  # Whole-matrix distribution: the same population the "average gene effect
  # above -0.5" common-essential filter is computed over, so a gene's subtype
  # boxes can be read against every screened line.
  dist_parts[[length(dist_parts) + 1L]] <- summarise_cols(
    mat, seq_len(nrow(mat)), gcols[keep], genes_ds[keep], "All lines", ds, NA_character_
  )

  for (st in unique(dep_ds$subtype)) {
    idx_s <- which(meta_a$OncotreeSubtype == st)
    idx_o <- which(meta_a$OncotreeSubtype != st &
                     meta_a$OncotreePrimaryDisease != "Non-Cancerous")
    gk <- unique(dep_ds$gene_key[dep_ds$subtype == st])
    gc <- match(gk, col_key)
    ok <- !is.na(gc)
    if (!any(ok)) next
    dist_parts[[length(dist_parts) + 1L]] <- summarise_cols(
      mat, idx_s, gc[ok], gk[ok], "Cancer of Interest", ds, st)
    dist_parts[[length(dist_parts) + 1L]] <- summarise_cols(
      mat, idx_o, gc[ok], gk[ok], "Other Cancers", ds, st)

    # Individual cell lines behind the subtype box, so the plot can show the
    # points and name them. Stored for the subtype group only: the reference
    # groups run to ~1,100 lines per gene, which is neither plottable as jitter
    # nor worth the size.
    if (length(idx_s)) {
      vals <- mat[idx_s, gc[ok], drop = FALSE]
      line_names <- if ("StrippedCellLineName" %in% names(meta_a)) {
        as.character(meta_a$StrippedCellLineName[idx_s])
      } else {
        rownames(mat)[idx_s]
      }
      line_names[is.na(line_names) | !nzchar(line_names)] <- rownames(mat)[idx_s][is.na(line_names) | !nzchar(line_names)]
      keep_pt <- !is.na(vals)
      if (any(keep_pt)) {
        wi <- which(keep_pt, arr.ind = TRUE)
        point_parts[[length(point_parts) + 1L]] <- data.frame(
          gene_key = gk[ok][wi[, "col"]],
          dataset = ds,
          subtype = st,
          cell_line = line_names[wi[, "row"]],
          value = vals[keep_pt],
          stringsAsFactors = FALSE, row.names = NULL
        )
      }
    }
  }
  rm(mat); invisible(gc())
  msg("  ", ds, "distributions done")
}
effect_distributions <- dplyr::bind_rows(dist_parts)
msg("  distribution rows:", format(nrow(effect_distributions), big.mark = ","))
effect_points <- dplyr::bind_rows(point_parts)
msg("  per-cell-line points:", format(nrow(effect_points), big.mark = ","))

out <- list(
  built_at          = Sys.time(),
  effect_distributions = effect_distributions,
  effect_points     = effect_points,
  effect_max        = EFFECT_MAX,
  pvalue_max        = PVALUE_MAX,
  genes             = genes,
  dependency        = dependency,
  dep_summary       = dep_summary,
  cys_sites_summary = cys_sites_summary,
  engaged_summary   = engaged_summary,
  engaged_sites     = if (is.null(engaged_tiered)) NULL else
                        engaged_tiered[, intersect(names(engaged_tiered),
                          c("gene_name", "cysteineid", "probe_name", "CR", "n_targets",
                            "cys_site_in_atlas", "cys_functional", "cys_ligandable",
                            "evidence_tier", "evidence_tier_label",
                            "cys_join_status", "cys_residue_mapping_status",
                            "cys_site_id", "cys_clinvar_pathogenic",
                            "cys_study_context", "cys_editor_support")), drop = FALSE],
  probes            = probes,
  probe_summary     = probe_summary,
  adme              = adme
)

out_path <- file.path(data_dir, "gene_index.rds")
# gzip, not xz. The index is read at session start and is small enough that
# decompression is the whole cost: 2.2 s as xz against 1.0 s as gzip, for four
# megabytes more on disk.
saveRDS(out, out_path, compress = "gzip")
msg("Wrote", out_path, sprintf("(%.1f MB)", file.size(out_path) / 1024^2))
