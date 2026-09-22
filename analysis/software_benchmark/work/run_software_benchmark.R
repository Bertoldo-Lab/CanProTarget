#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
benchmark_dir <- if (length(args)) {
  normalizePath(args[[1]], mustWork = TRUE)
} else {
  normalizePath(file.path(getwd(), "manuscript_revision_2026-08-03", "software_benchmark"),
                mustWork = TRUE)
}
project_root <- normalizePath(file.path(benchmark_dir, "..", ".."), mustWork = TRUE)
output_dir <- file.path(benchmark_dir, "outputs")
figure_dir <- file.path(benchmark_dir, "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# Provenance records paths relative to project_root so the manifest is
# portable and does not carry the operator's home directory into the repo.
rel_path <- function(p) {
  p <- normalizePath(p, mustWork = FALSE)
  if (startsWith(p, project_root)) sub("^/", "", substring(p, nchar(project_root) + 1L)) else basename(p)
}

# The application moved from the development tree (CanProTarget_dev/, with the
# manuscript inputs beside it) into the release repository, where the app sits
# at the root and the 2026-07-13 dependency snapshot is outside the repo. Both
# layouts are resolved here so the run is reproducible from either. The second
# argument, or CPT_BENCHMARK_INPUTS, points at the directory holding
# rna_all_gene_statistics.rds.
first_existing <- function(...) {
  cand <- c(...)
  hit <- cand[file.exists(cand)]
  if (length(hit)) hit[[1]] else cand[[1]]
}
app_root <- if (dir.exists(file.path(project_root, "CanProTarget_dev", "R"))) {
  file.path(project_root, "CanProTarget_dev")
} else {
  project_root
}
score_file   <- file.path(app_root, "R", "canprotarget_score.R")
atlas_file   <- file.path(app_root, "data", "cys_editing_atlas.rds")
# ADME is built from the chemoproteomics binding table and the SwissADME
# descriptors, exactly as the application builds it at startup.
binding_file <- first_existing(
  file.path(app_root, "data", "protein_binding_lookup_factored.rds"),
  file.path(app_root, "data", "protein_binding_lookup_preprocessed.rds")
)
adme_file    <- file.path(app_root, "data", "swissadme_preprocessed.rds")

archived_benchmark_dir <- if (length(args) >= 2) {
  normalizePath(args[[2]], mustWork = TRUE)
} else if (nzchar(Sys.getenv("CPT_BENCHMARK_INPUTS"))) {
  normalizePath(Sys.getenv("CPT_BENCHMARK_INPUTS"), mustWork = TRUE)
} else {
  file.path(project_root, "Original_CanProTarget",
            "manuscript_revision_2026-07-13", "data", "benchmark")
}
all_gene_stats_file <- file.path(archived_benchmark_dir, "rna_all_gene_statistics.rds")
site_summary_file <- file.path(archived_benchmark_dir, "gene_level_vs_exact_summary.csv")
site_cohort_file <- file.path(archived_benchmark_dir, "cohort_workflow_counts.csv")
covpdb_file <- file.path(benchmark_dir, "data_external", "covpdb_cysteine_targets_human.csv")
uniprot_file <- file.path(benchmark_dir, "data_external", "uniprot_covpdb_gene_mapping.tsv")

required_files <- c(score_file, atlas_file, binding_file, adme_file, all_gene_stats_file,
                    site_summary_file, site_cohort_file, covpdb_file, uniprot_file)
if (!all(file.exists(required_files))) {
  stop("Missing benchmark input(s): ",
       paste(required_files[!file.exists(required_files)], collapse = ", "), call. = FALSE)
}

source(score_file, local = globalenv())

clean_gene <- function(x) toupper(sub(" \\(\\d+\\)$", "", trimws(as.character(x))))
finite_or_na <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x[!is.finite(x)] <- NA_real_
  x
}
write_output <- function(x, name) {
  write.csv(x, file.path(output_dir, name), row.names = FALSE, na = "")
}

# ---- Inputs and independent positive set -----------------------

all_gene_stats <- readRDS(all_gene_stats_file)
atlas <- readRDS(atlas_file)

# The 2026-08-03 run predated the ADME dimension by one day and pinned
# adme_druggability to NA, so every ADME weight it recorded was inert. The
# dimension is now computed the way the application computes it: mean fragment
# developability over the covalent probes engaging the gene at CR >= 4.
adme_gene_scores <- local({
  binding <- readRDS(binding_file)
  descriptors <- readRDS(adme_file)
  out <- cpt_build_adme_gene_scores(binding, descriptors, min_cr = 4)
  rm(binding, descriptors); gc(verbose = FALSE)
  out
})
cat(sprintf("Built ADME gene scores for %d genes.\n", nrow(adme_gene_scores)))
covpdb <- read.csv(covpdb_file, stringsAsFactors = FALSE, check.names = FALSE)
uniprot <- read.delim(uniprot_file, stringsAsFactors = FALSE, check.names = FALSE)
names(uniprot)[names(uniprot) == "Entry"] <- "uniprot_accession"
names(uniprot)[names(uniprot) == "Gene Names (primary)"] <- "gene_symbol"

covpdb_map <- merge(covpdb, uniprot[c("uniprot_accession", "gene_symbol")],
                    by = "uniprot_accession", all.x = TRUE, sort = FALSE)
covpdb_map$gene_symbol <- trimws(covpdb_map$gene_symbol)
covpdb_map$gene_key <- clean_gene(covpdb_map$gene_symbol)
covpdb_map$mapping_status <- ifelse(nzchar(covpdb_map$gene_key), "mapped", "unmapped")
covpdb_map <- covpdb_map[order(covpdb_map$covpdb_index), ]
write_output(covpdb_map, "covpdb_human_target_gene_mapping.csv")

known_positive_genes <- unique(covpdb_map$gene_key[covpdb_map$mapping_status == "mapped"])
if (length(known_positive_genes) < 80L) {
  stop("Unexpectedly low CovPDB-to-gene mapping coverage.", call. = FALSE)
}

# ---- Vectorized cysteine dimensions faithful to current code ---

atlas$gene_key <- clean_gene(atlas$gene_symbol)
atlas_groups <- split(seq_len(nrow(atlas)), atlas$gene_key, drop = TRUE)

summarise_atlas_gene <- function(idx) {
  x <- atlas[idx, , drop = FALSE]

  ligandability <- NA_real_
  if ("ligandability_score" %in% names(x) && any(!is.na(x$ligandability_score))) {
    ligandability <- max(finite_or_na(x$ligandability_score), na.rm = TRUE)
  } else if ("ligandable" %in% names(x)) {
    vals <- x$ligandable
    if (!all(is.na(vals))) ligandability <- if (any(vals, na.rm = TRUE)) 50 else 0
  } else if ("functional_ligandable" %in% names(x)) {
    vals <- x$functional_ligandable
    if (!all(is.na(vals))) ligandability <- if (any(vals, na.rm = TRUE)) 75 else 0
  }

  conservation <- NA_real_
  if ("conservation_score" %in% names(x)) {
    vals <- finite_or_na(x$conservation_score)
    if (any(!is.na(vals))) conservation <- max(vals, na.rm = TRUE) * 100
  }

  clinical <- NA_real_
  if ("clinvar_pathogenic" %in% names(x)) {
    vals <- x$clinvar_pathogenic
    if (!all(is.na(vals))) clinical <- if (any(vals, na.rm = TRUE)) 100 else 0
  }

  c(ligandability = ligandability, conservation = conservation, clinical = clinical)
}

atlas_gene_matrix <- t(vapply(atlas_groups, summarise_atlas_gene, numeric(3)))
atlas_gene_dims <- data.frame(
  gene_key = rownames(atlas_gene_matrix),
  dim_ligandability = atlas_gene_matrix[, "ligandability"],
  dim_conservation = atlas_gene_matrix[, "conservation"],
  dim_clinical = atlas_gene_matrix[, "clinical"],
  stringsAsFactors = FALSE,
  row.names = NULL
)

compute_composite <- function(dim_df, weights, min_dimensions = 2L) {
  canonical <- c("dependency_strength", "cancer_selectivity", "cysteine_ligandability",
                 "conservation", "clinical_evidence", "adme_druggability")
  values <- as.matrix(dim_df[, canonical, drop = FALSE])
  storage.mode(values) <- "double"
  w <- unlist(weights)[canonical]
  keep <- is.finite(w) & w > 0
  values <- values[, keep, drop = FALSE]
  w <- w[keep]
  available <- is.finite(values)
  safe_values <- values
  safe_values[!available] <- 0
  numerator <- rowSums(sweep(safe_values, 2, w, `*`))
  denominator <- rowSums(sweep(available, 2, w, `*`))
  n_used <- rowSums(available)
  score <- rep(NA_real_, nrow(dim_df))
  valid <- n_used >= min_dimensions & denominator > 0
  score[valid] <- numerator[valid] / denominator[valid]
  score
}

active_dimensions <- c("dependency_strength", "cancer_selectivity",
                       "cysteine_ligandability", "conservation", "clinical_evidence")
all_dimensions <- c(active_dimensions, "adme_druggability")

# The reference configuration is the shipped default, which weights ADME 0.
# Two scenarios exist to answer questions the manuscript asks and the previous
# run could not: what activating ADME at the proposed 0.5 does, and what equal
# weighting means when it includes the sixth dimension rather than excluding it.
equal_five <- as.list(setNames(c(rep(1, length(active_dimensions)), 0), all_dimensions))
equal_six  <- as.list(setNames(rep(1, length(all_dimensions)), all_dimensions))
adme_active <- CPT_DEFAULT_WEIGHTS
adme_active$adme_druggability <- CPT_ADME_PROPOSED_WEIGHT

scenarios <- list(default = CPT_DEFAULT_WEIGHTS,
                  adme_active = adme_active,
                  equal_six = equal_six,
                  equal_five = equal_five)
for (dimension in active_dimensions) {
  w_leave <- CPT_DEFAULT_WEIGHTS
  w_leave[[dimension]] <- 0
  scenarios[[paste0("leave_out_", dimension)]] <- w_leave
  for (multiplier in c(0.5, 2)) {
    w_perturb <- CPT_DEFAULT_WEIGHTS
    w_perturb[[dimension]] <- w_perturb[[dimension]] * multiplier
    suffix <- if (multiplier == 0.5) "half" else "double"
    scenarios[[paste0(suffix, "_", dimension)]] <- w_perturb
  }
}

weight_table <- do.call(rbind, lapply(names(scenarios), function(nm) {
  data.frame(scenario = nm, as.list(unlist(scenarios[[nm]])), check.names = FALSE)
}))
write_output(weight_table, "cpt_weight_scenarios.csv")

# ---- Score every cohort in the manuscript dependency universe --

cohort_keys <- unique(all_gene_stats[c("panel", "cohort_code", "cohort")])
cohort_keys <- cohort_keys[order(cohort_keys$panel, cohort_keys$cohort_code), ]
baseline_parts <- vector("list", nrow(cohort_keys))

for (i in seq_len(nrow(cohort_keys))) {
  key <- cohort_keys[i, ]
  x <- all_gene_stats[
    all_gene_stats$panel == key$panel & all_gene_stats$cohort_code == key$cohort_code,
    , drop = FALSE
  ]
  x$gene_key <- clean_gene(x$gene_name)
  x$dependency_strength <- score_dependency_strength(finite_or_na(x$Cancer_Avg))
  x$cancer_selectivity <- score_cancer_selectivity(finite_or_na(x$EffectSize))
  atlas_match <- match(x$gene_key, atlas_gene_dims$gene_key)
  x$cysteine_ligandability <- atlas_gene_dims$dim_ligandability[atlas_match]
  x$conservation <- atlas_gene_dims$dim_conservation[atlas_match]
  x$clinical_evidence <- atlas_gene_dims$dim_clinical[atlas_match]
  x$adme_druggability <- score_adme_druggability(x$gene_key, adme_gene_scores)

  dim_df <- x[c("dependency_strength", "cancer_selectivity", "cysteine_ligandability",
                "conservation", "clinical_evidence", "adme_druggability")]
  for (scenario in names(scenarios)) {
    x[[paste0("score_", scenario)]] <- compute_composite(dim_df, scenarios[[scenario]])
  }

  baseline <- !is.na(x$EffectSize) & x$EffectSize < -0.1 &
    !is.na(x$p_value) & x$p_value < 0.05 &
    !is.na(x$Avg) & x$Avg > -0.5
  x <- x[baseline, , drop = FALSE]
  x$known_covpdb_target <- x$gene_key %in% known_positive_genes
  baseline_parts[[i]] <- x
  cat(sprintf("Scored %s/%s: %d baseline dependencies.\n",
              key$panel, key$cohort_code, nrow(x)))
}

baseline_scores <- do.call(rbind, baseline_parts)
rownames(baseline_scores) <- NULL

expected_counts <- c(Adult = 5570L, Paediatric = 9307L)
observed_counts <- table(baseline_scores$panel)
if (!identical(as.integer(observed_counts[names(expected_counts)]), as.integer(expected_counts))) {
  stop("Dependency universe did not reproduce the validated adult/paediatric counts.", call. = FALSE)
}

score_columns <- paste0("score_", names(scenarios))
export_columns <- c(
  "panel", "cohort_code", "cohort", "gene_name", "gene_key", "EffectSize", "p_value",
  "Avg", "Cancer_Avg", "dependency_strength", "cancer_selectivity",
  "cysteine_ligandability", "conservation", "clinical_evidence", "adme_druggability",
  "known_covpdb_target", "score_default", "score_adme_active", "score_equal_six",
  "score_equal_five"
)
write_output(baseline_scores[export_columns], "cpt_baseline_dependency_scores.csv")

# Independently call the application annotation function to confirm that the
# vectorized benchmark implementation reproduces every displayed CPT score.
validation_rows <- lapply(seq_len(nrow(cohort_keys)), function(i) {
  key <- cohort_keys[i, ]
  background <- all_gene_stats[
    all_gene_stats$panel == key$panel & all_gene_stats$cohort_code == key$cohort_code,
    , drop = FALSE
  ]
  gene_df <- background[
    !is.na(background$EffectSize) & background$EffectSize < -0.1 &
      !is.na(background$p_value) & background$p_value < 0.05 &
      !is.na(background$Avg) & background$Avg > -0.5,
    , drop = FALSE
  ]
  direct <- cpt_annotate_gene_table(gene_df, background, atlas,
                                    adme_data = adme_gene_scores)
  saved <- baseline_scores[
    baseline_scores$panel == key$panel & baseline_scores$cohort_code == key$cohort_code,
    c("gene_name", "score_default"), drop = FALSE
  ]
  comparison <- merge(direct[c("gene_name", "CPT_Score")], saved, by = "gene_name")
  equal_after_display_rounding <- round(comparison$score_default, 2) == comparison$CPT_Score
  data.frame(
    panel = key$panel,
    cohort_code = key$cohort_code,
    n_compared = nrow(comparison),
    n_equal_after_display_rounding = sum(equal_after_display_rounding, na.rm = TRUE),
    maximum_absolute_unrounded_difference = max(
      abs(comparison$CPT_Score - comparison$score_default), na.rm = TRUE
    ),
    stringsAsFactors = FALSE
  )
})
implementation_validation <- do.call(rbind, validation_rows)
write_output(implementation_validation, "cpt_implementation_validation.csv")
if (any(implementation_validation$n_compared !=
        implementation_validation$n_equal_after_display_rounding)) {
  stop("Vectorized benchmark scores do not reproduce displayed application scores.",
       call. = FALSE)
}

# ---- Weight robustness ----------------------------------------

top_overlap <- function(default_score, scenario_score, ids, k) {
  if (!length(ids)) return(c(overlap = NA_real_, jaccard = NA_real_, n_top = 0))
  d <- default_score
  s <- scenario_score
  d[!is.finite(d)] <- -Inf
  s[!is.finite(s)] <- -Inf
  k <- min(as.integer(k), length(ids))
  top_d <- ids[order(-d, ids)][seq_len(k)]
  top_s <- ids[order(-s, ids)][seq_len(k)]
  intersection <- length(intersect(top_d, top_s))
  union <- length(union(top_d, top_s))
  c(overlap = intersection / k, jaccard = intersection / union, n_top = k)
}

robustness_rows <- list()
rank_shift_rows <- list()
row_index <- 0L
shift_index <- 0L

for (i in seq_len(nrow(cohort_keys))) {
  key <- cohort_keys[i, ]
  x <- baseline_scores[
    baseline_scores$panel == key$panel & baseline_scores$cohort_code == key$cohort_code,
    , drop = FALSE
  ]
  default_score <- x$score_default
  default_rank <- rank(-default_score, ties.method = "average", na.last = "keep")
  top_default <- is.finite(default_rank) & default_rank <= min(50L, sum(is.finite(default_rank)))

  for (scenario in setdiff(names(scenarios), "default")) {
    scenario_score <- x[[paste0("score_", scenario)]]
    valid <- is.finite(default_score) & is.finite(scenario_score)
    spearman <- if (sum(valid) >= 3L) {
      suppressWarnings(cor(default_score[valid], scenario_score[valid], method = "spearman"))
    } else NA_real_

    fixed_metrics <- lapply(c(10L, 25L, 50L), function(k) {
      overlap <- top_overlap(default_score, scenario_score, x$gene_key, k)
      data.frame(top_definition = paste0("top_", k), top_n = overlap[["n_top"]],
                 overlap_fraction = overlap[["overlap"]], jaccard = overlap[["jaccard"]])
    })
    fractional_k <- max(1L, ceiling(sum(valid) * 0.10))
    overlap <- top_overlap(default_score, scenario_score, x$gene_key, fractional_k)
    fixed_metrics[[length(fixed_metrics) + 1L]] <- data.frame(
      top_definition = "top_10_percent", top_n = overlap[["n_top"]],
      overlap_fraction = overlap[["overlap"]], jaccard = overlap[["jaccard"]]
    )
    metrics <- do.call(rbind, fixed_metrics)
    metrics$panel <- key$panel
    metrics$cohort_code <- key$cohort_code
    metrics$cohort <- key$cohort
    metrics$scenario <- scenario
    metrics$n_candidates <- nrow(x)
    metrics$n_compared <- sum(valid)
    metrics$spearman <- spearman
    row_index <- row_index + 1L
    robustness_rows[[row_index]] <- metrics[c(
      "panel", "cohort_code", "cohort", "scenario", "n_candidates", "n_compared", "spearman",
      "top_definition", "top_n", "overlap_fraction", "jaccard"
    )]

    scenario_rank <- rank(-scenario_score, ties.method = "average", na.last = "keep")
    shifts <- abs(scenario_rank[top_default] - default_rank[top_default])
    finite_shifts <- shifts[is.finite(shifts)]
    shift_index <- shift_index + 1L
    rank_shift_rows[[shift_index]] <- data.frame(
      panel = key$panel,
      cohort_code = key$cohort_code,
      scenario = scenario,
      n_default_top50 = sum(top_default),
      n_default_top50_scored_under_scenario = length(finite_shifts),
      median_absolute_rank_shift = if (length(finite_shifts)) median(finite_shifts) else NA_real_,
      p95_absolute_rank_shift = if (length(finite_shifts)) {
        as.numeric(quantile(finite_shifts, 0.95, names = FALSE))
      } else NA_real_,
      maximum_absolute_rank_shift = if (length(finite_shifts)) max(finite_shifts) else NA_real_,
      stringsAsFactors = FALSE
    )
  }
}

robustness <- do.call(rbind, robustness_rows)
rank_shifts <- do.call(rbind, rank_shift_rows)
write_output(robustness, "cpt_weight_robustness_by_cohort.csv")
write_output(rank_shifts, "cpt_top50_rank_shift_by_cohort.csv")

robustness_summary <- do.call(rbind, lapply(split(robustness, robustness$scenario), function(x) {
  x25 <- x[x$top_definition == "top_25", ]
  x10p <- x[x$top_definition == "top_10_percent", ]
  data.frame(
    scenario = unique(x$scenario),
    n_cohorts = length(unique(x$cohort_code)),
    median_score_coverage = median(x$n_compared / x$n_candidates, na.rm = TRUE),
    minimum_score_coverage = min(x$n_compared / x$n_candidates, na.rm = TRUE),
    median_spearman = median(x$spearman, na.rm = TRUE),
    minimum_spearman = min(x$spearman, na.rm = TRUE),
    median_top25_overlap = median(x25$overlap_fraction, na.rm = TRUE),
    minimum_top25_overlap = min(x25$overlap_fraction, na.rm = TRUE),
    median_top10pct_overlap = median(x10p$overlap_fraction, na.rm = TRUE),
    minimum_top10pct_overlap = min(x10p$overlap_fraction, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
robustness_summary <- robustness_summary[order(-robustness_summary$median_spearman), ]
write_output(robustness_summary, "cpt_weight_robustness_summary.csv")

# ---- Independent CovPDB positive-set enrichment ---------------

baseline_scores$ranker_CPT_default <- baseline_scores$score_default
baseline_scores$ranker_CPT_equal <- baseline_scores$score_equal_five
baseline_scores$ranker_CPT_adme_active <- baseline_scores$score_adme_active
baseline_scores$ranker_Dependency <- baseline_scores$dependency_strength
baseline_scores$ranker_Selectivity <- baseline_scores$cancer_selectivity
baseline_scores$ranker_Dependency_selectivity_mean <- rowMeans(
  baseline_scores[c("dependency_strength", "cancer_selectivity")], na.rm = TRUE
)
baseline_scores$ranker_Cysteine_ligandability <- baseline_scores$cysteine_ligandability

rankers <- c(
  CPT_default = "ranker_CPT_default",
  CPT_adme_active = "ranker_CPT_adme_active",
  CPT_equal = "ranker_CPT_equal",
  Dependency = "ranker_Dependency",
  Selectivity = "ranker_Selectivity",
  Dependency_selectivity_mean = "ranker_Dependency_selectivity_mean",
  Cysteine_ligandability = "ranker_Cysteine_ligandability"
)

enrichment_rows <- list()
enrichment_index <- 0L
fractions <- c(0.05, 0.10, 0.20)

for (ranker in names(rankers)) {
  score_col <- rankers[[ranker]]
  for (fraction in fractions) {
    for (i in seq_len(nrow(cohort_keys))) {
      key <- cohort_keys[i, ]
      x <- baseline_scores[
        baseline_scores$panel == key$panel & baseline_scores$cohort_code == key$cohort_code,
        , drop = FALSE
      ]
      score <- x[[score_col]]
      ids <- x$gene_key
      n_top <- max(1L, ceiling(nrow(x) * fraction))
      ordering_score <- score
      ordering_score[!is.finite(ordering_score)] <- -Inf
      top_idx <- order(-ordering_score, ids)[seq_len(n_top)]
      n_positive <- sum(x$known_covpdb_target)
      top_positive <- sum(x$known_covpdb_target[top_idx])
      enrichment_index <- enrichment_index + 1L
      enrichment_rows[[enrichment_index]] <- data.frame(
        panel = key$panel,
        cohort_code = key$cohort_code,
        cohort = key$cohort,
        ranker = ranker,
        top_fraction = fraction,
        n_candidates = nrow(x),
        n_scored = sum(is.finite(score)),
        n_known_positive = n_positive,
        n_known_positive_scored = sum(x$known_covpdb_target & is.finite(score)),
        n_top = n_top,
        n_top_known_positive = top_positive,
        stringsAsFactors = FALSE
      )
    }
  }
}

enrichment_by_cohort <- do.call(rbind, enrichment_rows)
enrichment_by_cohort$overall_prevalence <- with(
  enrichment_by_cohort, n_known_positive / n_candidates
)
enrichment_by_cohort$top_prevalence <- with(
  enrichment_by_cohort, n_top_known_positive / n_top
)
enrichment_by_cohort$fold_enrichment <- with(
  enrichment_by_cohort, top_prevalence / overall_prevalence
)
enrichment_by_cohort$known_positive_recall <- with(
  enrichment_by_cohort, ifelse(n_known_positive > 0, n_top_known_positive / n_known_positive, NA_real_)
)
write_output(enrichment_by_cohort, "covpdb_positive_set_enrichment_by_cohort.csv")

set.seed(20260803)
bootstrap_enrichment <- function(x, iterations = 2000L) {
  cohort_ids <- unique(x$cohort_code)
  ratios <- replicate(iterations, {
    sampled <- sample(cohort_ids, length(cohort_ids), replace = TRUE)
    sampled_rows <- do.call(rbind, lapply(sampled, function(id) x[x$cohort_code == id, ]))
    top_prevalence <- sum(sampled_rows$n_top_known_positive) / sum(sampled_rows$n_top)
    overall_prevalence <- sum(sampled_rows$n_known_positive) / sum(sampled_rows$n_candidates)
    top_prevalence / overall_prevalence
  })
  as.numeric(quantile(ratios, c(0.025, 0.975), na.rm = TRUE, names = FALSE))
}

enrichment_summary <- do.call(rbind, lapply(
  split(enrichment_by_cohort,
        interaction(enrichment_by_cohort$ranker, enrichment_by_cohort$top_fraction, drop = TRUE)),
  function(x) {
    ci <- bootstrap_enrichment(x)
    n_candidates <- sum(x$n_candidates)
    n_positive <- sum(x$n_known_positive)
    n_top <- sum(x$n_top)
    n_top_positive <- sum(x$n_top_known_positive)
    overall_prevalence <- n_positive / n_candidates
    top_prevalence <- n_top_positive / n_top
    data.frame(
      ranker = unique(x$ranker),
      top_fraction = unique(x$top_fraction),
      n_cohorts = nrow(x),
      n_candidate_associations = n_candidates,
      n_scored_associations = sum(x$n_scored),
      n_known_positive_associations = n_positive,
      n_known_positive_scored = sum(x$n_known_positive_scored),
      n_top_associations = n_top,
      n_top_known_positive = n_top_positive,
      overall_positive_prevalence = overall_prevalence,
      top_positive_prevalence = top_prevalence,
      fold_enrichment = top_prevalence / overall_prevalence,
      bootstrap_ci_low = ci[[1]],
      bootstrap_ci_high = ci[[2]],
      known_positive_recall = n_top_positive / n_positive,
      score_coverage = sum(x$n_scored) / n_candidates,
      known_positive_score_coverage = sum(x$n_known_positive_scored) / n_positive,
      stringsAsFactors = FALSE
    )
  }
))
enrichment_summary <- enrichment_summary[order(enrichment_summary$top_fraction,
                                               -enrichment_summary$fold_enrichment), ]
write_output(enrichment_summary, "covpdb_positive_set_enrichment_summary.csv")

# ---- Exact-site evidence-resolution benchmark -----------------

site_summary <- read.csv(site_summary_file, stringsAsFactors = FALSE, check.names = FALSE)
site_summary$exact_site_retention_fraction <- with(
  site_summary, exact_site_functional_gene_associations / gene_level_combined_associations
)
site_summary$gene_level_precision_if_exact_site_is_reference <- site_summary$exact_site_retention_fraction
write_output(site_summary, "exact_site_discrimination_summary.csv")

site_cohort <- read.csv(site_cohort_file, stringsAsFactors = FALSE, check.names = FALSE)
site_cohort$exact_site_retention_fraction <- with(
  site_cohort, ifelse(gene_level_combined_associations > 0,
                      exact_site_functional_gene_associations / gene_level_combined_associations,
                      NA_real_)
)
write_output(site_cohort, "exact_site_discrimination_by_cohort.csv")

# ---- Coverage and provenance ----------------------------------

dimension_coverage <- data.frame(
  dimension = c("dependency_strength", "cancer_selectivity", "cysteine_ligandability",
                "conservation", "clinical_evidence", "adme_druggability", "CPT_default"),
  n_available = c(
    sum(is.finite(baseline_scores$dependency_strength)),
    sum(is.finite(baseline_scores$cancer_selectivity)),
    sum(is.finite(baseline_scores$cysteine_ligandability)),
    sum(is.finite(baseline_scores$conservation)),
    sum(is.finite(baseline_scores$clinical_evidence)),
    sum(is.finite(baseline_scores$adme_druggability)),
    sum(is.finite(baseline_scores$score_default))
  ),
  n_total = nrow(baseline_scores),
  stringsAsFactors = FALSE
)
dimension_coverage$fraction_available <- dimension_coverage$n_available / dimension_coverage$n_total
write_output(dimension_coverage, "cpt_dimension_coverage.csv")

dimension_count <- rowSums(is.finite(as.matrix(baseline_scores[c(
  "dependency_strength", "cancer_selectivity", "cysteine_ligandability",
  "conservation", "clinical_evidence", "adme_druggability"
)])))
dimension_count_distribution <- as.data.frame(table(dimension_count), stringsAsFactors = FALSE)
names(dimension_count_distribution) <- c("n_available_dimensions", "n_associations")
dimension_count_distribution$n_available_dimensions <- as.integer(
  as.character(dimension_count_distribution$n_available_dimensions)
)
dimension_count_distribution$n_associations <- as.integer(
  dimension_count_distribution$n_associations
)
dimension_count_distribution$fraction_associations <-
  dimension_count_distribution$n_associations / nrow(baseline_scores)
write_output(dimension_count_distribution, "cpt_dimension_count_distribution.csv")

covpdb_coverage <- data.frame(
  measure = c("CovPDB human cysteine targets", "Mapped to primary gene symbol",
              "Present in RNAi background", "Observed among baseline dependency associations"),
  n = c(
    nrow(covpdb_map),
    sum(covpdb_map$mapping_status == "mapped"),
    sum(known_positive_genes %in% unique(clean_gene(all_gene_stats$gene_name))),
    length(intersect(known_positive_genes,
                     unique(baseline_scores$gene_key[baseline_scores$known_covpdb_target])))
  ),
  stringsAsFactors = FALSE
)
write_output(covpdb_coverage, "covpdb_reference_coverage.csv")

resource_comparison <- data.frame(
  resource = c("CanProTarget", "Open Targets Platform", "CysDB", "DrugMap",
               "CovPDB", "CovalentInDB 2.0"),
  primary_scope = c(
    "Cancer dependency-to-cysteine prioritisation",
    "Target-disease evidence integration",
    "Aggregated human cysteine chemoproteomics",
    "Pan-cancer, cell-line-resolved cysteine ligandability",
    "Experimentally resolved covalent protein-ligand structures",
    "Covalent inhibitors, targets and structural annotations"
  ),
  cohort_specific_dependency = c("Yes", "Disease-level rather than this cohort workflow",
                                  "No", "No", "No", "No"),
  exact_residue_evidence = c("Yes", "Limited", "Yes", "Yes", "Yes", "Yes"),
  residue_function_evidence = c("Yes, through Cys_editing", "No", "Annotations vary by source",
                                "No direct editing evidence", "No", "No"),
  configurable_integrated_rank = c("Yes", "Association score", "No", "No", "No", "No"),
  comparison_type = c("Index workflow", rep("Complementary resource", 5)),
  source = c(
    "Local implementation benchmark",
    "https://platform-docs.opentargets.org/associations",
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC10510411/",
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC11143475/",
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC8728183/",
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC11701572/"
  ),
  stringsAsFactors = FALSE
)
write_output(resource_comparison, "software_resource_capability_comparison.csv")

provenance <- data.frame(
  item = c("run_date", "R_version", "score_implementation", "score_implementation_md5",
           "all_gene_statistics", "all_gene_statistics_md5", "cysteine_atlas",
           "cysteine_atlas_md5", "protein_binding", "protein_binding_md5",
           "swissadme_descriptors", "swissadme_descriptors_md5",
           "CovPDB_source", "UniProt_mapping_endpoint",
           "random_seed"),
  value = c(
    as.character(Sys.Date()), R.version.string,
    rel_path(score_file), unname(tools::md5sum(score_file)),
    rel_path(all_gene_stats_file), unname(tools::md5sum(all_gene_stats_file)),
    rel_path(atlas_file), unname(tools::md5sum(atlas_file)),
    rel_path(binding_file), unname(tools::md5sum(binding_file)),
    rel_path(adme_file), unname(tools::md5sum(adme_file)),
    "https://drug-discovery.vm.uni-freiburg.de/covpdb/proteins_list_by_id/search_type=by_residue_idsearch_id=3",
    "https://rest.uniprot.org/uniprotkb/search", "20260803"
  ),
  stringsAsFactors = FALSE
)
write_output(provenance, "benchmark_provenance.csv")

# The project carries no lockfile, so the run records its own environment.
# docs/R_ENVIRONMENT.md holds the curated version; this is the whole session,
# transitive dependencies included.
writeLines(utils::capture.output(utils::sessionInfo()),
           file.path(output_dir, "session_info.txt"))

# ---- Publication-ready figures --------------------------------

scenario_labels <- c(
  adme_active = "ADME active (0.5)",
  equal_six = "Equal weights, six dimensions",
  equal_five = "Equal weights, five dimensions",
  leave_out_dependency_strength = "Leave out dependency",
  leave_out_cancer_selectivity = "Leave out selectivity",
  leave_out_cysteine_ligandability = "Leave out ligandability",
  leave_out_conservation = "Leave out conservation",
  leave_out_clinical_evidence = "Leave out clinical evidence",
  half_dependency_strength = "Dependency ×0.5",
  double_dependency_strength = "Dependency ×2",
  half_cancer_selectivity = "Selectivity ×0.5",
  double_cancer_selectivity = "Selectivity ×2",
  half_cysteine_ligandability = "Ligandability ×0.5",
  double_cysteine_ligandability = "Ligandability ×2",
  half_conservation = "Conservation ×0.5",
  double_conservation = "Conservation ×2",
  half_clinical_evidence = "Clinical evidence ×0.5",
  double_clinical_evidence = "Clinical evidence ×2"
)
heatmap_data <- robustness[robustness$top_definition == "top_25", ]
heatmap_data$scenario_label <- scenario_labels[heatmap_data$scenario]
heatmap_data$scenario_label <- factor(heatmap_data$scenario_label,
                                      levels = rev(unname(scenario_labels)))
heatmap_data$cohort_code <- factor(heatmap_data$cohort_code,
                                   levels = cohort_keys$cohort_code)

p_robustness <- ggplot(heatmap_data, aes(cohort_code, scenario_label, fill = overlap_fraction)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradient(limits = c(0, 1), low = "#D1495B", high = "#264653",
                      oob = scales::squish, name = "Top-25 overlap") +
  labs(x = "Cancer cohort", y = NULL,
       title = "CPT ranking robustness to alternative evidence weights",
       subtitle = "Fraction of the default top 25 recovered within each cohort") +
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold"), legend.position = "right")
ggsave(file.path(figure_dir, "cpt_weight_robustness_heatmap.png"), p_robustness,
       width = 10.5, height = 6.8, dpi = 300)
ggsave(file.path(figure_dir, "cpt_weight_robustness_heatmap.pdf"), p_robustness,
       width = 10.5, height = 6.8)

ranker_labels <- c(
  CPT_default = "CPT Score (default)",
  CPT_adme_active = "CPT Score (ADME active, 0.5)",
  CPT_equal = "CPT Score (equal weights, five dimensions)",
  Dependency = "Dependency only",
  Selectivity = "Selectivity only",
  Dependency_selectivity_mean = "Dependency + selectivity",
  Cysteine_ligandability = "Cysteine ligandability only"
)
enrichment_plot_data <- enrichment_summary[enrichment_summary$top_fraction == 0.10, ]
enrichment_plot_data$ranker_label <- factor(
  ranker_labels[enrichment_plot_data$ranker],
  levels = rev(unname(ranker_labels))
)
p_enrichment <- ggplot(enrichment_plot_data,
                       aes(fold_enrichment, ranker_label, color = ranker == "CPT_default")) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey55") +
  geom_errorbar(aes(xmin = bootstrap_ci_low, xmax = bootstrap_ci_high),
                orientation = "y", width = 0.16, linewidth = 0.7) +
  geom_point(size = 2.8) +
  scale_color_manual(values = c(`TRUE` = "#D1495B", `FALSE` = "#386FA4"), guide = "none") +
  labs(x = "Fold enrichment of known CovPDB targets in the top 10%", y = NULL,
       title = "Recovery of structurally observed human cysteine-covalent targets",
       subtitle = "Points show pooled enrichment; intervals are cohort-bootstrap 95% CIs") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.major.y = element_blank(), plot.title = element_text(face = "bold"))
ggsave(file.path(figure_dir, "covpdb_positive_set_enrichment.png"), p_enrichment,
       width = 8.6, height = 4.6, dpi = 300)
ggsave(file.path(figure_dir, "covpdb_positive_set_enrichment.pdf"), p_enrichment,
       width = 8.6, height = 4.6)

site_plot_data <- site_summary[site_summary$scope == "Molecules with <=20 targets", ]
site_plot_data <- rbind(
  data.frame(panel = site_plot_data$panel, evidence = "Retained by exact-site evidence",
             n = site_plot_data$exact_site_functional_gene_associations),
  data.frame(panel = site_plot_data$panel, evidence = "Unsupported gene-level upgrades",
             n = site_plot_data$gene_level_false_upgrades)
)
site_plot_data$evidence <- factor(
  site_plot_data$evidence,
  levels = c("Unsupported gene-level upgrades", "Retained by exact-site evidence")
)
p_site <- ggplot(site_plot_data, aes(panel, n, fill = evidence)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5), color = "white",
            fontface = "bold", size = 4) +
  scale_fill_manual(values = c("Unsupported gene-level upgrades" = "#D1495B",
                               "Retained by exact-site evidence" = "#2A9D8F")) +
  labs(x = NULL, y = "Gene-cohort associations",
       fill = NULL, title = "Exact-site matching resolves gene-level evidence upgrades",
       subtitle = "Covalent molecules with no more than 20 annotated targets") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.major.x = element_blank(), plot.title = element_text(face = "bold"),
        legend.position = "bottom")
ggsave(file.path(figure_dir, "exact_site_discrimination.png"), p_site,
       width = 7.4, height = 4.8, dpi = 300)
ggsave(file.path(figure_dir, "exact_site_discrimination.pdf"), p_site,
       width = 7.4, height = 4.8)

cat(sprintf(
  paste0("Benchmark complete: %d cohorts, %d dependency associations, %d CovPDB targets, ",
         "%d weight scenarios.\n"),
  nrow(cohort_keys), nrow(baseline_scores), length(known_positive_genes), length(scenarios)
))
