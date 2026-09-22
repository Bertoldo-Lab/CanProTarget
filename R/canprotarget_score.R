# ============================================================
# Script:   canprotarget_score.R
# Purpose:  CanProTarget Score: composite target prioritization metric.
#           Combines multiple evidence axes into a single ranked score
#           for identifying druggable cancer dependencies.
#
# Method:   Rank-based aggregation (percentile rank per dimension,
#           weighted sum). Robust to different scales across data types.
#
# Usage:    Called by MCP API, Shiny modules, and report templates.
#           No Shiny dependencies; pure R functions.
# ============================================================
# Sections:
#   Default weights
#     cpt_coerce_weights()
#   Percentile rank helper
#     percentile_rank()
#   Dimension scoring functions
#     score_dependency_strength(), score_cancer_selectivity(), score_cysteine_ligandability(), score_conservation(), +1 more
#   ADME / developability
#     cpt_window_score(), cpt_probe_developability(), cpt_build_adme_gene_scores(), score_adme_druggability()
#   Main scoring function
#     cpt_score()
#   Subtype gene-effect background
#     cpt_subtype_background(), cpt_dimension_report(), cpt_weighted_composite(), cpt_priority_label(), +6 more
#   Single-gene convenience function
#     cpt_score_single(), cpt_rank_targets()
#   Residue-resolved engagement, evidence tiers, Site CPT
#     cpt_evidence_tier(), cpt_annotate_engaged_tiers(), cpt_evidence_tier_label(), cpt_build_smcl_index(), +8 more
# ============================================================

# ---- Default weights -------------------------------------------
# Each dimension contributes proportionally. Weights are normalized
# internally so they sum to 1. Users can override via the weights param.

CPT_DEFAULT_WEIGHTS <- list(
  dependency_strength = 3.0,
  cancer_selectivity  = 3.0,
  cysteine_ligandability = 2.0,
  conservation        = 1.5,
  clinical_evidence   = 1.0,
  # ADME is COMPUTED but weighted 0 by default, so composite scores are
  # unchanged from the pre-implementation baseline. The value is still reported
  # (dim_adme, radar, MCP JSON) so it can be inspected before it counts.
  # Raising this changes every score, ranking and report: at 0.5 it is
  # enough to move VCP above KRAS for top rank in PDAC. Treat the default as a
  # methods decision requiring sign-off, not a display preference.
  # See docs/CPT_SCORE.md section 4.6.
  adme_druggability   = 0.0
)

#' Suggested ADME weight once the formula is signed off (was the original default).
CPT_ADME_PROPOSED_WEIGHT <- 0.5

#' Human-readable labels for CPT weight controls (UI / docs).
CPT_WEIGHT_LABELS <- c(
  dependency_strength = "Dependency strength",
  cancer_selectivity = "Cancer selectivity",
  cysteine_ligandability = "Cysteine ligandability",
  conservation = "Conservation",
  clinical_evidence = "Clinical evidence (ClinVar)",
  adme_druggability = "ADME / drug-likeness"
)

#' Coerce a named list/vector of weights to the canonical CPT weight list.
#' Missing keys filled from defaults; non-finite / negative values clamped to 0.
#' @param weights Named list or numeric vector (optional)
#' @return Named list matching CPT_DEFAULT_WEIGHTS keys
cpt_coerce_weights <- function(weights = NULL) {
  out <- CPT_DEFAULT_WEIGHTS
  if (is.null(weights)) return(out)
  w <- unlist(weights)
  if (is.null(names(w)) || !length(w)) return(out)
  for (nm in names(out)) {
    if (nm %in% names(w)) {
      val <- suppressWarnings(as.numeric(w[[nm]]))
      if (length(val) && is.finite(val[1]) && val[1] >= 0) {
        out[[nm]] <- val[1]
      }
    }
  }
  # If every weight is 0, fall back to defaults (avoid all-NA composites)
  if (sum(unlist(out)) <= 0) return(CPT_DEFAULT_WEIGHTS)
  out
}

# ---- Percentile rank helper ------------------------------------

#' Convert a numeric vector to percentile ranks (0-100).
#' Higher percentile = more favorable for targeting.
#' @param x Numeric vector
#' @param higher_is_better Logical. If FALSE, lower raw values get higher ranks.
#' @return Numeric vector of percentile ranks (0-100), NA preserved.
percentile_rank <- function(x, higher_is_better = TRUE) {
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  valid <- !is.na(x)
  n_valid <- sum(valid)
  ranks <- rep(NA_real_, length(x))
  # Single value: return 50 (median, no ranking possible)
  if (n_valid == 1L) {
    ranks[valid] <- 50
    return(ranks)
  }
  if (!higher_is_better) {
    # Lower values are better (e.g., gene effect: more negative = more essential)
    ranks[valid] <- (rank(-x[valid], ties.method = "average") - 1) / (n_valid - 1) * 100
  } else {
    ranks[valid] <- (rank(x[valid], ties.method = "average") - 1) / (n_valid - 1) * 100
  }
  ranks
}

# ---- Dimension scoring functions -------------------------------

#' Score dependency strength.
#' More negative gene effect = stronger dependency = higher score.
#' @param cancer_mean Numeric vector of mean gene effect in the cancer subtype
#' @return Percentile rank (0-100)
score_dependency_strength <- function(cancer_mean) {
  percentile_rank(cancer_mean, higher_is_better = FALSE)
}

#' Score cancer selectivity.
#' More negative effect size (cancer vs others) = more selective = higher score.
#' @param effect_size Numeric vector of effect sizes (cancer mean - other mean)
#' @return Percentile rank (0-100)
score_cancer_selectivity <- function(effect_size) {
  percentile_rank(effect_size, higher_is_better = FALSE)
}

#' Score cysteine ligandability.
#' Uses the best (max) ligandability score for any cysteine on the gene.
#' @param gene_symbols Character vector of gene symbols
#' @param cys_atlas Data frame: cysteine editing atlas
#' @return Named numeric vector of ligandability scores (0-100), NA if no cysteine data
score_cysteine_ligandability <- function(gene_symbols, cys_atlas) {
  if (is.null(cys_atlas) || nrow(cys_atlas) == 0) {
    return(rep(NA_real_, length(gene_symbols)))
  }

  # For each gene, take the max ligandability score across all its cysteine sites
  gene_upper <- toupper(gene_symbols)
  atlas_upper <- toupper(cys_atlas$gene_symbol)

  scores <- vapply(gene_upper, function(g) {
    hits <- cys_atlas[atlas_upper == g, , drop = FALSE]
    if (nrow(hits) == 0) return(NA_real_)  # not in atlas: missing, drop dimension

    # Composite: ligandability_score if available, otherwise binary ligandable flag.
    # Explicit non-ligandable (FALSE) → 0; all-NA flags → NA (no evidence, not penalty).
    if ("ligandability_score" %in% colnames(hits) && any(!is.na(hits$ligandability_score))) {
      max(hits$ligandability_score, na.rm = TRUE)
    } else if ("ligandable" %in% colnames(hits)) {
      lig_vals <- hits$ligandable
      if (all(is.na(lig_vals))) {
        NA_real_
      } else if (any(lig_vals, na.rm = TRUE)) {
        50
      } else {
        0
      }
    } else if ("functional_ligandable" %in% colnames(hits)) {
      fl_vals <- hits$functional_ligandable
      if (all(is.na(fl_vals))) {
        NA_real_
      } else if (any(fl_vals, na.rm = TRUE)) {
        75
      } else {
        0
      }
    } else {
      # Gene is in atlas but no ligandability columns at all
      NA_real_
    }
  }, numeric(1))

  # Already on 0-100 scale from ligandability_score; just return as-is
  scores
}

#' Score evolutionary conservation.
#' Higher conservation score = more conserved = higher score.
#' @param gene_symbols Character vector of gene symbols
#' @param cys_atlas Data frame: cysteine editing atlas (has conservation_score)
#' @return Numeric vector (0-100 scale)
score_conservation <- function(gene_symbols, cys_atlas) {
  if (is.null(cys_atlas) || nrow(cys_atlas) == 0 ||
      !"conservation_score" %in% colnames(cys_atlas)) {
    return(rep(NA_real_, length(gene_symbols)))
  }

  gene_upper <- toupper(gene_symbols)
  atlas_upper <- toupper(cys_atlas$gene_symbol)

  scores <- vapply(gene_upper, function(g) {
    hits <- cys_atlas[atlas_upper == g, , drop = FALSE]
    if (nrow(hits) == 0) return(NA_real_)
    cons <- hits$conservation_score[!is.na(hits$conservation_score)]
    if (!length(cons)) return(NA_real_)
    # Max conservation across cysteines for this gene
    max(cons)
  }, numeric(1))

  # Conservation is typically 0-1; scale to 0-100
  scores * 100
}

#' Score clinical evidence (ClinVar pathogenic annotations).
#' Binary: does any cysteine on this gene have ClinVar pathogenic annotation?
#' @param gene_symbols Character vector of gene symbols
#' @param cys_atlas Data frame: cysteine editing atlas
#' @return Numeric vector: 100 if pathogenic annotation exists, 0 otherwise, NA if no data
score_clinical_evidence <- function(gene_symbols, cys_atlas) {
  if (is.null(cys_atlas) || nrow(cys_atlas) == 0 ||
      !"clinvar_pathogenic" %in% colnames(cys_atlas)) {
    return(rep(NA_real_, length(gene_symbols)))
  }

  gene_upper <- toupper(gene_symbols)
  atlas_upper <- toupper(cys_atlas$gene_symbol)

  vapply(gene_upper, function(g) {
    hits <- cys_atlas[atlas_upper == g, , drop = FALSE]
    if (nrow(hits) == 0) return(NA_real_)  # not in atlas → missing
    cv <- hits$clinvar_pathogenic
    if (all(is.na(cv))) return(NA_real_)   # no ClinVar annotation → missing
    if (any(cv, na.rm = TRUE)) 100 else 0  # in atlas, annotated, no pathogenic → 0
  }, numeric(1))
}

# ---- ADME / developability ------------------------------------
#
# PROVISIONAL FORMULA — pending methods sign-off. See docs/CPT_SCORE.md 4.6.
#
# Design decisions, all deliberate:
#
#  1. Reactivity filters are EXCLUDED. Brenk flags 971/1000 probes and PAINS is
#     near-constant; the Brenk alert *is* the covalent warhead (acrylamide /
#     chloroacetamide). Rewarding "fewer structural alerts" would penalise
#     compounds for being covalent probes, inverting the platform's purpose.
#
#  2. Pass/fail drug-likeness rules are EXCLUDED. Over this fragment library
#     (median MW 235) they carry no signal: 991/1000 probes have zero Lipinski
#     violations, 993/1000 share one bioavailability score, 989/1000 are "High"
#     GI absorption. Including them would add a near-constant term that dilutes
#     the biological dimensions.
#
#  3. Only continuous physicochemical descriptors are used, scored against
#     FRAGMENT-appropriate windows (rule of three: MW <= 300, cLogP <= 3) rather
#     than Lipinski's rule of five, which fragments pass trivially.
#
#  4. Aggregation is the MEAN across qualifying probes, not the max. Taking the
#     best probe looks intuitive but is confounded by assay coverage: it saturates
#     at 100 for 75.8% of genes and correlates +0.52 (Spearman) with the number of
#     probes hitting the gene, so it measures how well-screened a target is rather
#     than how developable its chemistry is. The mean decouples from coverage
#     (+0.025) and reads as "how well does the engaging chemotype sit in
#     fragment-appropriate property space".
#
#  5. Qualifying = CR >= 4, reusing the engagement threshold already used for
#     `n_targets` in the preprocessing pipeline. Genes with no probe reaching it
#     score NA (no credible engagement, therefore no ADME claim) rather than 0.
#
#  KNOWN LIMITATION: the resulting spread is narrow (IQR ~92.6-98.7, 19.4% at
#  100). Fragments really are physicochemically benign, so this dimension shifts
#  scores more than it reorders them. Do not expect it to discriminate strongly.
#  An alternative that spreads slightly wider, weighting the mean toward selective
#  probes via `n_targets`, is left unimplemented because it mixes pharmacology into
#  a pharmacokinetics axis. That is a call for the methods review.

#' Taper a value to 0-100 across an ideal window with linear shoulders.
#' Full marks inside [ideal_lo, ideal_hi]; falls linearly to 0 at zero_lo / zero_hi.
#' @return Numeric vector 0-100, NA preserved
cpt_window_score <- function(x, ideal_lo, ideal_hi, zero_lo, zero_hi) {
  x <- suppressWarnings(as.numeric(x))
  out <- rep(NA_real_, length(x))
  ok <- !is.na(x)
  if (!any(ok)) return(out)
  v <- x[ok]
  s <- rep(100, length(v))
  # low shoulder
  lo <- v < ideal_lo
  s[lo] <- 100 * (v[lo] - zero_lo) / (ideal_lo - zero_lo)
  # high shoulder
  hi <- v > ideal_hi
  s[hi] <- 100 * (zero_hi - v[hi]) / (zero_hi - ideal_hi)
  out[ok] <- pmin(100, pmax(0, s))
  out
}

#' Per-probe fragment developability score (0-100).
#' Mean of the available physicochemical components; NA if none available.
#' @param adme Data frame of SwissADME descriptors (one row per probe)
#' @return Numeric vector aligned to rows of `adme`
cpt_probe_developability <- function(adme) {
  if (is.null(adme) || !nrow(adme)) return(numeric(0))

  comp <- list(
    # Fragment rule of three: <=300 ideal, useless by 500
    mw  = if ("MW" %in% names(adme)) {
      cpt_window_score(adme$MW, 0, 300, -Inf, 500)
    } else NULL,
    # Lipophilicity: 0-3 ideal. Too polar cannot cross membranes, too greasy
    # brings solubility and promiscuity problems.
    logp = if ("Consensus.Log.P" %in% names(adme)) {
      cpt_window_score(adme$Consensus.Log.P, 0, 3, -3, 6)
    } else NULL,
    # Polar surface area: <=90 comfortable for oral, degraded by 140 (Veber)
    tpsa = if ("TPSA" %in% names(adme)) {
      cpt_window_score(adme$TPSA, 0, 90, -Inf, 140)
    } else NULL,
    # Synthetic accessibility (1 easy - 10 hard)
    sa  = if ("Synthetic.Accessibility" %in% names(adme)) {
      cpt_window_score(adme$Synthetic.Accessibility, 0, 3, -Inf, 7)
    } else NULL
  )
  comp <- comp[!vapply(comp, is.null, logical(1))]
  if (!length(comp)) return(rep(NA_real_, nrow(adme)))

  m <- do.call(cbind, comp)
  apply(m, 1, function(r) if (all(is.na(r))) NA_real_ else mean(r, na.rm = TRUE))
}

#' Build the per-gene ADME lookup from chemoproteomics + SwissADME.
#'
#' Collapses the ~10.6M-row binding table into one row per gene so scoring stays
#' fast. Build once at startup and pass the result to `score_adme_druggability()`.
#'
#' @param binding Data frame: protein_binding_lookup_preprocessed.rds
#' @param adme Data frame: swissadme_preprocessed.rds
#' @param min_cr Numeric: minimum competition ratio to count as engagement (default 4)
#' @return Data frame: gene_key, adme_score (mean developability across qualifying
#'         probes), n_probes, best_probe (most developable single probe, reported
#'         for interpretability only — it does not drive adme_score)
#' Canonical probe name.
#'
#' The chemoproteomics tables carry SwissADME's short codes (AC5, CL174); the
#' SwissADME file on disk carries the publication names (ACRYL_5, CL_174). The
#' two have to be spelled the same way before anything joins them. The Shiny
#' loader normalises on read, but callers that read the RDS directly do not, so
#' the scorer canonicalises its own inputs rather than trusting them.
cpt_canonical_probe_name <- function(x) {
  x <- as.character(x)
  x <- sub("^CL_([0-9]+)$", "CL\\1", x)
  x <- sub("^ACRYL_([0-9]+)$", "AC\\1", x)
  x
}

cpt_build_adme_gene_scores <- function(binding, adme, min_cr = 4) {
  empty <- data.frame(gene_key = character(0), adme_score = numeric(0),
                      n_probes = integer(0), best_probe = character(0),
                      stringsAsFactors = FALSE)
  if (is.null(binding) || is.null(adme) || !nrow(binding) || !nrow(adme)) return(empty)
  if (!all(c("probe_name", "CR") %in% names(binding))) return(empty)
  if (!"probe_name" %in% names(adme)) return(empty)

  gene_col <- if ("gene_name_key" %in% names(binding)) "gene_name_key" else "gene_name"
  if (!gene_col %in% names(binding)) return(empty)

  probe_score <- data.frame(
    probe_name = cpt_canonical_probe_name(adme$probe_name),
    dev = cpt_probe_developability(adme),
    stringsAsFactors = FALSE
  )
  probe_score <- probe_score[!is.na(probe_score$dev), , drop = FALSE]
  probe_score <- probe_score[!duplicated(probe_score$probe_name), , drop = FALSE]
  if (!nrow(probe_score)) return(empty)

  keep <- !is.na(binding$CR) & binding$CR >= min_cr
  b <- binding[keep, c(gene_col, "probe_name"), drop = FALSE]
  names(b)[1] <- "gene_key"
  b$probe_name <- cpt_canonical_probe_name(b$probe_name)
  b$gene_key <- toupper(as.character(b$gene_key))
  b <- b[!is.na(b$gene_key) & nzchar(b$gene_key), , drop = FALSE]
  if (!nrow(b)) return(empty)

  b <- merge(unique(b), probe_score, by = "probe_name")
  if (!nrow(b)) return(empty)

  # adme_score = MEAN developability across qualifying probes (see note 4 above).
  agg <- stats::aggregate(dev ~ gene_key, data = b, FUN = mean)
  names(agg)[2] <- "adme_score"
  cnt <- stats::aggregate(probe_name ~ gene_key, data = b,
                          FUN = function(x) length(unique(x)))
  names(cnt)[2] <- "n_probes"

  # best_probe is reported for interpretability only; it does not drive the score.
  ord <- order(b$gene_key, -b$dev)
  bb <- b[ord, , drop = FALSE]
  top <- bb[!duplicated(bb$gene_key), c("gene_key", "probe_name"), drop = FALSE]
  names(top)[2] <- "best_probe"

  out <- merge(merge(agg, cnt, by = "gene_key"), top, by = "gene_key")
  out$n_probes <- as.integer(out$n_probes)
  rownames(out) <- NULL
  out[, c("gene_key", "adme_score", "n_probes", "best_probe")]
}

#' Score ADME / developability of the best covalent probe engaging each gene.
#'
#' @param gene_symbols Character vector of gene symbols
#' @param adme_data Per-gene lookup from `cpt_build_adme_gene_scores()`.
#'        NULL, or any frame lacking `gene_key`/`adme_score`, yields all-NA so
#'        callers without the chemoproteomics tables keep working unchanged.
#' @return Numeric vector (0-100, or NA where no probe reaches the CR threshold)
score_adme_druggability <- function(gene_symbols, adme_data = NULL) {
  n <- length(gene_symbols)
  if (is.null(adme_data) || !is.data.frame(adme_data) || !nrow(adme_data) ||
      !all(c("gene_key", "adme_score") %in% names(adme_data))) {
    return(rep(NA_real_, n))
  }
  idx <- match(toupper(sub(" \\(\\d+\\)$", "", as.character(gene_symbols))),
               adme_data$gene_key)
  as.numeric(adme_data$adme_score[idx])
}

# ---- Main scoring function -------------------------------------

#' Compute CanProTarget Score for genes in a cancer subtype.
#'
#' Combines multiple evidence dimensions into a single composite score
#' using rank-based aggregation with configurable weights.
#'
#' @param gene_stats Data frame with columns: gene_name, cancer_mean (or Cancer_Avg),
#'        effect_size (or EffectSize). From ge_analysis() or api_query_dependency().
#' @param cys_atlas Data frame: cysteine editing atlas (from data/cys_editing_atlas.rds)
#' @param adme_data Data frame: SwissADME data (NULL if unavailable)
#' @param weights Named list of dimension weights (defaults to CPT_DEFAULT_WEIGHTS)
#' @param min_dimensions Integer: minimum non-NA dimensions required for a valid score (default 2)
#' @return Data frame with gene_name, dimension scores, composite CPT score, and rank
cpt_score <- function(gene_stats,
                      cys_atlas = NULL,
                      adme_data = NULL,
                      weights = CPT_DEFAULT_WEIGHTS,
                      min_dimensions = 2L) {

  if (is.null(gene_stats) || nrow(gene_stats) == 0) {
    return(data.frame(
      gene_name = character(0),
      cpt_score = numeric(0),
      cpt_rank = integer(0),
      stringsAsFactors = FALSE
    ))
  }

  # Normalize column names
  df <- gene_stats
  if ("Cancer_Avg" %in% colnames(df) && !"cancer_mean" %in% colnames(df)) {
    df$cancer_mean <- df$Cancer_Avg
  }
  if ("EffectSize" %in% colnames(df) && !"effect_size" %in% colnames(df)) {
    df$effect_size <- df$EffectSize
  }

  genes <- df$gene_name
  n <- length(genes)

  # Compute each dimension
  dim_dep <- if ("cancer_mean" %in% colnames(df)) {
    score_dependency_strength(df$cancer_mean)
  } else {
    rep(NA_real_, n)
  }

  dim_sel <- if ("effect_size" %in% colnames(df)) {
    score_cancer_selectivity(df$effect_size)
  } else {
    rep(NA_real_, n)
  }

  dim_lig <- score_cysteine_ligandability(genes, cys_atlas)
  dim_con <- score_conservation(genes, cys_atlas)
  dim_cli <- score_clinical_evidence(genes, cys_atlas)
  dim_adme <- score_adme_druggability(genes, adme_data)

  # Build dimension matrix
  dim_matrix <- data.frame(
    dependency_strength = dim_dep,
    cancer_selectivity = dim_sel,
    cysteine_ligandability = dim_lig,
    conservation = dim_con,
    clinical_evidence = dim_cli,
    adme_druggability = dim_adme,
    stringsAsFactors = FALSE
  )

  # Normalize weights (only for dimensions that have ANY data)
  w <- unlist(weights)
  dim_names <- colnames(dim_matrix)
  w <- w[dim_names]

  # For each gene, compute weighted average of available dimensions
  composite <- vapply(seq_len(n), function(i) {
    vals <- as.numeric(dim_matrix[i, ])
    available <- !is.na(vals)
    n_avail <- sum(available)
    if (n_avail < min_dimensions) return(NA_real_)
    w_avail <- w[available]
    w_norm <- w_avail / sum(w_avail)
    sum(vals[available] * w_norm)
  }, numeric(1))

  # Build result
  result <- data.frame(
    gene_name = genes,
    cpt_score = round(composite, 2),
    dim_dependency = round(dim_dep, 1),
    dim_selectivity = round(dim_sel, 1),
    dim_ligandability = round(dim_lig, 1),
    dim_conservation = round(dim_con, 1),
    dim_clinical = round(dim_cli, 1),
    dim_adme = round(dim_adme, 1),
    n_dimensions = rowSums(!is.na(dim_matrix)),
    stringsAsFactors = FALSE
  )

  # Rank by composite score (highest = best target)
  result$cpt_rank <- NA_integer_
  scored <- !is.na(result$cpt_score)
  if (any(scored)) {
    result$cpt_rank[scored] <- rank(-result$cpt_score[scored], ties.method = "min")
  }

  # Sort by rank

  result <- result[order(result$cpt_rank, na.last = TRUE), ]
  rownames(result) <- NULL
  result
}

# ---- Subtype gene-effect background ----------------------------

#' Mean gene effect and effect size for every gene in a subtype.
#' Used so percentile ranks are relative to the full gene universe,
#' not a single gene (which always ranks at the median).
#'
#' @param matrix Gene effect matrix (models x genes)
#' @param meta Cancer model metadata with ModelID + OncotreeSubtype
#' @param subtype Character: OncotreeSubtype label
#' @return List with gene_names, cancer_means, effect_sizes, n_cell_lines
cpt_subtype_background <- function(matrix, meta, subtype) {
  meta_aligned <- meta[match(rownames(matrix), meta$ModelID), ]
  cancer_idx <- which(meta_aligned$OncotreeSubtype %in% subtype)
  other_idx <- which(!meta_aligned$OncotreeSubtype %in% subtype)

  if (length(cancer_idx) < 3L) {
    stop("Subtype '", subtype, "' has <3 cell lines.", call. = FALSE)
  }

  cancer_means <- colMeans(matrix[cancer_idx, , drop = FALSE], na.rm = TRUE)
  other_means <- colMeans(matrix[other_idx, , drop = FALSE], na.rm = TRUE)

  list(
    gene_names = sub(" \\(\\d+\\)$", "", colnames(matrix)),
    cancer_means = as.numeric(cancer_means),
    effect_sizes = as.numeric(cancer_means - other_means),
    n_cell_lines = length(cancer_idx)
  )
}

#' Build dimension metadata for JSON (active vs missing, with caveats).
#' @param dim_scores Named list of dimension scores (may be NA)
#' @param raw_values Named list of optional raw values
#' @return List with dimensions, active, missing, n_used
cpt_dimension_report <- function(dim_scores, raw_values = list()) {
  descriptions <- list(
    dependency_strength = "Mean gene effect in subtype (more negative = more essential). Score is percentile among all genes in this subtype.",
    cancer_selectivity = "Effect size vs other cancers (more negative = more selective). Score is percentile among all genes in this subtype.",
    cysteine_ligandability = "Best ligandability evidence across the gene's cysteine sites (atlas).",
    conservation = "Evolutionary conservation of targetable cysteine (atlas).",
    clinical_evidence = "ClinVar pathogenic annotations at cysteine sites (atlas).",
    adme_druggability = paste(
      "Provisional developability: mean fragment physicochemical score (MW, cLogP,",
      "TPSA, synthetic accessibility) across covalent probes engaging the gene at",
      "CR >= 4. Reactivity filters deliberately excluded. Weight 0 by default."
    )
  )

  dims <- list()
  active <- character(0)
  missing <- character(0)

  for (nm in names(dim_scores)) {
    sc <- dim_scores[[nm]]
    entry <- list(
      score = if (is.null(sc) || (length(sc) == 1 && is.na(sc))) NULL else round(as.numeric(sc), 1),
      description = descriptions[[nm]] %||% nm,
      active = !(is.null(sc) || (length(sc) == 1 && is.na(sc)))
    )
    if (!is.null(raw_values[[nm]]) && !is.na(raw_values[[nm]])) {
      entry$raw_value <- round(as.numeric(raw_values[[nm]]), 4)
    }
    dims[[nm]] <- entry
    if (isTRUE(entry$active)) {
      active <- c(active, nm)
    } else {
      missing <- c(missing, nm)
    }
  }

  list(
    dimensions = dims,
    dimensions_active = active,
    dimensions_missing = missing,
    n_dimensions_used = length(active)
  )
}

# Null-coalesce (local; mcp_worker also defines one)
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Weighted composite from available dimension scores.
#' @param dim_scores Named numeric (NA allowed)
#' @param weights Named weights list
#' @param min_dimensions Minimum non-NA dims required
#' @return Scalar composite or NA
cpt_weighted_composite <- function(dim_scores, weights = CPT_DEFAULT_WEIGHTS, min_dimensions = 2L) {
  w <- unlist(weights)
  vals <- unlist(dim_scores)
  # Align names
  common <- intersect(names(vals), names(w))
  vals <- vals[common]
  w <- w[common]
  available <- !is.na(vals)
  if (sum(available) < min_dimensions) return(NA_real_)
  w_avail <- w[available]
  w_norm <- w_avail / sum(w_avail)
  sum(vals[available] * w_norm)
}

#' Priority band label for a CPT score (heuristic bands, not clinical cutoffs).
cpt_priority_label <- function(score) {
  if (is.na(score)) return("unscored")
  if (score >= 70) return("high")
  if (score >= 40) return("moderate")
  "low"
}

#' Apply non-dependency demotion rules to a priority label (shared by MCP + reports).
#' Non-dependencies are never "high"; non-dep + non-selective forced to "low".
cpt_apply_priority_guards <- function(priority,
                                      is_dependency,
                                      is_selective_by_effect_size = TRUE) {
  p <- as.character(priority)[1]
  if (is.na(p) || !nzchar(p)) p <- "unscored"
  dep <- isTRUE(is_dependency)
  sel <- isTRUE(is_selective_by_effect_size)
  if (!dep && identical(p, "high")) p <- "moderate"
  if (!dep && !sel && identical(p, "moderate")) p <- "low"
  p
}

#' Human-readable priority string for reports from a cpt_score_single result.
cpt_priority_display <- function(score_result) {
  if (is.null(score_result)) return("Unscored")
  p <- score_result$priority
  if (is.null(p) || is.na(p) || !nzchar(as.character(p)[1])) {
    # Fallback: recompute with same guards as cpt_score_single
    p <- cpt_priority_label(score_result$cpt_score)
    p <- cpt_apply_priority_guards(
      p,
      is_dependency = isTRUE(score_result$is_dependency),
      is_selective_by_effect_size = isTRUE(score_result$is_selective_by_effect_size)
    )
  }
  switch(as.character(p),
    high = "High priority",
    moderate = "Moderate priority",
    low = "Low priority",
    unscored = "Unscored",
    paste0(toupper(substring(p, 1, 1)), substring(p, 2), " priority")
  )
}

#' Annotate a filtered gene table with CPT scores ranked against a full background.
#' Percentiles for dependency/selectivity use the full gene universe; cysteine dims
#' are computed only for the filtered genes (faster for Shiny tables).
#'
#' @param gene_df Filtered table (e.g. cancer_gene_df) with gene_name + Cancer_Avg/EffectSize
#' @param background_df Full gene table used for percentile ranks (all_gene_ge_df)
#' @param cys_atlas Cysteine atlas or NULL
#' @param weights CPT weight list
#' @return gene_df with CPT_Score, CPT_Rank, and dimension columns prepended after gene_name
cpt_annotate_gene_table <- function(gene_df,
                                    background_df,
                                    cys_atlas = NULL,
                                    weights = CPT_DEFAULT_WEIGHTS,
                                    adme_data = NULL) {
  if (is.null(gene_df) || !nrow(gene_df) || is.null(background_df) || !nrow(background_df)) {
    return(gene_df)
  }

  bg <- background_df
  if ("Cancer_Avg" %in% colnames(bg) && !"cancer_mean" %in% colnames(bg)) {
    bg$cancer_mean <- bg$Cancer_Avg
  }
  if ("EffectSize" %in% colnames(bg) && !"effect_size" %in% colnames(bg)) {
    bg$effect_size <- bg$EffectSize
  }
  if (!all(c("gene_name", "cancer_mean", "effect_size") %in% colnames(bg))) {
    return(gene_df)
  }

  bg$gene_key <- toupper(sub(" \\(\\d+\\)$", "", as.character(bg$gene_name)))
  dim_dep_all <- score_dependency_strength(bg$cancer_mean)
  dim_sel_all <- score_cancer_selectivity(bg$effect_size)
  names(dim_dep_all) <- bg$gene_key
  names(dim_sel_all) <- bg$gene_key

  out <- gene_df
  out$gene_key <- toupper(sub(" \\(\\d+\\)$", "", as.character(out$gene_name)))
  genes_clean <- sub(" \\(\\d+\\)$", "", as.character(out$gene_name))

  dim_dep <- unname(dim_dep_all[out$gene_key])
  dim_sel <- unname(dim_sel_all[out$gene_key])
  dim_lig <- as.numeric(score_cysteine_ligandability(genes_clean, cys_atlas))
  dim_con <- as.numeric(score_conservation(genes_clean, cys_atlas))
  dim_cli <- as.numeric(score_clinical_evidence(genes_clean, cys_atlas))
  dim_adme <- as.numeric(score_adme_druggability(genes_clean, adme_data))

  n <- nrow(out)
  dim_matrix <- data.frame(
    dependency_strength = dim_dep,
    cancer_selectivity = dim_sel,
    cysteine_ligandability = dim_lig,
    conservation = dim_con,
    clinical_evidence = dim_cli,
    adme_druggability = dim_adme,
    stringsAsFactors = FALSE
  )
  w <- unlist(weights)[colnames(dim_matrix)]
  composite <- vapply(seq_len(n), function(i) {
    vals <- as.numeric(dim_matrix[i, ])
    available <- !is.na(vals)
    # Align with MCP single-gene / rank (min 2 dims)
    if (sum(available) < 2L) return(NA_real_)
    w_a <- w[available]
    sum(vals[available] * (w_a / sum(w_a)))
  }, numeric(1))

  out$CPT_Score <- round(composite, 2)
  out$dim_dependency <- round(dim_dep, 1)
  out$dim_selectivity <- round(dim_sel, 1)
  out$dim_ligandability <- round(dim_lig, 1)
  out$dim_conservation <- round(dim_con, 1)
  out$dim_clinical <- round(dim_cli, 1)
  scored <- !is.na(out$CPT_Score)
  out$CPT_Rank <- NA_integer_
  if (any(scored)) {
    out$CPT_Rank[scored] <- rank(-out$CPT_Score[scored], ties.method = "min")
  }
  out$gene_key <- NULL

  # Put CPT columns near the front
  front <- c("gene_name", "CPT_Score", "CPT_Rank",
             "dim_dependency", "dim_selectivity", "dim_ligandability",
             "dim_conservation", "dim_clinical")
  front <- intersect(front, colnames(out))
  rest <- setdiff(colnames(out), front)
  out[, c(front, rest), drop = FALSE]
}

#' Dimension scores as a named numeric vector for radar / sensitivity.
#' @param score_result List from cpt_score_single / api_canprotarget_score
#' @return Named numeric (NA allowed)
cpt_dimension_vector <- function(score_result) {
  dims <- score_result$dimensions
  if (is.null(dims)) return(numeric(0))
  vapply(names(dims), function(nm) {
    sc <- dims[[nm]]$score
    if (is.null(sc)) NA_real_ else as.numeric(sc)
  }, numeric(1))
}

#' CPT dimension radar (ggplot polar).
#' Only non-NA dimensions are drawn so missing evidence is not shown as zero.
#' @param dim_scores Named numeric 0-100
#' @param title Plot title
#' @return ggplot
cpt_dimension_radar_gg <- function(dim_scores, title = "CPT dimensions") {
  if (!length(dim_scores)) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No dimensions") +
        ggplot2::theme_void()
    )
  }
  labels_all <- names(dim_scores)
  vals_all <- as.numeric(dim_scores)
  keep <- !is.na(vals_all)
  n_missing <- sum(!keep)
  labels <- labels_all[keep]
  vals <- vals_all[keep]
  if (!length(vals)) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5,
                          label = "All CPT dimensions missing") +
        ggplot2::theme_void()
    )
  }
  n <- length(vals)
  # close polygon
  ang <- seq(0, 2 * pi, length.out = n + 1)
  r <- c(vals, vals[1]) / 100
  df <- data.frame(
    x = r * cos(ang - pi / 2),
    y = r * sin(ang - pi / 2)
  )
  # axis labels
  lab_r <- 1.15
  lab_ang <- seq(0, 2 * pi, length.out = n + 1)[seq_len(n)]
  lab_df <- data.frame(
    x = lab_r * cos(lab_ang - pi / 2),
    y = lab_r * sin(lab_ang - pi / 2),
    label = paste0(labels, "\n", round(vals, 0))
  )
  # unit circle guides
  circ <- function(radius) {
    th <- seq(0, 2 * pi, length.out = 100)
    data.frame(x = radius * cos(th), y = radius * sin(th))
  }
  sub <- if (n_missing > 0) {
    paste0("Active dims only (0-100); ", n_missing, " missing dim(s) omitted")
  } else {
    "Scores 0-100 (all dimensions active)"
  }
  ggplot2::ggplot() +
    ggplot2::geom_path(data = circ(0.5), ggplot2::aes(x, y),
                       colour = "grey85", linewidth = 0.3) +
    ggplot2::geom_path(data = circ(1), ggplot2::aes(x, y),
                       colour = "grey75", linewidth = 0.4) +
    ggplot2::geom_polygon(data = df, ggplot2::aes(x, y),
                          fill = "#4a7fb5", alpha = 0.35, colour = "#1a3a5c") +
    ggplot2::geom_text(data = lab_df, ggplot2::aes(x, y, label = label),
                       size = 2.8, colour = "#1a3a5c") +
    ggplot2::coord_equal(xlim = c(-1.45, 1.45), ylim = c(-1.45, 1.45)) +
    ggplot2::labs(title = title, subtitle = sub) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, colour = "#1a3a5c"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, colour = "#666666", size = 9)
    )
}

#' Leave-one-weight-perturbed sensitivity of the composite score.
#' Multiplies each weight by factors in turn; reports composite for each scenario.
#'
#' @param dim_scores Named numeric dimension scores
#' @param weights Named weight list (default CPT_DEFAULT_WEIGHTS)
#' @param factors Numeric multipliers (default 0.5, 1, 2)
#' @return data.frame: dimension, factor, weight_used, cpt_score
cpt_weight_sensitivity <- function(dim_scores,
                                   weights = CPT_DEFAULT_WEIGHTS,
                                   factors = c(0.5, 1, 2)) {
  base_w <- unlist(weights)
  vals <- unlist(dim_scores)
  common <- intersect(names(vals), names(base_w))
  vals <- vals[common]
  base_w <- base_w[common]
  available <- !is.na(vals)
  if (sum(available) < 1L) {
    return(data.frame(
      dimension = character(0), factor = numeric(0),
      weight_used = numeric(0), cpt_score = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  score_with <- function(w) {
    w <- w[available]
    v <- vals[available]
    sum(v * (w / sum(w)))
  }

  rows <- list()
  # baseline
  rows[[1]] <- data.frame(
    dimension = "(baseline)",
    factor = 1,
    weight_used = NA_real_,
    cpt_score = round(score_with(base_w), 2),
    stringsAsFactors = FALSE
  )
  k <- 2L
  for (dm in names(base_w)[available]) {
    for (f in factors) {
      if (identical(as.numeric(f), 1) && dm == names(base_w)[available][1]) {
        # baseline already covers factor=1 overall; still record per-dim at 1 for clarity
      }
      w2 <- base_w
      w2[dm] <- base_w[dm] * f
      rows[[k]] <- data.frame(
        dimension = dm,
        factor = f,
        weight_used = unname(w2[dm]),
        cpt_score = round(score_with(w2), 2),
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }
  dplyr::bind_rows(rows)
}

# ---- Single-gene convenience function -------------------------

#' Compute CPT Score for a single gene in a specific subtype.
#' Dependency and selectivity percentiles are ranked against ALL genes
#' in the subtype (not a single-row table, which would always return 50).
#'
#' @param gene Character: gene symbol
#' @param subtype Character: OncotreeSubtype
#' @param dataset Character: "CRISPR" or "RNAi"
#' @param data_env Environment/list with loaded data (from mcp_worker or app)
#' @param weights Named list of dimension weights (defaults to CPT_DEFAULT_WEIGHTS)
#' @return List with score breakdown suitable for JSON serialization
cpt_score_single <- function(gene, subtype, dataset = "CRISPR", data_env,
                             weights = CPT_DEFAULT_WEIGHTS) {
  weights <- cpt_coerce_weights(weights)
  # Prefer shared normalizer when api_functions.R is sourced; fall back for unit tests
  if (exists("cpt_normalize_dataset", mode = "function")) {
    dataset <- cpt_normalize_dataset(dataset)
  } else {
    dataset <- toupper(as.character(dataset)[1])
    if (dataset == "RNAI") dataset <- "RNAi"
  }
  matrix <- if (dataset == "CRISPR") data_env$crispr_matrix else data_env$rnai_matrix
  meta <- data_env$cancer_model_data
  cys_atlas <- data_env$cys_atlas

  if (is.null(matrix)) stop("Dataset '", dataset, "' not loaded.", call. = FALSE)
  if (is.null(meta)) stop("Cancer model metadata not loaded.", call. = FALSE)

  bg <- cpt_subtype_background(matrix, meta, subtype)
  col_idx <- which(toupper(bg$gene_names) == toupper(gene))
  if (!length(col_idx)) stop("Gene '", gene, "' not found in ", dataset, ".", call. = FALSE)
  col_idx <- col_idx[1]
  actual_gene <- bg$gene_names[col_idx]

  cancer_mean <- bg$cancer_means[col_idx]
  effect_size <- bg$effect_sizes[col_idx]

  # Percentiles against full background (fixes single-gene rank = 50 bug)
  dim_dep <- score_dependency_strength(bg$cancer_means)[col_idx]
  dim_sel <- score_cancer_selectivity(bg$effect_sizes)[col_idx]
  dim_lig <- score_cysteine_ligandability(actual_gene, cys_atlas)
  dim_con <- score_conservation(actual_gene, cys_atlas)
  dim_cli <- score_clinical_evidence(actual_gene, cys_atlas)
  dim_adme <- score_adme_druggability(actual_gene, data_env$adme_gene_scores)

  dim_scores <- list(
    dependency_strength = dim_dep,
    cancer_selectivity = dim_sel,
    cysteine_ligandability = dim_lig[1],
    conservation = dim_con[1],
    clinical_evidence = dim_cli[1],
    adme_druggability = dim_adme[1]
  )

  # Require 2 dims when both dep+sel exist (always for genes in matrix); matches rank_targets
  composite <- cpt_weighted_composite(dim_scores, weights = weights, min_dimensions = 2L)
  report <- cpt_dimension_report(
    dim_scores,
    raw_values = list(
      dependency_strength = cancer_mean,
      cancer_selectivity = effect_size
    )
  )

  # Raw thresholds (same as api_query_dependency) so we do not over-label
  # mid-pack non-essentials just because percentiles sit near the median.
  is_dependency <- !is.na(cancer_mean) && cancer_mean < -0.5
  is_selective_raw <- !is.na(effect_size) && effect_size < -0.2

  caveats <- character(0)
  if ("adme_druggability" %in% report$dimensions_missing) {
    caveats <- c(caveats, paste(
      "No ADME/developability score: no covalent probe engages this gene at CR >= 4",
      "in the bundled chemoproteomic screens. Absence of probe data is not evidence",
      "against ligandability."
    ))
  } else {
    caveats <- c(caveats, paste(
      "ADME/developability is provisional: mean fragment physicochemical score of",
      "covalent probes engaging this gene at CR >= 4. It carries weight 0 by default,",
      "so it is reported but does not affect the composite. See docs/CPT_SCORE.md 4.6."
    ))
  }
  if ("cysteine_ligandability" %in% report$dimensions_missing) {
    caveats <- c(caveats, "No cysteine ligandability data for this gene in the atlas.")
  }
  if (report$n_dimensions_used < 3L) {
    caveats <- c(caveats, paste0(
      "Only ", report$n_dimensions_used,
      " scoring dimension(s) active; compare carefully with genes scored on more axes."
    ))
  }
  if (!is_dependency) {
    caveats <- c(caveats,
      "Mean gene effect is not below the -0.5 dependency threshold; percentile ranks alone can look middling for non-essentials."
    )
  }
  caveats <- c(caveats,
    "CPT Score is a research prioritization aid, not a clinical recommendation.",
    "Dependency percentiles are within this subtype gene universe only."
  )

  # Priority bands use composite score, but non-dependencies are never "high"
  priority <- cpt_apply_priority_guards(
    cpt_priority_label(composite),
    is_dependency = is_dependency,
    is_selective_by_effect_size = is_selective_raw
  )

  interpretation <- if (is.na(composite)) {
    paste0("Insufficient data to score ", actual_gene, " in ", subtype)
  } else if (!is_dependency) {
    paste0(
      actual_gene, " has CPT Score ", round(composite, 2),
      "/100 in ", subtype, " but is NOT a dependency by mean gene effect (",
      round(cancer_mean, 3), "); treat priority as ", toupper(priority),
      " at best (", report$n_dimensions_used, " of 6 dimensions active)."
    )
  } else {
    paste0(
      actual_gene, " is a ", toupper(priority), "-PRIORITY target in ", subtype,
      " (CPT Score: ", round(composite, 2), "/100 using ",
      report$n_dimensions_used, " of 6 dimensions)"
    )
  }

  list(
    gene = actual_gene,
    subtype = subtype,
    dataset = dataset,
    n_cell_lines_in_subtype = bg$n_cell_lines,
    # Keep NA (not NULL) so report templates can use is.na(); jsonlite maps NA to null
    cpt_score = if (is.na(composite)) NA_real_ else round(composite, 2),
    priority = priority,
    is_dependency = is_dependency,
    is_selective_by_effect_size = is_selective_raw,
    dimensions = report$dimensions,
    dimensions_active = report$dimensions_active,
    dimensions_missing = report$dimensions_missing,
    n_dimensions_used = report$n_dimensions_used,
    weights = lapply(weights, function(x) as.numeric(x)),
    caveats = caveats,
    interpretation = interpretation
  )
}

#' Rank genes in a subtype by CPT Score.
#' Starts from the most selective dependencies, then re-ranks by composite score.
#'
#' @param subtype Character: OncotreeSubtype
#' @param dataset Character: "CRISPR" or "RNAi"
#' @param n Integer: number of results to return
#' @param pool Integer: how many selective deps to score before cutting to n
#' @param require_ligandable Logical: if TRUE, drop genes with no ligandability score
#' @param data_env Environment/list with loaded data
#' @return List suitable for JSON
cpt_rank_targets <- function(subtype,
                             dataset = "CRISPR",
                             n = 20L,
                             pool = 100L,
                             require_ligandable = FALSE,
                             data_env,
                             weights = CPT_DEFAULT_WEIGHTS) {
  weights <- cpt_coerce_weights(weights)
  if (exists("cpt_normalize_dataset", mode = "function")) {
    dataset <- cpt_normalize_dataset(dataset)
  } else {
    dataset <- toupper(as.character(dataset)[1])
    if (dataset == "RNAI") dataset <- "RNAi"
  }
  matrix <- if (dataset == "CRISPR") data_env$crispr_matrix else data_env$rnai_matrix
  meta <- data_env$cancer_model_data
  cys_atlas <- data_env$cys_atlas

  if (is.null(matrix)) stop("Dataset '", dataset, "' not loaded.", call. = FALSE)
  if (is.null(meta)) stop("Cancer model metadata not loaded.", call. = FALSE)

  n <- as.integer(max(1L, min(as.integer(n)[1], 100L)))
  if (is.na(n) || n < 1L) n <- 20L
  pool <- as.integer(max(n, min(as.integer(pool)[1], 500L)))
  if (is.na(pool)) pool <- 100L

  bg <- cpt_subtype_background(matrix, meta, subtype)

  # Full-background percentiles (stable ranks)
  dim_dep_all <- score_dependency_strength(bg$cancer_means)
  dim_sel_all <- score_cancer_selectivity(bg$effect_sizes)

  # Candidate pool: most selective by effect size among genes that look like
  # dependencies (mean < -0.5). Fall back to ES-only pool if too few deps.
  dep_mask <- !is.na(bg$cancer_means) & bg$cancer_means < -0.5
  dep_idx <- which(dep_mask)
  if (length(dep_idx) >= 10L) {
    ord <- dep_idx[order(bg$effect_sizes[dep_idx])]
  } else {
    ord <- order(bg$effect_sizes)  # most negative first
  }
  cand_idx <- ord[seq_len(min(pool, length(ord)))]

  cand_genes <- bg$gene_names[cand_idx]
  dim_lig <- score_cysteine_ligandability(cand_genes, cys_atlas)
  dim_con <- score_conservation(cand_genes, cys_atlas)
  dim_cli <- score_clinical_evidence(cand_genes, cys_atlas)
  dim_adme <- score_adme_druggability(cand_genes, data_env$adme_gene_scores)

  scores <- vapply(seq_along(cand_idx), function(i) {
    idx <- cand_idx[i]
    ds <- list(
      dependency_strength = dim_dep_all[idx],
      cancer_selectivity = dim_sel_all[idx],
      cysteine_ligandability = dim_lig[i],
      conservation = dim_con[i],
      clinical_evidence = dim_cli[i],
      adme_druggability = dim_adme[i]
    )
    if (isTRUE(require_ligandable) && is.na(dim_lig[i])) return(NA_real_)
    cpt_weighted_composite(ds, weights = weights, min_dimensions = 2L)
  }, numeric(1))

  n_dims <- vapply(seq_along(cand_idx), function(i) {
    sum(!is.na(c(
      dim_dep_all[cand_idx[i]],
      dim_sel_all[cand_idx[i]],
      dim_lig[i],
      dim_con[i],
      dim_cli[i],
      dim_adme[i]
    )))
  }, integer(1))

  is_dep <- !is.na(bg$cancer_means[cand_idx]) & bg$cancer_means[cand_idx] < -0.5
  is_sel_es <- !is.na(bg$effect_sizes[cand_idx]) & bg$effect_sizes[cand_idx] < -0.2
  priority <- vapply(seq_along(cand_idx), function(i) {
    cpt_apply_priority_guards(
      cpt_priority_label(scores[i]),
      is_dependency = is_dep[i],
      is_selective_by_effect_size = is_sel_es[i]
    )
  }, character(1))

  tab <- data.frame(
    gene = cand_genes,
    cpt_score = round(scores, 2),
    priority = priority,
    is_dependency = is_dep,
    cancer_mean = round(bg$cancer_means[cand_idx], 4),
    effect_size = round(bg$effect_sizes[cand_idx], 4),
    dim_dependency = round(dim_dep_all[cand_idx], 1),
    dim_selectivity = round(dim_sel_all[cand_idx], 1),
    dim_ligandability = round(dim_lig, 1),
    dim_conservation = round(dim_con, 1),
    dim_clinical = round(dim_cli, 1),
    dim_adme = round(dim_adme, 1),
    n_dimensions_used = n_dims,
    stringsAsFactors = FALSE
  )

  tab <- tab[!is.na(tab$cpt_score), , drop = FALSE]
  # If enough true dependencies scored, rank CPT among deps only (strict score order).
  n_dep_scored <- sum(tab$is_dependency, na.rm = TRUE)
  if (n_dep_scored >= 10L) {
    tab <- tab[tab$is_dependency %in% TRUE, , drop = FALSE]
  }
  tab <- tab[order(-tab$cpt_score, tab$effect_size), , drop = FALSE]
  tab <- utils::head(tab, n)
  rownames(tab) <- NULL
  tab$cpt_rank <- seq_len(nrow(tab))

  records <- lapply(seq_len(nrow(tab)), function(i) as.list(tab[i, ]))

  list(
    subtype = subtype,
    dataset = dataset,
    n_cell_lines = bg$n_cell_lines,
    pool_size = length(cand_idx),
    pool_prefer_dependencies = length(dep_idx) >= 10L,
    require_ligandable = isTRUE(require_ligandable),
    n_returned = nrow(tab),
    caveats = c(
      "Ranked from a selective-dependency pool (prefer cancer_mean < -0.5), then re-scored with CPT dimensions.",
      "When >=10 scored dependencies exist, ranking is CPT score among dependencies only (descending).",
      paste("ADME dimension is computed where probe data exists but carries weight 0 by",
            "default, so it does not affect this ranking unless the weight is raised."),
      "Research prioritization aid only; not a clinical recommendation."
    ),
    genes = records
  )
}

# ---- Residue-resolved engagement, evidence tiers, Site CPT ------
# Gene-level CPT above is unchanged. These helpers annotate *engaged*
# cysteines (CR >= 4) and score the exact residue. Rankings never let a
# higher Site CPT promote a worse evidence tier.

if (!exists("cpt_gene_match_key", mode = "function")) {
  cpt_gene_match_key <- function(x) {
    toupper(trimws(sub("\\s*\\(\\d+\\)$", "", as.character(x))))
  }
}

#' Evidence tier for one cysteine (paper Tiers 1-4).
#'
#' Defined on the residue, independent of whether a probe is attached.
#' Engaged (CR >= 4) sites are the ones the paper ranks; untested sites
#' that are missing from the atlas are Tier 4.
#'
#' @param in_atlas Logical
#' @param functional Logical (atlas functional flag)
#' @param atlas_ligandable Logical (atlas ligandable flag)
#' @return Integer 1-4, same length as inputs
cpt_evidence_tier <- function(in_atlas, functional = FALSE, atlas_ligandable = FALSE) {
  # A zero-length argument makes the whole answer zero-length, the way
  # `TRUE & logical(0)` does. Taking the max instead let a scalar TRUE return
  # one tier for an empty table, and assigning that back as a column threw.
  lens <- c(length(in_atlas), length(functional), length(atlas_ligandable))
  if (any(lens == 0L)) return(integer(0))
  n <- max(lens)
  in_atlas <- rep(!is.na(in_atlas) & as.logical(in_atlas), length.out = n)
  functional <- rep(!is.na(functional) & as.logical(functional), length.out = n)
  atlas_ligandable <- rep(
    !is.na(atlas_ligandable) & as.logical(atlas_ligandable),
    length.out = n
  )
  out <- rep(4L, n)
  out[in_atlas & functional & atlas_ligandable] <- 1L
  out[in_atlas & functional & !atlas_ligandable] <- 2L
  out[in_atlas & !functional] <- 3L
  out
}

#' Annotate a probe table with exact-site evidence tiers (paper Tiers 1-4).
cpt_annotate_engaged_tiers <- function(probe_df, atlas) {
  if (is.null(probe_df) || !nrow(probe_df)) return(probe_df)
  if (!exists("cpt_cys_annotate_probe_table", mode = "function")) {
    return(probe_df)
  }
  ann <- cpt_cys_annotate_probe_table(probe_df, atlas)
  ann$evidence_tier <- cpt_evidence_tier(
    ann$cys_site_in_atlas,
    ann$cys_functional,
    ann$cys_ligandable
  )
  ann$evidence_tier_label <- cpt_evidence_tier_label(ann$evidence_tier)
  ann
}

cpt_evidence_tier_label <- function(tier) {
  vapply(as.integer(tier), function(t) {
    switch(as.character(t),
      "1" = "exact functional site with atlas ligandability",
      "2" = "exact functional site without atlas ligandability",
      "3" = "tested, not functional in Cys_editing contexts",
      "4" = "absent from atlas (untested)",
      "unknown"
    )
  }, character(1))
}

#' Compact CR >= min_cr probe-cysteine index for MCP / Site CPT.
#'
#' The full binding table is ~10.6M rows. Engagement at CR >= 4 is ~30k rows
#' and is what the paper treats as a chemoproteomic record.
cpt_build_smcl_index <- function(binding, min_cr = 4) {
  empty <- data.frame(
    gene_key = character(0),
    gene_name = character(0),
    proteinid = character(0),
    cysteineid = character(0),
    residue_number = integer(0),
    probe_name = character(0),
    CR = numeric(0),
    n_targets = numeric(0),
    ligandable = character(0),
    stringsAsFactors = FALSE
  )
  if (is.null(binding) || !nrow(binding)) return(empty)
  if (!all(c("probe_name", "CR", "cysteineid") %in% names(binding))) return(empty)

  cr <- suppressWarnings(as.numeric(binding$CR))
  keep <- !is.na(cr) & cr >= min_cr
  if (!any(keep)) return(empty)

  cols <- intersect(
    c("gene_name", "gene_name_key", "proteinid", "cysteineid",
      "probe_name", "CR", "n_targets", "ligandable"),
    names(binding)
  )
  x <- as.data.frame(binding[keep, cols, drop = FALSE], stringsAsFactors = FALSE)
  x$CR <- suppressWarnings(as.numeric(x$CR))
  if ("n_targets" %in% names(x)) {
    x$n_targets <- suppressWarnings(as.numeric(x$n_targets))
  } else {
    x$n_targets <- NA_real_
  }
  if ("gene_name" %in% names(x)) {
    x$gene_key <- cpt_gene_match_key(x$gene_name)
  } else if ("gene_name_key" %in% names(x)) {
    x$gene_key <- toupper(as.character(x$gene_name_key))
    x$gene_name <- x$gene_key
  } else {
    return(empty)
  }
  x$residue_number <- suppressWarnings(
    as.integer(sub(".*_C", "", as.character(x$cysteineid)))
  )
  x <- x[!is.na(x$residue_number) & nzchar(x$gene_key), , drop = FALSE]
  rownames(x) <- NULL
  x
}

cpt_parse_site_id <- function(site_id) {
  s <- as.character(site_id)[1]
  parts <- strsplit(s, "_", fixed = TRUE)[[1]]
  if (length(parts) < 2L) {
    return(list(gene = NA_character_, pos = NA_integer_))
  }
  pos <- suppressWarnings(as.integer(parts[length(parts)]))
  gene <- paste(parts[-length(parts)], collapse = "_")
  list(gene = gene, pos = pos)
}

cpt_smcl_rows_to_list <- function(rows, n_max = 15L) {
  if (is.null(rows) || !nrow(rows)) return(list())
  rows <- rows[order(-rows$CR, rows$n_targets, na.last = TRUE), , drop = FALSE]
  n_max <- as.integer(n_max)[1]
  if (!is.finite(n_max) || n_max < 1L) n_max <- 15L
  rows <- utils::head(rows, n_max)
  lapply(seq_len(nrow(rows)), function(i) {
    nt <- rows$n_targets[i]
    list(
      probe_name = as.character(rows$probe_name[i]),
      CR = round(rows$CR[i], 3),
      n_targets = if (is.na(nt)) NULL else as.integer(nt),
      cysteineid = as.character(rows$cysteineid[i]),
      proteinid = if ("proteinid" %in% names(rows)) as.character(rows$proteinid[i]) else NULL,
      prioritised = !is.na(nt) && nt <= 20
    )
  })
}

#' Site-level cysteine dimensions (this residue only, not max-across-gene).
cpt_atlas_site_dims <- function(row) {
  lig <- NA_real_
  con <- NA_real_
  cli <- NA_real_
  if (is.null(row) || !nrow(row)) {
    return(list(lig = lig, con = con, cli = cli))
  }
  if ("ligandability_score" %in% names(row) && !is.na(row$ligandability_score[1])) {
    lig <- as.numeric(row$ligandability_score[1])
  } else if ("ligandable" %in% names(row) && !is.na(row$ligandable[1])) {
    lig <- if (isTRUE(as.logical(row$ligandable[1]))) 50 else 0
  }
  if ("conservation_score" %in% names(row) && !is.na(row$conservation_score[1])) {
    con <- as.numeric(row$conservation_score[1])
    if (is.finite(con) && con <= 1) con <- con * 100
  }
  if ("clinvar_pathogenic" %in% names(row) && !all(is.na(row$clinvar_pathogenic))) {
    cli <- if (any(as.logical(row$clinvar_pathogenic), na.rm = TRUE)) 100 else 0
  }
  list(lig = lig, con = con, cli = cli)
}

#' Combined atlas + engaged-SMCL site records for one gene.
#'
#' @return NULL if the gene is in neither the atlas nor the SMCL index
cpt_gene_site_records <- function(gene, data_env, min_cr = 4,
                                  max_targets = NULL, smcl_cap = 15L) {
  gene_key <- cpt_gene_match_key(gene)[1]
  atlas <- data_env$cys_atlas
  smcl <- data_env$smcl_index

  atlas_hits <- data.frame()
  if (!is.null(atlas) && is.data.frame(atlas) && nrow(atlas) > 0) {
    atlas_hits <- atlas[toupper(as.character(atlas$gene_symbol)) == gene_key, , drop = FALSE]
  }
  smcl_hits <- if (!is.null(smcl) && is.data.frame(smcl) && nrow(smcl) > 0) {
    smcl[smcl$gene_key == gene_key, , drop = FALSE]
  } else {
    NULL
  }
  if (!is.null(smcl_hits) && nrow(smcl_hits) > 0 && !is.null(min_cr)) {
    smcl_hits <- smcl_hits[!is.na(smcl_hits$CR) & smcl_hits$CR >= min_cr, , drop = FALSE]
  }
  # Apply the selectivity ceiling on the FULL CR>=4 list, then cap for display.
  # Filtering after a CR-sorted cap would drop a site whose best probes are
  # promiscuous but which still has a restricted SMCL further down the list.
  if (!is.null(smcl_hits) && nrow(smcl_hits) > 0 &&
      !is.null(max_targets) && is.finite(as.numeric(max_targets)[1])) {
    mt <- as.numeric(max_targets)[1]
    keep_nt <- is.na(smcl_hits$n_targets) | smcl_hits$n_targets <= mt
    smcl_hits <- smcl_hits[keep_nt, , drop = FALSE]
  }
  if ((is.null(atlas_hits) || nrow(atlas_hits) == 0) &&
      (is.null(smcl_hits) || nrow(smcl_hits) == 0)) {
    return(NULL)
  }

  display_gene <- if (nrow(atlas_hits) > 0) {
    as.character(atlas_hits$gene_symbol[1])
  } else if (!is.null(smcl_hits) && nrow(smcl_hits) > 0 && "gene_name" %in% names(smcl_hits)) {
    as.character(smcl_hits$gene_name[1])
  } else {
    gene_key
  }
  uniprot <- if (nrow(atlas_hits) > 0) atlas_hits$uniprot_accession[1] else {
    if (!is.null(smcl_hits) && nrow(smcl_hits) > 0 && "proteinid" %in% names(smcl_hits)) {
      sub("-.*$", "", sub("_C.*$", "", as.character(smcl_hits$proteinid[1])))
    } else {
      NULL
    }
  }
  protein_name <- if (nrow(atlas_hits) > 0) atlas_hits$protein_name[1] else NULL

  atlas_pos <- if (nrow(atlas_hits) > 0) as.integer(atlas_hits$cysteine_position) else integer(0)
  smcl_pos <- if (!is.null(smcl_hits) && nrow(smcl_hits) > 0) {
    unique(as.integer(smcl_hits$residue_number))
  } else {
    integer(0)
  }
  smcl_pos <- smcl_pos[!is.na(smcl_pos)]
  all_pos <- sort(unique(c(atlas_pos, smcl_pos)))

  records <- lapply(all_pos, function(pos) {
    a_row <- NULL
    if (nrow(atlas_hits) > 0) {
      w <- which(as.integer(atlas_hits$cysteine_position) == pos)
      if (length(w)) a_row <- atlas_hits[w[1], , drop = FALSE]
    }
    in_atlas <- !is.null(a_row)
    p_rows <- data.frame()
    if (!is.null(smcl_hits) && nrow(smcl_hits) > 0) {
      p_rows <- smcl_hits[
        !is.na(smcl_hits$residue_number) & smcl_hits$residue_number == pos,
        ,
        drop = FALSE
      ]
    }
    engaged <- is.data.frame(p_rows) && nrow(p_rows) > 0
    functional <- if (in_atlas) as.logical(a_row$functional[1]) else NA
    ligandable <- if (in_atlas) as.logical(a_row$ligandable[1]) else NA
    fl <- if (in_atlas) as.logical(a_row$functional_ligandable[1]) else NA
    tier <- cpt_evidence_tier(in_atlas, isTRUE(functional), isTRUE(ligandable))
    dims <- cpt_atlas_site_dims(a_row)
    site_id <- if (in_atlas && !is.na(a_row$site_id[1]) && nzchar(a_row$site_id[1])) {
      as.character(a_row$site_id[1])
    } else {
      paste0(display_gene, "_", pos)
    }
    n_prior <- if (engaged && "n_targets" %in% names(p_rows)) {
      sum(!is.na(p_rows$n_targets) & p_rows$n_targets <= 20)
    } else {
      0L
    }
    list(
      site_id = site_id,
      gene_symbol = display_gene,
      cysteine_position = as.integer(pos),
      uniprot_accession = if (in_atlas) a_row$uniprot_accession[1] else uniprot,
      in_atlas = in_atlas,
      engaged = engaged,
      evidence_tier = as.integer(tier),
      evidence_tier_label = cpt_evidence_tier_label(tier),
      functional = if (in_atlas) as.logical(functional) else FALSE,
      ligandable = if (in_atlas) as.logical(ligandable) else FALSE,
      functional_ligandable = if (in_atlas) as.logical(fl) else FALSE,
      editor_support = if (in_atlas) a_row$editor_support[1] else NULL,
      study_context = if (in_atlas) a_row$study_context[1] else NULL,
      ligandability_score = if (in_atlas && !is.na(a_row$ligandability_score[1])) {
        round(a_row$ligandability_score[1], 1)
      } else {
        NULL
      },
      conservation_score = if (in_atlas && !is.na(a_row$conservation_score[1])) {
        round(a_row$conservation_score[1], 2)
      } else {
        NULL
      },
      clinvar_pathogenic = if (in_atlas) as.logical(a_row$clinvar_pathogenic[1]) else NA,
      n_smcls = as.integer(nrow(p_rows)),
      n_smcls_prioritised = as.integer(n_prior),
      smcls = cpt_smcl_rows_to_list(p_rows, smcl_cap),
      dim_ligandability = dims$lig,
      dim_conservation = dims$con,
      dim_clinical = dims$cli
    )
  })

  ord <- order(
    vapply(records, function(r) as.integer(r$evidence_tier), integer(1)),
    -vapply(records, function(r) as.integer(isTRUE(r$engaged)), integer(1)),
    -vapply(records, function(r) as.integer(isTRUE(r$functional_ligandable)), integer(1)),
    -vapply(records, function(r) as.integer(isTRUE(r$functional)), integer(1)),
    vapply(records, function(r) r$cysteine_position, integer(1))
  )
  records <- records[ord]

  list(
    gene = display_gene,
    protein_name = protein_name,
    uniprot_accession = uniprot,
    n_atlas = nrow(atlas_hits),
    n_engaged = sum(vapply(records, function(r) isTRUE(r$engaged), logical(1))),
    sites = records
  )
}

#' Site CPT for one residue: gene dep/sel percentiles + this site's cysteine dims.
#'
#' Does not modify gene-level CPT. Missing site dimensions are dropped and
#' remaining weights renormalised, same as the gene composite.
cpt_site_cpt_value <- function(site_rec, dim_dep, dim_sel, dim_adme = NA_real_,
                               weights = CPT_DEFAULT_WEIGHTS) {
  ds <- list(
    dependency_strength = dim_dep,
    cancer_selectivity = dim_sel,
    cysteine_ligandability = site_rec$dim_ligandability,
    conservation = site_rec$dim_conservation,
    clinical_evidence = site_rec$dim_clinical,
    adme_druggability = dim_adme
  )
  cpt_weighted_composite(ds, weights = weights, min_dimensions = 2L)
}

#' Attach Site CPT to gene site records using a precomputed gene score.
cpt_attach_site_cpt <- function(site_bundle, gene_score, weights = CPT_DEFAULT_WEIGHTS) {
  if (is.null(site_bundle) || !length(site_bundle$sites)) return(site_bundle)
  dim_dep <- gene_score$dimensions$dependency_strength$score
  dim_sel <- gene_score$dimensions$cancer_selectivity$score
  dim_adme <- if (!is.null(gene_score$dimensions$adme_druggability$score)) {
    gene_score$dimensions$adme_druggability$score
  } else {
    NA_real_
  }
  site_bundle$sites <- lapply(site_bundle$sites, function(s) {
    val <- cpt_site_cpt_value(s, dim_dep, dim_sel, dim_adme, weights)
    s$site_cpt_score <- if (is.na(val)) NA_real_ else round(val, 2)
    s
  })
  site_bundle
}

#' Restrict a gene list to the manuscript discovery filters.
#'
#' Paper RNAi defaults: effect size < -0.1, one-sided Welch p < 0.05,
#' whole-matrix mean gene effect > -0.5 (not a common essential). This is
#' NOT the CPT dependency pool (subtype mean < -0.5). Site CPT is still
#' computed with platform weights on the genes that pass.
cpt_apply_manuscript_gene_filters <- function(gene_rank,
                                              subtype,
                                              dataset,
                                              data_env,
                                              effect_size_max = NULL,
                                              p_max = NULL,
                                              exclude_common_essentials = FALSE,
                                              pool = 100L,
                                              weights = CPT_DEFAULT_WEIGHTS) {
  weights <- cpt_coerce_weights(weights)
  if (exists("cpt_normalize_dataset", mode = "function")) {
    dataset <- cpt_normalize_dataset(dataset)
  }
  matrix <- if (dataset == "CRISPR") data_env$crispr_matrix else data_env$rnai_matrix
  meta <- data_env$cancer_model_data
  cys_atlas <- data_env$cys_atlas
  if (is.null(matrix) || is.null(meta)) {
    return(gene_rank)
  }

  bg <- cpt_subtype_background(matrix, meta, subtype)
  gene_rank$n_cell_lines <- bg$n_cell_lines
  keep <- !is.na(bg$effect_sizes)
  es_cut <- if (is.null(effect_size_max)) NA_real_ else as.numeric(effect_size_max)[1]
  if (is.finite(es_cut)) {
    keep <- keep & bg$effect_sizes < es_cut
  }
  if (isTRUE(exclude_common_essentials)) {
    genome_means <- as.numeric(colMeans(matrix, na.rm = TRUE))
    keep <- keep & !is.na(genome_means) & genome_means > -0.5
  }
  idx <- which(keep)
  p_cut <- if (is.null(p_max)) NA_real_ else as.numeric(p_max)[1]
  pvals <- rep(NA_real_, length(bg$gene_names))
  if (is.finite(p_cut) && length(idx)) {
    meta_aligned <- meta[match(rownames(matrix), meta$ModelID), ]
    cancer_idx <- which(meta_aligned$OncotreeSubtype %in% subtype)
    other_idx <- which(!meta_aligned$OncotreeSubtype %in% subtype)
    pvals[idx] <- vapply(idx, function(j) {
      x <- matrix[cancer_idx, j]
      y <- matrix[other_idx, j]
      x <- x[!is.na(x)]
      y <- y[!is.na(y)]
      if (length(x) < 3L || length(y) < 3L) return(NA_real_)
      tryCatch(
        stats::t.test(x, y, alternative = "less")$p.value,
        error = function(e) NA_real_
      )
    }, numeric(1))
    idx <- idx[!is.na(pvals[idx]) & pvals[idx] < p_cut]
  }
  if (!length(idx)) {
    gene_rank$genes <- list()
    gene_rank$n_returned <- 0L
    gene_rank$pool_size <- 0L
    gene_rank$filter <- list(
      effect_size_max = es_cut,
      p_max = p_cut,
      exclude_common_essentials = isTRUE(exclude_common_essentials)
    )
    return(gene_rank)
  }

  idx <- idx[order(bg$effect_sizes[idx])]
  pool <- as.integer(max(1L, min(as.integer(pool)[1], 500L)))
  if (is.na(pool)) pool <- 100L
  idx <- idx[seq_len(min(pool, length(idx)))]

  dim_dep_all <- score_dependency_strength(bg$cancer_means)
  dim_sel_all <- score_cancer_selectivity(bg$effect_sizes)
  cand_genes <- bg$gene_names[idx]
  dim_lig <- score_cysteine_ligandability(cand_genes, cys_atlas)
  dim_con <- score_conservation(cand_genes, cys_atlas)
  dim_cli <- score_clinical_evidence(cand_genes, cys_atlas)
  dim_adme <- score_adme_druggability(cand_genes, data_env$adme_gene_scores)

  genes <- lapply(seq_along(idx), function(i) {
    j <- idx[i]
    ds <- list(
      dependency_strength = dim_dep_all[j],
      cancer_selectivity = dim_sel_all[j],
      cysteine_ligandability = dim_lig[i],
      conservation = dim_con[i],
      clinical_evidence = dim_cli[i],
      adme_druggability = dim_adme[i]
    )
    sc <- cpt_weighted_composite(ds, weights = weights, min_dimensions = 2L)
    list(
      gene = cand_genes[i],
      cpt_score = if (is.na(sc)) NA_real_ else round(sc, 2),
      cancer_mean = round(bg$cancer_means[j], 4),
      effect_size = round(bg$effect_sizes[j], 4),
      p_value = if (is.na(pvals[j])) NULL else signif(pvals[j], 4),
      dim_dependency = round(dim_dep_all[j], 1),
      dim_selectivity = round(dim_sel_all[j], 1),
      dim_adme = round(dim_adme[i], 1)
    )
  })
  gene_rank$genes <- genes
  gene_rank$n_returned <- length(genes)
  gene_rank$pool_size <- length(genes)
  gene_rank$pool_prefer_dependencies <- FALSE
  gene_rank$filter <- list(
    effect_size_max = es_cut,
    p_max = p_cut,
    exclude_common_essentials = isTRUE(exclude_common_essentials)
  )
  gene_rank
}

#' Rank engaged sites in a subtype. Order is evidence tier (1 first), then Site CPT.
#'
#' Gene-level `cpt_rank_targets` is unchanged. A high Site CPT on a Tier 4
#' residue cannot rank above a Tier 1/2 site.
cpt_rank_site_targets <- function(subtype,
                                  dataset = "CRISPR",
                                  n = 20L,
                                  pool = 100L,
                                  min_cr = 4,
                                  max_targets = 20,
                                  effect_size_max = NULL,
                                  p_max = NULL,
                                  exclude_common_essentials = FALSE,
                                  data_env,
                                  weights = CPT_DEFAULT_WEIGHTS) {
  weights <- cpt_coerce_weights(weights)
  n <- as.integer(max(1L, min(as.integer(n)[1], 100L)))
  if (is.na(n) || n < 1L) n <- 20L
  max_targets <- if (is.null(max_targets)) NULL else as.numeric(max_targets)[1]
  manuscript_pool <- (!is.null(effect_size_max) && is.finite(as.numeric(effect_size_max)[1])) ||
    (!is.null(p_max) && is.finite(as.numeric(p_max)[1])) ||
    isTRUE(exclude_common_essentials)

  if (isTRUE(manuscript_pool)) {
    gene_rank <- cpt_apply_manuscript_gene_filters(
      list(subtype = subtype, dataset = dataset, genes = list()),
      subtype = subtype,
      dataset = dataset,
      data_env = data_env,
      effect_size_max = effect_size_max,
      p_max = p_max,
      exclude_common_essentials = exclude_common_essentials,
      pool = pool,
      weights = weights
    )
  } else {
    gene_rank <- cpt_rank_targets(
      subtype = subtype,
      dataset = dataset,
      n = pool,
      pool = pool,
      require_ligandable = FALSE,
      data_env = data_env,
      weights = weights
    )
  }

  rows <- list()
  k <- 1L
  for (g in gene_rank$genes) {
    gene <- g$gene
    bundle <- cpt_gene_site_records(
      gene, data_env,
      min_cr = min_cr,
      max_targets = max_targets,
      smcl_cap = 8L
    )
    if (is.null(bundle)) next
    dim_dep <- g$dim_dependency
    dim_sel <- g$dim_selectivity
    dim_adme <- g$dim_adme
    for (s in bundle$sites) {
      if (!isTRUE(s$engaged)) next
      site_cpt <- cpt_site_cpt_value(s, dim_dep, dim_sel, dim_adme, weights)
      top <- if (length(s$smcls)) s$smcls[[1]] else NULL
      rows[[k]] <- data.frame(
        gene = gene,
        site_id = s$site_id,
        cysteine_position = s$cysteine_position,
        evidence_tier = s$evidence_tier,
        evidence_tier_label = s$evidence_tier_label,
        site_cpt_score = if (is.na(site_cpt)) NA_real_ else round(site_cpt, 2),
        gene_cpt_score = g$cpt_score,
        cancer_mean = g$cancer_mean,
        effect_size = g$effect_size,
        n_smcls = s$n_smcls,
        top_probe = if (is.null(top)) NA_character_ else top$probe_name,
        top_CR = if (is.null(top)) NA_real_ else top$CR,
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }

  if (!length(rows)) {
    return(list(
      subtype = subtype,
      dataset = dataset,
      n_cell_lines = gene_rank$n_cell_lines,
      pool_size = gene_rank$pool_size,
      min_cr = min_cr,
      max_targets = max_targets,
      n_returned = 0L,
      ranking = "evidence_tier then site_cpt_score within tier",
      sites = list(),
      caveats = c(
        "No engaged cysteines (CR >= min_cr) in the ranked gene pool.",
        "Site CPT never crosses evidence-tier boundaries."
      )
    ))
  }

  tab <- do.call(rbind, rows)
  tab <- tab[!is.na(tab$site_cpt_score), , drop = FALSE]
  tab$rank_in_tier <- NA_integer_
  for (t in sort(unique(tab$evidence_tier))) {
    idx <- which(tab$evidence_tier == t)
    o <- order(-tab$site_cpt_score[idx], tab$effect_size[idx])
    tab$rank_in_tier[idx[o]] <- seq_along(idx)
  }
  tab <- tab[order(tab$evidence_tier, -tab$site_cpt_score, tab$effect_size), , drop = FALSE]
  tab <- utils::head(tab, n)
  rownames(tab) <- NULL
  tab$rank <- seq_len(nrow(tab))
  records <- lapply(seq_len(nrow(tab)), function(i) as.list(tab[i, ]))

  list(
    subtype = subtype,
    dataset = dataset,
    n_cell_lines = gene_rank$n_cell_lines,
    pool_size = gene_rank$pool_size,
    min_cr = min_cr,
    max_targets = max_targets,
    effect_size_max = if (is.null(effect_size_max)) NULL else as.numeric(effect_size_max)[1],
    p_max = if (is.null(p_max)) NULL else as.numeric(p_max)[1],
    exclude_common_essentials = isTRUE(exclude_common_essentials),
    gene_pool = if (isTRUE(manuscript_pool)) "manuscript_filters" else "cpt_dependency_pool",
    n_returned = nrow(tab),
    ranking = "evidence_tier (1 first) then site_cpt_score within tier",
    caveats = c(
      "Site CPT uses gene-level dependency/selectivity percentiles plus cysteine dimensions from the exact residue.",
      "Candidates are ordered by evidence tier first, so a high-scoring untested (Tier 4) site cannot outrank a functional Tier 1 or 2 site.",
      "Engagement requires CR >= min_cr; max_targets is applied to all engaging probes before the display cap.",
      "Default gene pool is the CPT dependency ranking. Pass effect_size_max / p_max / exclude_common_essentials to use the manuscript RNAi filters (ES, Welch p, whole-matrix mean > -0.5).",
      "Engagement is from the bundled chemoproteomic table (CysDB-indexed studies), not a live CysDB query.",
      "Research prioritization aid only; not a clinical recommendation."
    ),
    sites = records
  )
}

