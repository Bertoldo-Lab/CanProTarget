# ============================================================
# Script:   functions.R
# Purpose:  Pure (non-Shiny) helper functions for CanProTarget:
#           - Gene effect statistics (DepMap-style, uses cdsrmodels)
#           - Probe analysis (chemoproteomic linking, CR + targets)
#           - SwissADME data processing
# Inputs:   Called from app.R and Shiny modules via shared_data
# Outputs:  Analysis data frames; no Shiny side-effects
# ============================================================
# Sections:
#   1. Gene effect analysis
#        1.0 feasible_oncotree_subtypes() — subtypes with enough lines in matrix
#        1.1 create_group_membership()  — binary cancer label vector
#        1.2 lin_ass_pval()             — linear association p-values
#        1.3 run_lm_ge()                — limma moderated t-test
#        1.4 ge_analysis()              — combine + filter results
#   2. (Removed: probe_analysis superseded by preprocessed RDS pipeline)
#   3. SwissADME analysis
#        3.1 process_swissadme_data()         — swissadme_preprocessed.rds only
#        3.2 filter_swissadme_data()          — filter to one probe
#        3.3 overlay_canonical_smiles_from_swissadme() — align lookup SMILES with SwissADME
# ============================================================


# ---- 1. Gene effect analysis -----------------------------------

# 1.0 Subtype feasibility (DepMap matrix ∩ metadata)
# Count cell lines per OncotreeSubtype among models that appear in the gene-effect matrix.
n_cell_lines_by_oncotree_subtype <- function(meta, model_ids_in_matrix) {
  meta <- meta[!is.na(meta$ModelID) & meta$ModelID %in% model_ids_in_matrix, , drop = FALSE]
  if (!nrow(meta)) {
    return(setNames(integer(0), character(0)))
  }
  tab <- table(meta$OncotreeSubtype, useNA = "no")
  sort(tab, decreasing = TRUE)
}

# Minimum cell lines per group for linear GE (picker + validation); cdsrmodels n.min aligned.
MIN_GE_CELL_LINES <- 3L

# Subtypes with at least min_lines cell lines in the matrix (aligned with lin_associations n.min).
feasible_oncotree_subtypes <- function(meta, model_ids_in_matrix, min_lines = MIN_GE_CELL_LINES) {
  tab <- n_cell_lines_by_oncotree_subtype(meta, model_ids_in_matrix)
  if (!length(tab)) {
    return(character(0))
  }
  names(tab)[tab >= min_lines]
}

# Before on-the-fly lin_ass_pval / limma: ensure binary Y has enough 1s and 0s.
validate_linear_ge_sample_size <- function(effect_size_matrix, cell_line_info, selected_subtype,
                                          min_in_subtype = MIN_GE_CELL_LINES,
                                          min_other = MIN_GE_CELL_LINES) {
  Y <- create_group_membership(cell_line_info, selected_subtype, effect_size_matrix)
  n1 <- sum(Y == 1, na.rm = TRUE)
  n0 <- sum(Y == 0, na.rm = TRUE)
  if (n1 < min_in_subtype) {
    stop(
      "Only ", n1, " cell line(s) in the gene effect matrix match the selected cancer subtype(s). ",
      "At least ", min_in_subtype, " are required (overlap between matrix rows and DepMap ModelID / OncotreeSubtype). ",
      "Choose another subtype from the list.",
      call. = FALSE
    )
  }
  if (n0 < min_other) {
    stop(
      "Only ", n0, " cell line(s) fall outside the selected subtype(s) in this matrix; ",
      "at least ", min_other, " are required for the comparison. Check matrix vs cancer_model_data.rds alignment.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# 1.1 create_group_membership
# Builds a binary matrix: 1 = selected cancer subtype, 0 = everything else.
# Used as the response vector Y for both lin_ass_pval and run_lm_ge.
create_group_membership <- function(cell_line_info, selected_subtype, effect_size_matrix) {
  # Order cell line info to match the rows of the effect size matrix
  cell_line_info <- cell_line_info[match(rownames(effect_size_matrix), cell_line_info$ModelID), ]

  # Binary label: 1 = cancer of interest, 0 = all other cell lines
  group_membership <- ifelse(cell_line_info$OncotreeSubtype %in% selected_subtype, 1, 0)
  Y <- matrix(group_membership, ncol = 1)
  Y
}

# 1.2 lin_ass_pval
# Computes per-gene linear-association p-values using cdsrmodels::lin_associations.
# Returns the result table from lin_associations (one row per gene).
lin_ass_pval <- function(effect_size_matrix, cell_line_info, selected_subtype) {
  Y <- create_group_membership(cell_line_info, selected_subtype, effect_size_matrix)

  results <- lin_associations(
    X             = as.matrix(effect_size_matrix),
    Y             = Y,
    W             = NULL,
    n.min         = MIN_GE_CELL_LINES,
    shrinkage     = TRUE,
    a             = 0,
    MHC_direction = NULL
  )
  results$res.table
}

# 1.3 run_lm_ge
# Runs limma-moderated t-test (cancer of interest vs all others) per gene.
# Returns a data frame with EffectSize (logFC), Avg, adj.P.Val, etc.
run_lm_ge <- function(effect_size_matrix, cell_line_info, selected_subtype) {
  Y <- create_group_membership(cell_line_info, selected_subtype, effect_size_matrix)

  # run_lm_stats_limma from cdsrmodels: genes in columns, cell lines in rows
  run_lm_stats_limma(
    mat         = as.matrix(effect_size_matrix),
    vec         = Y,
    covars      = NULL,
    weights     = NULL,
    target_type = "feature",
    limma_trend = FALSE
  )
}

# ---- 1.35 Precomputed index + static cancer subtype list files ----------------
# Subtype dropdown reads plain text under data/ (see app_config$cancer_subtype_list_files).
# Offline: docs/scripts/preprocess_data.R calls cancer_subtype_choices_build() and writeLines().

dataset_short_name <- function(label) {
  if (grepl("CRISPR", label, ignore.case = TRUE)) {
    "CRISPR"
  } else if (grepl("RNAi", label, ignore.case = TRUE)) {
    "RNAi"
  } else {
    gsub("[^A-Za-z0-9]+", "", label)
  }
}

# Index of precomputed per-subtype effect sizes: TSV filenames under data/ (no CSV logs).
precomputed_index_dataframe <- function(data_dir) {
  parts <- list()
  for (subdir in "precomputed_effectsizes") {
    d <- file.path(data_dir, subdir)
    if (!dir.exists(d)) {
      next
    }
    files <- list.files(d, pattern = "\\.(rds|tsv)$", full.names = FALSE)
    if (!length(files)) {
      next
    }
    parts[[length(parts) + 1L]] <- data.frame(
      dataset = ifelse(grepl("^CRISPR_", files), "CRISPR", "RNAi"),
      subtype = sub("^(CRISPR|RNAi)_", "", sub("\\.(rds|tsv)$", "", files)),
      n_lines = NA_integer_,
      n_genes = NA_integer_,
      status  = "SUCCESS",
      time_sec = NA_real_,
      stringsAsFactors = FALSE
    )
  }
  if (!length(parts)) {
    return(NULL)
  }
  out <- dplyr::bind_rows(parts)
  dplyr::distinct(out, .data$dataset, .data$subtype, .keep_all = TRUE)
}

#' Normalise a gene label for matching ACROSS data sources.
#'
#' DepMap gene-effect matrix columns carry an Entrez suffix (`"KRAS (3845)"`),
#' which flows into every `gene_name` produced by `ge_analysis()`. The
#' chemoproteomics binding table and the cysteine atlas use bare symbols
#' (`"KRAS"`). Comparing the two raw strings silently matches nothing, so any
#' join between a dependency table and a probe/atlas table must go through here.
#'
#' @param x Character vector of gene labels, with or without the Entrez suffix
#' @return Upper-case bare symbols
cpt_gene_match_key <- function(x) {
  toupper(trimws(sub("\\s*\\(\\d+\\)$", "", as.character(x))))
}

cpt_sanitize_subtype_key <- function(subtype) {
  out <- gsub("[^A-Za-z0-9]+", "_", subtype)
  out <- gsub("_+", "_", out)
  out <- gsub("^_|_$", "", out)
  out
}

# Alias used by dependencies_module.R and precompute scripts
sanitize_subtype <- cpt_sanitize_subtype_key

#' Build ordered cancer subtype choices (same rules as the legacy Shiny picker).
#' Used offline in docs/scripts/preprocess_data.R to write cancer_subtypes_*.txt.
cancer_subtype_choices_build <- function(meta, model_ids, data_dir, dataset_short) {
  meta_subtypes <- sort(unique(meta$OncotreeSubtype))
  meta_subtypes <- meta_subtypes[!is.na(meta_subtypes) & meta_subtypes != ""]

  feasible <- if (!is.null(model_ids) && length(model_ids)) {
    feasible_oncotree_subtypes(meta, model_ids, min_lines = MIN_GE_CELL_LINES)
  } else {
    character(0)
  }

  # Every statistically feasible subtype is selectable. A precomputed TSV is
  # an acceleration path, not an availability requirement: ge_results()
  # explicitly falls back to on-the-fly analysis when no TSV exists.
  if (length(feasible)) {
    choices <- feasible
  } else if (is.null(model_ids) || !length(model_ids)) {
    choices <- meta_subtypes
  } else {
    choices <- character(0)
  }

  sort(unique(choices))
}

read_cancer_subtype_list_file <- function(path) {
  if (!file.exists(path)) {
    return(character(0))
  }
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  lines <- trimws(lines)
  lines <- lines[nzchar(lines) & !grepl("^#", lines)]
  sort(unique(lines))
}

# 1.4 ge_analysis
# Combines lin_ass_pval and run_lm_ge results, computes per-group averages,
# runs one-sided t-tests, and applies user-selected filters to produce:
#   - all_gene_ge_df : every gene with all statistics
#   - cancer_gene_df : filtered to cancer-specific dependencies only
ge_analysis <- function(pval_results,
                        lm_results,
                        apply_pvalue_filter       = FALSE,
                        min_effect_size,
                        exclude_common_essentials = FALSE,
                        filter_pval_vs_noncancer  = FALSE,
                        effect_size_matrix        = NULL,
                        cancer_model_data         = NULL,
                        selected_subtype          = NULL) {

  # Standardise column names from upstream functions
  lm_results   <- dplyr::rename(lm_results, gene_name = feature) %>% as.data.frame()
  pval_results <- dplyr::select(pval_results, ind.var, p.val) %>%
    dplyr::rename(gene_name = ind.var, p_value = p.val) %>%
    as.data.frame()

  # Merge limma results with linear-association p-values
  all_gene_ge_df <- merge(lm_results, pval_results, by = "gene_name", all.x = TRUE)

  # Add -log10(p_value) for volcano plots
  all_gene_ge_df <- all_gene_ge_df %>%
    dplyr::mutate(
      neg_log10_p_value = ifelse(
        is.na(p_value) | !is.numeric(p_value), NA, -log10(p_value)
      )
    )

  # Per-group averages and t-tests
  if (!is.null(effect_size_matrix) &&
      !is.null(cancer_model_data)  &&
      !is.null(selected_subtype)) {

    model_data <- cancer_model_data[
      match(rownames(effect_size_matrix), cancer_model_data$ModelID), ]

    cancer_indices       <- which(model_data$OncotreeSubtype %in% selected_subtype)
    non_cancer_indices   <- which(model_data$OncotreePrimaryDisease == "Non-Cancerous")
    other_indices        <- which(!(model_data$OncotreeSubtype %in% selected_subtype))
    other_cancer_indices <- which(
      !(model_data$OncotreeSubtype %in% selected_subtype) &
        model_data$OncotreePrimaryDisease != "Non-Cancerous"
    )

    # Other_Avg: mean across all cell lines EXCEPT the cancer of interest
    if (length(other_indices) > 0) {
      other_avg <- colMeans(effect_size_matrix[other_indices, , drop = FALSE], na.rm = TRUE)
      all_gene_ge_df <- all_gene_ge_df %>%
        dplyr::left_join(
          data.frame(gene_name = names(other_avg),
                     Other_Avg = as.numeric(other_avg),
                     stringsAsFactors = FALSE),
          by = "gene_name"
        )
    } else { all_gene_ge_df$Other_Avg <- NA }

    # NonCancer_Avg: mean across non-cancerous cell lines
    if (length(non_cancer_indices) > 0) {
      nc_avg <- colMeans(effect_size_matrix[non_cancer_indices, , drop = FALSE], na.rm = TRUE)
      all_gene_ge_df <- all_gene_ge_df %>%
        dplyr::left_join(
          data.frame(gene_name     = names(nc_avg),
                     NonCancer_Avg = as.numeric(nc_avg),
                     stringsAsFactors = FALSE),
          by = "gene_name"
        )
    } else { all_gene_ge_df$NonCancer_Avg <- NA }

    # Cancer_Avg: mean across the cancer of interest cell lines
    if (length(cancer_indices) > 0) {
      c_avg <- colMeans(effect_size_matrix[cancer_indices, , drop = FALSE], na.rm = TRUE)
      all_gene_ge_df <- all_gene_ge_df %>%
        dplyr::left_join(
          data.frame(gene_name  = names(c_avg),
                     Cancer_Avg = as.numeric(c_avg),
                     stringsAsFactors = FALSE),
          by = "gene_name"
        )
    } else { all_gene_ge_df$Cancer_Avg <- NA }

    # pval_vs_Other_Avg: one-sided t-test, cancer vs ALL other cell lines
    if (length(cancer_indices) >= 3 && length(other_indices) >= 3) {
      pvals_vs_other <- sapply(colnames(effect_size_matrix), function(gene) {
        x <- effect_size_matrix[cancer_indices,  gene]
        y <- effect_size_matrix[other_indices,   gene]
        x <- x[!is.na(x)]; y <- y[!is.na(y)]
        if (length(x) < 3 || length(y) < 3) return(NA_real_)
        tryCatch(t.test(x, y, alternative = "less")$p.value, error = function(e) NA_real_)
      })
      all_gene_ge_df <- all_gene_ge_df %>%
        dplyr::left_join(
          data.frame(gene_name         = names(pvals_vs_other),
                     pval_vs_Other_Avg = as.numeric(pvals_vs_other),
                     stringsAsFactors  = FALSE),
          by = "gene_name"
        )
    } else { all_gene_ge_df$pval_vs_Other_Avg <- NA }

    # pval_vs_NonCancer: one-sided t-test, cancer vs non-cancerous cell lines
    if (length(cancer_indices) >= 3 && length(non_cancer_indices) >= 3) {
      pvals_vs_nc <- sapply(colnames(effect_size_matrix), function(gene) {
        x <- effect_size_matrix[cancer_indices,     gene]
        y <- effect_size_matrix[non_cancer_indices, gene]
        x <- x[!is.na(x)]; y <- y[!is.na(y)]
        if (length(x) < 3 || length(y) < 3) return(NA_real_)
        tryCatch(t.test(x, y, alternative = "less")$p.value, error = function(e) NA_real_)
      })
      all_gene_ge_df <- all_gene_ge_df %>%
        dplyr::left_join(
          data.frame(gene_name         = names(pvals_vs_nc),
                     pval_vs_NonCancer = as.numeric(pvals_vs_nc),
                     stringsAsFactors  = FALSE),
          by = "gene_name"
        )
    } else { all_gene_ge_df$pval_vs_NonCancer <- NA }

    # pval_vs_OtherCancers: one-sided t-test, cancer of interest vs other cancers
    if (length(cancer_indices) >= 3 && length(other_cancer_indices) >= 3) {
      pvals_vs_oc <- sapply(colnames(effect_size_matrix), function(gene) {
        x <- effect_size_matrix[cancer_indices,       gene]
        y <- effect_size_matrix[other_cancer_indices, gene]
        x <- x[!is.na(x)]; y <- y[!is.na(y)]
        if (length(x) < 3 || length(y) < 3) return(NA_real_)
        tryCatch(t.test(x, y, alternative = "less")$p.value, error = function(e) NA_real_)
      })
      all_gene_ge_df <- all_gene_ge_df %>%
        dplyr::left_join(
          data.frame(gene_name            = names(pvals_vs_oc),
                     pval_vs_OtherCancers = as.numeric(pvals_vs_oc),
                     stringsAsFactors     = FALSE),
          by = "gene_name"
        )
    } else { all_gene_ge_df$pval_vs_OtherCancers <- NA }

  } else {
    all_gene_ge_df$Other_Avg           <- NA
    all_gene_ge_df$NonCancer_Avg       <- NA
    all_gene_ge_df$Cancer_Avg          <- NA
    all_gene_ge_df$pval_vs_Other_Avg   <- NA
    all_gene_ge_df$pval_vs_NonCancer   <- NA
    all_gene_ge_df$pval_vs_OtherCancers <- NA
  }

  # ---- Apply filters to produce cancer_gene_df ----

  # Core filter: effect size threshold (more negative = stronger dependency)
  cancer_gene_df <- all_gene_ge_df %>%
    dplyr::filter(!is.na(EffectSize) & EffectSize <= min_effect_size)

  # Optional: linear-association p-value < 0.05
  if (apply_pvalue_filter) {
    cancer_gene_df <- dplyr::filter(cancer_gene_df, p_value < 0.05)
  }

  # Optional: exclude broadly essential genes (Avg < -0.5 across all cell lines)
  if (exclude_common_essentials) {
    essential_genes <- cancer_gene_df %>%
      dplyr::filter(Avg < -0.5) %>%
      dplyr::pull(gene_name)
    cancer_gene_df <- dplyr::filter(cancer_gene_df, !gene_name %in% essential_genes)
  }

  # Optional: significantly more essential vs non-cancer (one-sided t-test p < 0.05)
  if (filter_pval_vs_noncancer) {
    cancer_gene_df <- cancer_gene_df %>%
      dplyr::filter(!is.na(pval_vs_NonCancer) & pval_vs_NonCancer < 0.05)
  }

  list(all_gene_ge_df = all_gene_ge_df, cancer_gene_df = cancer_gene_df)
}

# 1.5 Group comparison data (pure; used by Dependencies tab + tests)
#' Build long-format gene-effect scores for cancer / non-cancer / other groups.
#' @param mat Gene-effect matrix (rows = ModelID, cols = genes)
#' @param meta Cancer model metadata with ModelID, OncotreeSubtype,
#'   OncotreePrimaryDisease, StrippedCellLineName
#' @param gene_id Column name in mat
#' @param cancer_subtypes Character vector of OncotreeSubtype labels for "of interest"
#' @return data.frame with score, group, cell_line; or NULL if invalid / empty
cpt_group_comparison_df <- function(mat, meta, gene_id, cancer_subtypes) {
  if (is.null(gene_id) || !nzchar(as.character(gene_id)[1])) return(NULL)
  gene_id <- as.character(gene_id)[1]
  if (is.null(mat) || is.null(meta) || !(gene_id %in% colnames(mat))) return(NULL)
  if (!all(c("ModelID", "OncotreeSubtype", "OncotreePrimaryDisease") %in% colnames(meta))) {
    return(NULL)
  }

  meta_aligned <- meta[match(rownames(mat), meta$ModelID), ]

  cancer_idx    <- which(meta_aligned$OncotreeSubtype %in% cancer_subtypes)
  noncancer_idx <- which(meta_aligned$OncotreePrimaryDisease == "Non-Cancerous")
  other_idx     <- which(
    !(meta_aligned$OncotreeSubtype %in% cancer_subtypes) &
      meta_aligned$OncotreePrimaryDisease != "Non-Cancerous"
  )

  cell_names <- if ("StrippedCellLineName" %in% colnames(meta_aligned)) {
    as.character(meta_aligned$StrippedCellLineName)
  } else {
    rownames(mat)
  }

  make_group_df <- function(idx, label) {
    if (!length(idx)) return(NULL)
    data.frame(
      score     = as.numeric(mat[idx, gene_id]),
      group     = label,
      cell_line = cell_names[idx],
      stringsAsFactors = FALSE
    )
  }

  gdf <- dplyr::bind_rows(
    make_group_df(cancer_idx,    "Cancer of Interest"),
    make_group_df(noncancer_idx, "Non-Cancer"),
    make_group_df(other_idx,     "Other Cancers")
  )
  if (is.null(gdf) || !nrow(gdf)) return(NULL)
  gdf
}


# ---- 2. Probe analysis -----------------------------------------
# NOTE: Legacy probe_analysis() removed (never called; superseded by
# preprocessed protein_binding_lookup_preprocessed.rds pipeline).
# Probe data is now handled directly in dependencies_module.R using
# the pre-built RDS file. See docs/DATA_PROVENANCE.md.


# ---- 3. SwissADME analysis -------------------------------------

#' Correct the spelling of the Kuljanin/Gygi study in the Dataset column.
#'
#' The source chemoproteomics table spells it "kunljanin". The value is shown
#' to users in the probe tables and the study is named in DATA_LICENSES.md, so
#' it is corrected on load rather than left to the file. Dataset is only ever
#' displayed, grouped or sorted, never matched against a literal, so renaming
#' the level changes nothing else.
#' @param x character or factor Dataset column
#' @return the same type, with the study name spelled correctly
cpt_fix_dataset_names <- function(x) {
  fix <- function(v) sub("^kunljanin_", "kuljanin_", as.character(v))
  if (is.factor(x)) {
    levels(x) <- fix(levels(x))
    x
  } else {
    fix(x)
  }
}

# One spelling of the rule, defined with the scorer that also has to apply it.
simplify_probe_name <- function(x) cpt_canonical_probe_name(x)

# 3.0 BOILED-Egg ggplot (pure; SwissADME geometry: x = WLOGP, y = TPSA)
#' Build a BOILED-Egg plot from a SwissADME-like table.
#' @param adme data.frame with probe_name and WLOGP/TPSA (or aliases)
#' @param selected_probe Character probe to highlight
#' @param show_mode "all" or "selected_only"
#' @param theme_fn Optional ggplot2 theme function (default theme_minimal)
#' @return ggplot object
cpt_boiled_egg_gg <- function(adme, selected_probe = NULL, show_mode = "all",
                              theme_fn = NULL) {
  if (is.null(theme_fn)) {
    theme_fn <- function() ggplot2::theme_minimal()
  }
  if (is.null(adme) || !nrow(adme) || !("probe_name" %in% colnames(adme))) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No SwissADME data") +
        ggplot2::theme_void()
    )
  }

  logp_col <- intersect(
    c("WLOGP", "wlogp", "Consensus.Log.Po.w", "consensus_logp", "XLogP3"),
    colnames(adme)
  )[1]
  tpsa_col <- intersect(
    c("TPSA", "tpsa", "Topological.Polar.Surface.Area"),
    colnames(adme)
  )[1]

  if (is.na(logp_col) || is.na(tpsa_col)) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate(
          "text", x = 0.5, y = 0.5,
          label = "WLOGP or TPSA column not found in data"
        ) +
        ggplot2::theme_void()
    )
  }

  df <- data.frame(
    probe = as.character(adme$probe_name),
    WLOGP = suppressWarnings(as.numeric(adme[[logp_col]])),
    TPSA  = suppressWarnings(as.numeric(adme[[tpsa_col]])),
    stringsAsFactors = FALSE
  )
  df$selected <- !is.null(selected_probe) &
    nzchar(as.character(selected_probe)[1]) &
    df$probe == as.character(selected_probe)[1]

  if (is.null(show_mode) || !nzchar(show_mode)) show_mode <- "all"

  # SwissADME BOILED-Egg geometry: x = WLOGP, y = TPSA
  theta <- seq(0, 2 * pi, length.out = 200)
  egg_outer <- data.frame(
    WLOGP = 2.5 + 3 * cos(theta),
    TPSA  = 90  + 90 * sin(theta)
  )
  yolk <- data.frame(
    WLOGP = 2.5 + 2 * cos(theta),
    TPSA  = 60  + 60 * sin(theta)
  )

  g <- ggplot2::ggplot() +
    ggplot2::geom_polygon(
      data = egg_outer,
      ggplot2::aes(x = WLOGP, y = TPSA),
      fill = "grey95", colour = "grey70"
    ) +
    ggplot2::geom_polygon(
      data = yolk,
      ggplot2::aes(x = WLOGP, y = TPSA),
      fill = "#F5D76E", colour = NA, alpha = 0.7
    ) +
    ggplot2::coord_cartesian(xlim = c(-3, 7), ylim = c(0, 180), expand = FALSE) +
    ggplot2::labs(
      x = "WLOGP",
      y = "TPSA (Å²)",
      # No in-plot title: the section heading directly above already names
      # the chart, and a second one clipped against the panel edge.
      title = NULL,
      # White = high GI absorption (HIA); yolk = BBB permeation (Daina & Zoete)
      subtitle = "White: high HIA predicted · Yolk: BBB permeation predicted"
    ) +
    theme_fn()

  if (identical(show_mode, "selected_only")) {
    sel <- df[df$selected, , drop = FALSE]
    if (!nrow(sel)) {
      g <- g +
        ggplot2::annotate(
          "text",
          x = 2, y = 90,
          label = "No row for the selected probe in SwissADME table",
          colour = "grey40", size = 3.5
        )
    } else {
      g <- g +
        # `text` is for ggplotly tooltips (ggplot2 warns, plotly uses it)
        suppressWarnings(
          ggplot2::geom_point(
            data = sel,
            ggplot2::aes(x = WLOGP, y = TPSA, text = probe),
            colour = CPT_PAL$ink, size = 3.4
          )
        )
    }
  } else {
    g <- g +
      suppressWarnings(
        ggplot2::geom_point(
          data = df,
          ggplot2::aes(x = WLOGP, y = TPSA, text = probe),
          colour = CPT_PAL$light, alpha = 0.6, size = 2
        )
      )
    sel <- df[df$selected, , drop = FALSE]
    if (nrow(sel)) {
      g <- g +
        suppressWarnings(
          ggplot2::geom_point(
            data = sel,
            ggplot2::aes(x = WLOGP, y = TPSA, text = probe),
            colour = CPT_PAL$ink, size = 3.4
          )
        )
    }
  }
  g
}

# 3.1 process_swissadme_data (swissadme_preprocessed.rds only; build offline — see docs/DATA_PROVENANCE.md)
process_swissadme_data <- function(data_dir) {
  rds_path <- file.path(data_dir, "swissadme_preprocessed.rds")
  if (!file.exists(rds_path)) {
    return(tibble::tibble())
  }
  adme <- tryCatch(
    readRDS(rds_path),
    error = function(e) {
      stop("Failed to read swissadme_preprocessed.rds: ", conditionMessage(e),
           call. = FALSE)
    }
  )
  if ("probe_name" %in% colnames(adme)) {
    adme$probe_name <- simplify_probe_name(adme$probe_name)
  }
  adme
}

# Column order for Full SwissADME table / export: key identifiers first.
order_swissadme_display_columns <- function(df) {
  if (is.null(df) || ncol(df) == 0L) {
    return(df)
  }
  cn <- names(df)
  first <- character()
  if ("probe_name" %in% cn) {
    first <- c(first, "probe_name")
  }
  comp <- cn[cn %in% c("compound_names", "compound.names")][1]
  if (!is.na(comp)) {
    first <- c(first, comp)
  }
  smi <- cn[cn == "Canonical.SMILES"][1]
  if (is.na(smi)) {
    smi <- cn[grepl("canonical", cn, ignore.case = TRUE) &
      grepl("smile", cn, ignore.case = TRUE)][1]
  }
  if (!is.na(smi) && nzchar(smi)) {
    first <- c(first, smi)
  }
  fo <- cn[tolower(cn) == "formula"][1]
  if (!is.na(fo)) {
    first <- c(first, fo)
  }
  first <- first[first %in% cn]
  rest <- setdiff(cn, first)
  df[, c(first, rest), drop = FALSE]
}

# 3.2 filter_swissadme_data
# Subsets the SwissADME table to a single selected probe name.
filter_swissadme_data <- function(swissadme_table, selected_probe) {
  if (is.null(selected_probe)) return(NULL)
  swissadme_table %>%
    dplyr::filter(probe_name == selected_probe)
}

# 3.3 overlay_canonical_smiles_from_swissadme
# Protein binding table joins SMILES from Compound Keys (table-s2.xlsx), which may differ
# from SwissADME's canonical representation. Replace per probe_name so structures match
# the SwissADME tab (Canonical.SMILES).
canonical_smiles_by_probe <- function(adme) {
  if (is.null(adme) || nrow(adme) == 0 || !"probe_name" %in% colnames(adme)) {
    return(NULL)
  }
  cn <- colnames(adme)
  sm_col <- cn[grepl("canonical", cn, ignore.case = TRUE) &
    grepl("smiles", cn, ignore.case = TRUE)]
  if (length(sm_col) == 0) {
    return(NULL)
  }
  col <- sm_col[1]
  adme %>%
    dplyr::mutate(
      sm_tmp = trimws(as.character(.data[[col]]))
    ) %>%
    dplyr::filter(!is.na(sm_tmp), nzchar(sm_tmp)) %>%
    dplyr::distinct(probe_name, .keep_all = TRUE) %>%
    dplyr::select(probe_name, canonical_smiles = sm_tmp)
}

overlay_canonical_smiles_from_swissadme <- function(lookup, adme) {
  if (is.null(lookup) || nrow(lookup) == 0) {
    return(lookup)
  }
  map_df <- canonical_smiles_by_probe(adme)
  if (is.null(map_df) || nrow(map_df) == 0) {
    return(lookup)
  }
  lookup %>%
    dplyr::mutate(probe_name = simplify_probe_name(.data$probe_name)) %>%
    dplyr::left_join(map_df, by = "probe_name") %>%
    dplyr::mutate(
      SMILES = dplyr::coalesce(
        canonical_smiles,
        suppressWarnings(as.character(SMILES))
      )
    ) %>%
    dplyr::select(-dplyr::any_of("canonical_smiles"))
}

# ---- ggplot theme (used by app + modules) ----------------------
# In-app theme wrapper. Uses the branded theme_canprotarget() with
# app_config sizing overrides for the Shiny dashboard context.
cpt_theme <- function() {
  cfg <- if (exists("app_config", inherits = TRUE) && !is.null(app_config$plot_theme)) {
    app_config$plot_theme
  } else {
    list(base_size = 16)
  }

  theme_canprotarget(base_size = cfg$base_size)
}

# Per-subtype effect sizes: one file per dataset x subtype under
# data/precomputed_effectsizes/.
#
# Stored as .rds rather than .tsv. The same 141 tables are 288 MB instead of
# 850 MB and read about three times faster, which matters twice: the directory
# is most of the deployment bundle, and one of these files is read every time
# the subtype changes. .tsv is still read when present so a half-converted or
# externally generated directory keeps working.
cpt_effectsize_path <- function(data_dir, dataset_short, subtype) {
  stub <- file.path(data_dir, "precomputed_effectsizes",
                    paste0(dataset_short, "_", sanitize_subtype(subtype)))
  rds <- paste0(stub, ".rds")
  if (file.exists(rds)) rds else paste0(stub, ".tsv")
}

cpt_read_effectsizes <- function(path) {
  if (!file.exists(path)) return(NULL)
  if (grepl("\\.rds$", path)) readRDS(path)
  else readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)
}
