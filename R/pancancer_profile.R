# ============================================================
# Script:   pancancer_profile.R
# Purpose:  Pan-cancer dependency profile for a single gene.
#           Mean gene effect per cancer subtype, sorted by dependency
#           strength.
#
# Usage:    cpt_pancancer_profile() and api_pancancer_profile() serve the
#           MCP pancancer_profile tool and assess_target.
#           No Shiny dependencies.
# ============================================================
# Sections:
#   Pan-cancer data extraction
#     cpt_pancancer_profile()
#   MCP API wrapper
#     api_pancancer_profile()
# ============================================================

# ---- Pan-cancer data extraction --------------------------------

#' Extract mean gene effect per cancer subtype for one gene.
#'
#' Groups cell lines by OncotreeSubtype, computes mean gene effect,
#' and returns a summary table sorted by dependency strength.
#'
#' @param gene Character: gene symbol
#' @param matrix Matrix: gene effect matrix (rows = cell lines, cols = genes)
#' @param meta Data frame: cancer model metadata with ModelID + OncotreeSubtype
#' @param min_cell_lines Integer: minimum cell lines per subtype to include (default 3)
#' @return Data frame with columns: subtype, mean_effect, sd_effect, n_cell_lines, is_dependency
cpt_pancancer_profile <- function(gene, matrix, meta, min_cell_lines = 3L) {
  # Find gene column
  gene_names <- sub(" \\(\\d+\\)$", "", colnames(matrix))
  col_idx <- which(toupper(gene_names) == toupper(gene))
  if (!length(col_idx)) {
    stop("Gene '", gene, "' not found in the gene effect matrix.", call. = FALSE)
  }
  col_idx <- col_idx[1]
  actual_gene <- gene_names[col_idx]

  # Align metadata to matrix rows
  meta_aligned <- meta[match(rownames(matrix), meta$ModelID), ]

  # Gene effect values
  gene_values <- matrix[, col_idx]

  # Group by subtype
  subtypes <- meta_aligned$OncotreeSubtype
  valid <- !is.na(subtypes) & !is.na(gene_values)

  df <- data.frame(
    subtype = subtypes[valid],
    effect = gene_values[valid],
    stringsAsFactors = FALSE
  )

  # Aggregate per subtype
  agg <- do.call(rbind, lapply(split(df, df$subtype), function(sub_df) {
    data.frame(
      subtype = sub_df$subtype[1],
      mean_effect = mean(sub_df$effect, na.rm = TRUE),
      sd_effect = sd(sub_df$effect, na.rm = TRUE),
      n_cell_lines = nrow(sub_df),
      stringsAsFactors = FALSE
    )
  }))
  rownames(agg) <- NULL

  # Filter by minimum sample size

  agg <- agg[agg$n_cell_lines >= min_cell_lines, , drop = FALSE]

  # Add dependency flag and sort
  agg$is_dependency <- agg$mean_effect < -0.5
  agg <- agg[order(agg$mean_effect), ]
  rownames(agg) <- NULL

  attr(agg, "gene") <- actual_gene
  agg
}

# ---- MCP API wrapper -------------------------------------------

#' Pan-cancer profile for MCP/API use.
#' Returns structured data (not a plot) suitable for JSON serialization.
#'
#' @param gene Character: gene symbol
#' @param dataset Character: "CRISPR" or "RNAi"
#' @param data_env Environment/list with loaded data
#' @param top_n Integer: number of subtypes to return (default 30)
#' @return List with gene info and ranked subtype dependency table
api_pancancer_profile <- function(gene, dataset = "CRISPR", data_env, top_n = 30L) {
  if (exists("cpt_normalize_dataset", mode = "function")) {
    dataset <- cpt_normalize_dataset(dataset)
  } else {
    dataset <- toupper(as.character(dataset)[1])
    if (dataset == "RNAI") dataset <- "RNAi"
  }
  matrix <- if (dataset == "CRISPR") data_env$crispr_matrix else data_env$rnai_matrix
  meta <- data_env$cancer_model_data

  if (is.null(matrix)) stop("Dataset '", dataset, "' not loaded.", call. = FALSE)
  if (is.null(meta)) stop("Cancer model metadata not loaded.", call. = FALSE)

  profile <- cpt_pancancer_profile(gene, matrix, meta, min_cell_lines = 3L)
  actual_gene <- attr(profile, "gene")

  # Limit output (clamp; avoid negative head() quirks)
  if (exists("cpt_clamp_int", mode = "function")) {
    top_n <- cpt_clamp_int(top_n, min_val = 1L, max_val = 200L, default = 30L)
  } else {
    top_n <- as.integer(max(1L, min(as.integer(top_n)[1], 200L)))
    if (is.na(top_n)) top_n <- 30L
  }
  top_n <- min(top_n, nrow(profile))
  top_df <- utils::head(profile, top_n)

  # Convert to list of records
  records <- lapply(seq_len(nrow(top_df)), function(i) {
    row <- top_df[i, ]
    list(
      subtype = row$subtype,
      mean_effect = round(row$mean_effect, 4),
      sd = round(row$sd_effect, 4),
      n_cell_lines = as.integer(row$n_cell_lines),
      is_dependency = row$is_dependency
    )
  })

  n_dep <- sum(profile$is_dependency)

  list(
    gene = actual_gene,
    dataset = dataset,
    total_subtypes = nrow(profile),
    n_dependency_subtypes = as.integer(n_dep),
    showing_top_n = top_n,
    subtypes = records,
    interpretation = paste0(
      actual_gene, " is a dependency (effect < -0.5) in ", n_dep, " of ",
      nrow(profile), " cancer subtypes in the ", dataset, " dataset.",
      if (n_dep == 0) " This gene is NOT broadly essential." else
      if (n_dep > 10) " This gene appears broadly essential (possible common essential)." else
      paste0(" Top dependency: ", top_df$subtype[1], " (mean effect: ", round(top_df$mean_effect[1], 3), ").")
    )
  )
}
