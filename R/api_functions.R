# ============================================================
# Script:   api_functions.R
# Purpose:  Pure query functions for MCP/API access to CanProTarget.
#           These are stateless functions that operate on pre-loaded data.
#           No Shiny dependencies. Called by R/mcp_worker.R.
#
# Note:     Core dependency analysis builds on the original cell_cpt
#           Shiny app (John-Paul Ong). MCP/API wrappers and composite
#           scoring are platform extensions layered on top of that work.
# ============================================================
# Sections:
#   Helpers
#     cpt_normalize_dataset(), cpt_clamp_int(), parse_gene_names(), find_gene_column(), +2 more
#   Tool: list_subtypes
#     api_list_subtypes()
#   Tool: list_genes
#     api_list_genes()
#   Tool: query_dependency
#     api_query_dependency()
#   Tool: top_dependencies
#     api_top_dependencies()
#   Tool: gene_cysteines
#     api_gene_cysteines()
#   Tool: cysteine_detail
#     api_cysteine_detail()
#   Tool: compare_subtypes
#     api_compare_subtypes()
#   Tool: platform_info
#     api_platform_info()
#   Tool: canprotarget_score
#     api_canprotarget_score()
#   Tool: rank_targets
#     api_rank_targets(), api_rank_site_targets()
#   Tool: assess_target
#     api_assess_target()
# ============================================================

# ---- Helpers ---------------------------------------------------

CPT_PLATFORM_VERSION <- "1.0.0"
CPT_CITATION <- "Ong JP, Martins D, Bertoldo JB. CanProTarget. Zenodo. DOI pending release."
CPT_REPO <- "https://github.com/Bertoldo-Lab/CanProTarget"

#' Normalize dataset label to "CRISPR" or "RNAi".
#' Rejects typos rather than falling through silently to the RNAi matrix.
cpt_normalize_dataset <- function(dataset) {
  d <- toupper(trimws(as.character(dataset)[1]))
  if (is.na(d) || !nzchar(d)) {
    stop("dataset must be 'CRISPR' or 'RNAi'", call. = FALSE)
  }
  # Accept display labels like "CRISPR (23Q4)"
  if (startsWith(d, "CRISPR")) return("CRISPR")
  if (d %in% c("RNAI", "RNAi") || startsWith(d, "RNAI") || startsWith(d, "DEMETER")) {
    return("RNAi")
  }
  stop("dataset must be 'CRISPR' or 'RNAi' (got '", dataset, "')", call. = FALSE)
}

#' Clamp integer params used by top_n / pool tools.
cpt_clamp_int <- function(x, min_val = 1L, max_val = 100L, default = 20L) {
  n <- suppressWarnings(as.integer(x)[1])
  if (is.na(n)) n <- as.integer(default)
  as.integer(max(min_val, min(n, max_val)))
}

# Parse gene names from matrix column format "GENE (ENTREZ_ID)" -> "GENE"
parse_gene_names <- function(col_names) {

  sub(" \\(\\d+\\)$", "", col_names)
}

# Find column index by gene name (case-insensitive)
find_gene_column <- function(matrix, gene_name) {
  gene_names <- parse_gene_names(colnames(matrix))
  idx <- which(toupper(gene_names) == toupper(gene_name))
  if (length(idx) == 0L) return(NULL)
  idx[1]
}

#' Shared provenance block for agent/report consumers.
#' Keeps outputs auditable without turning us into a general workbench.
cpt_provenance <- function(extra = NULL) {
  out <- list(
    platform = "CanProTarget",
    version = CPT_PLATFORM_VERSION,
    citation = CPT_CITATION,
    repository = CPT_REPO,
    original_app = "John-Paul Ong, cell_cpt (original Shiny foundation)",
    data_snapshot = list(
      crispr = "DepMap 23Q4",
      rnai = "DEMETER2 v6",
      cysteine_atlas = "Li et al. Nat Chem Biol 2023"
    )
  )
  if (!is.null(extra) && is.list(extra)) {
    out <- c(out, extra)
  }
  out
}

#' Follow-up links for structure / portals (no embedded viewers).
#' @param gene Gene symbol
#' @param uniprot Optional UniProt accession
#' @return Named list of URLs
cpt_external_resources <- function(gene, uniprot = NULL) {
  gene <- as.character(gene)[1]
  links <- list(
    depmap = paste0("https://depmap.org/portal/gene/", utils::URLencode(gene, reserved = TRUE), "?tab=overview"),
    genecards = paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", utils::URLencode(gene, reserved = TRUE)),
    oncokb = paste0("https://www.oncokb.org/gene/", utils::URLencode(gene, reserved = TRUE)),
    pdb_search = paste0("https://www.rcsb.org/search?q=", utils::URLencode(gene, reserved = TRUE))
  )
  if (!is.null(uniprot) && length(uniprot) > 0 && !is.na(uniprot[1]) && nzchar(uniprot[1])) {
    u <- as.character(uniprot[1])
    links$uniprot <- paste0("https://www.uniprot.org/uniprot/", u)
    links$alphafold <- paste0("https://alphafold.ebi.ac.uk/entry/", u)
  } else {
    links$uniprot_search <- paste0(
      "https://www.uniprot.org/uniprotkb?query=",
      utils::URLencode(gene, reserved = TRUE)
    )
  }
  links
}

# ---- Tool: list_subtypes ---------------------------------------

#' List available cancer subtypes for a given dataset.
#' @param dataset Character: "CRISPR" or "RNAi"
#' @param data_env Environment/list containing loaded data
#' @return List with dataset, count, and subtypes vector
api_list_subtypes <- function(dataset, data_env) {
  dataset <- cpt_normalize_dataset(dataset)

  file_key <- if (dataset == "CRISPR") "cancer_subtypes_CRISPR.txt" else "cancer_subtypes_RNAi.txt"
  path <- file.path(data_env$data_dir, file_key)

  if (!file.exists(path)) {
    stop("Subtype list file not found: ", file_key, call. = FALSE)
  }

  subtypes <- readLines(path, warn = FALSE, encoding = "UTF-8")
  subtypes <- trimws(subtypes)
  subtypes <- subtypes[nzchar(subtypes) & !grepl("^#", subtypes)]
  subtypes <- sort(unique(subtypes))

  list(
    dataset = dataset,
    count = length(subtypes),
    subtypes = subtypes
  )
}

# ---- Tool: list_genes ------------------------------------------

#' List available genes in a dataset.
#' @param dataset Character: "CRISPR" or "RNAi"
#' @param data_env Environment/list containing loaded data
#' @return List with dataset, count, and genes vector
api_list_genes <- function(dataset, data_env) {
  dataset <- cpt_normalize_dataset(dataset)
  matrix <- if (dataset == "CRISPR") data_env$crispr_matrix else data_env$rnai_matrix

  if (is.null(matrix)) {
    stop("Dataset '", dataset, "' is not loaded.", call. = FALSE)
  }

  genes <- sort(unique(parse_gene_names(colnames(matrix))))

  list(
    dataset = dataset,
    count = length(genes),
    genes = genes
  )
}

# ---- Tool: query_dependency ------------------------------------

#' Get dependency score for a specific gene in a cancer subtype.
#' Uses mean difference + one-sided Welch/Student t-test (fast single-gene path).
#' Note: the Shiny Dependencies table may use limma/precomputed TSVs; p-values
#' are not guaranteed to match that UI path for the same gene/subtype.
#' @param gene Character: gene symbol (e.g. "KRAS")
#' @param subtype Character: OncotreeSubtype (e.g. "Non-Small Cell Lung Cancer")
#' @param dataset Character: "CRISPR" or "RNAi"
#' @param data_env Environment/list containing loaded data
#' @return List with gene, subtype, effect_size, p_value, stats
api_query_dependency <- function(gene, subtype, dataset = "CRISPR", data_env) {
  dataset <- cpt_normalize_dataset(dataset)
  matrix <- if (dataset == "CRISPR") data_env$crispr_matrix else data_env$rnai_matrix
  meta <- data_env$cancer_model_data

  if (is.null(matrix)) stop("Dataset '", dataset, "' not loaded.", call. = FALSE)
  if (is.null(meta)) stop("Cancer model metadata not loaded.", call. = FALSE)

  # Find the gene column
  col_idx <- find_gene_column(matrix, gene)
  if (is.null(col_idx)) {
    stop("Gene '", gene, "' not found in ", dataset, " dataset.", call. = FALSE)
  }
  actual_gene <- parse_gene_names(colnames(matrix)[col_idx])

  # Validate subtype exists
  subtypes_in_meta <- unique(meta$OncotreeSubtype[!is.na(meta$OncotreeSubtype)])
  if (!subtype %in% subtypes_in_meta) {
    stop("Subtype '", subtype, "' not found in metadata.", call. = FALSE)
  }

  # Validate sample sizes (non-NA gene-effect cells, not just metadata counts)
  model_ids_in_matrix <- rownames(matrix)
  meta_aligned <- meta[match(model_ids_in_matrix, meta$ModelID), ]
  Y <- ifelse(meta_aligned$OncotreeSubtype %in% subtype, 1, 0)
  Y[is.na(Y)] <- 0

  gene_values <- matrix[, col_idx]
  cancer_values <- gene_values[Y == 1]
  other_values <- gene_values[Y == 0]
  cancer_values <- cancer_values[!is.na(cancer_values)]
  other_values <- other_values[!is.na(other_values)]

  n_in_subtype <- length(cancer_values)
  n_other <- length(other_values)

  if (n_in_subtype < 3L) {
    stop("Only ", n_in_subtype, " non-NA cell line(s) match subtype '", subtype,
         "'. Need at least 3.", call. = FALSE)
  }
  if (n_other < 3L) {
    stop("Only ", n_other, " non-NA cell line(s) outside subtype '", subtype,
         "'. Need at least 3 for comparison.", call. = FALSE)
  }

  cancer_mean <- mean(cancer_values)
  other_mean <- mean(other_values)
  effect_size <- cancer_mean - other_mean

  # t-test: cancer vs others (one-sided, more essential = more negative)
  t_result <- tryCatch(
    t.test(cancer_values, other_values, alternative = "less"),
    error = function(e) NULL
  )
  p_value <- if (!is.null(t_result)) t_result$p.value else NA_real_

  # Interpret
  is_dependency <- !is.na(cancer_mean) && cancer_mean < -0.5
  is_selective <- !is.na(effect_size) && effect_size < -0.1 &&
    !is.na(p_value) && p_value < 0.05

  list(
    gene = actual_gene,
    subtype = subtype,
    dataset = dataset,
    n_cell_lines_in_subtype = as.integer(n_in_subtype),
    n_cell_lines_other = as.integer(n_other),
    cancer_mean_effect = round(cancer_mean, 4),
    other_mean_effect = round(other_mean, 4),
    effect_size = round(effect_size, 4),
    p_value = signif(p_value, 4),
    is_dependency = unname(is_dependency),
    is_selective_dependency = unname(is_selective),
    interpretation = if (is_selective) {
      paste0(actual_gene, " is a SELECTIVE dependency in ", subtype,
             " (effect size: ", round(effect_size, 3), ", p=", signif(p_value, 3), ")")
    } else if (is_dependency) {
      paste0(actual_gene, " is a dependency in ", subtype,
             " but NOT selective vs other cancers (mean effect: ", round(cancer_mean, 3), ")")
    } else {
      paste0(actual_gene, " is NOT a dependency in ", subtype,
             " (mean effect: ", round(cancer_mean, 3), ")")
    },
    methods_note = paste0(
      "Mean gene effect in ", subtype, " cell lines vs all others; one-sided t-test ",
      "(cancer more essential). Not the Shiny limma/precomputed path. ",
      "Thresholds: dependency mean < -0.5; selective if effect size < -0.1 and p < 0.05."
    ),
    external_resources = cpt_external_resources(actual_gene),
    provenance = cpt_provenance(list(dataset = dataset, analysis = "query_dependency"))
  )
}

# ---- Tool: top_dependencies ------------------------------------

#' Top genes by mean effect-size (cancer - others).
#' Ranked by selectivity of the mean difference only — not p-values and not
#' CPT Score. Use rank_targets for composite prioritization; use query_dependency
#' for per-gene significance.
#' @param subtype Character: OncotreeSubtype
#' @param dataset Character: "CRISPR" or "RNAi"
#' @param n Integer: number of top genes to return (default 20, max 100)
#' @param data_env Environment/list containing loaded data
#' @return List with subtype, dataset, and ranked gene table
api_top_dependencies <- function(subtype, dataset = "CRISPR", n = 20L, data_env) {
  dataset <- cpt_normalize_dataset(dataset)
  matrix <- if (dataset == "CRISPR") data_env$crispr_matrix else data_env$rnai_matrix
  meta <- data_env$cancer_model_data

  if (is.null(matrix)) stop("Dataset '", dataset, "' not loaded.", call. = FALSE)
  if (is.null(meta)) stop("Cancer model metadata not loaded.", call. = FALSE)

  n <- cpt_clamp_int(n, min_val = 1L, max_val = 100L, default = 20L)

  # Align metadata
  model_ids_in_matrix <- rownames(matrix)
  meta_aligned <- meta[match(model_ids_in_matrix, meta$ModelID), ]
  cancer_idx <- which(meta_aligned$OncotreeSubtype %in% subtype)

  if (length(cancer_idx) < 3L) {
    stop("Only ", length(cancer_idx), " cell line(s) match subtype '", subtype,
         "'. Need at least 3.", call. = FALSE)
  }

  other_idx <- which(!meta_aligned$OncotreeSubtype %in% subtype)
  if (length(other_idx) < 3L) {
    stop("Only ", length(other_idx), " cell line(s) outside subtype '", subtype,
         "'. Need at least 3 for comparison.", call. = FALSE)
  }

  # Mean effect in cancer subtype
  cancer_means <- colMeans(matrix[cancer_idx, , drop = FALSE], na.rm = TRUE)
  other_means <- colMeans(matrix[other_idx, , drop = FALSE], na.rm = TRUE)
  effect_sizes <- cancer_means - other_means

  # Build results
  gene_names <- parse_gene_names(names(cancer_means))
  results <- data.frame(
    gene = gene_names,
    cancer_mean = round(as.numeric(cancer_means), 4),
    other_mean = round(as.numeric(other_means), 4),
    effect_size = round(as.numeric(effect_sizes), 4),
    is_dependency = as.numeric(cancer_means) < -0.5,
    stringsAsFactors = FALSE
  )

  # Sort by effect size (most negative = most selectively essential)
  results <- results[order(results$effect_size), ]
  top_results <- utils::head(results, n)
  rownames(top_results) <- NULL

  # Convert to list of records for JSON
  records <- lapply(seq_len(nrow(top_results)), function(i) {
    as.list(top_results[i, ])
  })

  list(
    subtype = subtype,
    dataset = dataset,
    n_cell_lines = length(cancer_idx),
    n_genes_total = ncol(matrix),
    top_n = n,
    ranking = "effect_size_mean_difference",
    methods_note = paste0(
      "Ranked by mean effect size (subtype mean - other mean) only. ",
      "No p-values and not CPT Score. is_dependency flags cancer_mean < -0.5. ",
      "For composite prioritization use rank_targets; for significance use query_dependency."
    ),
    genes = records
  )
}

# ---- Tool: gene_cysteines --------------------------------------

#' Get all cysteine sites for a gene from the Cys_editing atlas.
#' @param gene Character: gene symbol (e.g. "EGFR")
#' @param data_env Environment/list containing loaded data
#' @return List with gene info and cysteine sites
api_gene_cysteines <- function(gene, data_env) {
  bundle <- cpt_gene_site_records(gene, data_env, min_cr = 4, smcl_cap = 15L)
  if (is.null(bundle) || !length(bundle$sites)) {
    stop("Gene '", gene, "' not found in cysteine editing atlas or chemoproteomic table.",
         call. = FALSE)
  }

  actual_gene <- bundle$gene
  protein_name <- bundle$protein_name
  uniprot_acc <- bundle$uniprot_accession

  sites <- lapply(bundle$sites, function(s) {
    list(
      site_id = s$site_id,
      cysteine_position = s$cysteine_position,
      in_atlas = s$in_atlas,
      engaged = s$engaged,
      evidence_tier = s$evidence_tier,
      evidence_tier_label = s$evidence_tier_label,
      functional = s$functional,
      ligandable = s$ligandable,
      functional_ligandable = s$functional_ligandable,
      editor_support = s$editor_support,
      study_context = s$study_context,
      ligandability_score = s$ligandability_score,
      conservation_score = s$conservation_score,
      clinvar_pathogenic = s$clinvar_pathogenic,
      n_smcls = s$n_smcls,
      n_smcls_prioritised = s$n_smcls_prioritised,
      smcls = s$smcls
    )
  })

  n_functional <- sum(vapply(sites, function(s) isTRUE(s$functional), logical(1)))
  n_ligandable <- sum(vapply(sites, function(s) isTRUE(s$ligandable), logical(1)))
  n_func_lig <- sum(vapply(sites, function(s) isTRUE(s$functional_ligandable), logical(1)))
  n_engaged <- sum(vapply(sites, function(s) isTRUE(s$engaged), logical(1)))
  n_tier <- function(t) sum(vapply(sites, function(s) identical(as.integer(s$evidence_tier), as.integer(t)), logical(1)))

  list(
    gene = actual_gene,
    protein_name = protein_name,
    uniprot_accession = uniprot_acc,
    total_sites = length(sites),
    n_atlas = as.integer(bundle$n_atlas),
    n_engaged = as.integer(n_engaged),
    n_functional = as.integer(n_functional),
    n_ligandable = as.integer(n_ligandable),
    n_functional_ligandable = as.integer(n_func_lig),
    n_tier_1 = as.integer(n_tier(1L)),
    n_tier_2 = as.integer(n_tier(2L)),
    n_tier_3 = as.integer(n_tier(3L)),
    n_tier_4 = as.integer(n_tier(4L)),
    summary = paste0(
      actual_gene, " has ", length(sites), " cysteine site(s) (",
      bundle$n_atlas, " in Cys_editing atlas, ", n_engaged,
      " engaged at CR >= 4). Atlas: ", n_functional, " functional, ",
      n_ligandable, " ligandable, ", n_func_lig, " functional+ligandable. ",
      "Tiers 1-4: ", n_tier(1L), "/", n_tier(2L), "/", n_tier(3L), "/", n_tier(4L),
      ". Sites ordered engaged-by-tier, then functional+ligandable first."
    ),
    sites = sites,
    external_resources = cpt_external_resources(actual_gene, uniprot_acc),
    provenance = cpt_provenance()
  )
}

# ---- Tool: cysteine_detail -------------------------------------

#' Get full annotation for a specific cysteine site.
#' @param site_id Character: e.g. "EGFR_797" or gene + position
#' @param data_env Environment/list containing loaded data
#' @return List with all annotations for the site
api_cysteine_detail <- function(site_id, data_env) {
  atlas <- data_env$cys_atlas

  if (is.null(atlas)) {
    stop("Cysteine editing atlas not loaded.", call. = FALSE)
  }

  # Try exact match first
  hit <- atlas[atlas$site_id == site_id, , drop = FALSE]

  # If not found, try parsing as GENE_POS
  if (nrow(hit) == 0L) {
    parts <- strsplit(toupper(site_id), "_")[[1]]
    if (length(parts) >= 2) {
      gene <- paste(parts[-length(parts)], collapse = "_")
      pos <- suppressWarnings(as.integer(parts[length(parts)]))
      if (!is.na(pos)) {
        hit <- atlas[toupper(atlas$gene_symbol) == gene &
                     atlas$cysteine_position == pos, , drop = FALSE]
      }
    }
  }

  if (nrow(hit) == 0L) {
    parsed <- cpt_parse_site_id(site_id)
    bundle <- if (!is.na(parsed$pos) && nzchar(parsed$gene)) {
      tryCatch(
        cpt_gene_site_records(parsed$gene, data_env, min_cr = 4, smcl_cap = 40L),
        error = function(e) NULL
      )
    } else {
      NULL
    }
    rec <- NULL
    if (!is.null(bundle) && length(bundle$sites)) {
      rec <- Filter(function(s) {
        identical(as.integer(s$cysteine_position), as.integer(parsed$pos)) ||
          identical(s$site_id, site_id)
      }, bundle$sites)
      if (length(rec)) rec <- rec[[1]]
    }
    if (is.null(rec) || !length(rec)) {
      stop("Site '", site_id, "' not found in cysteine editing atlas or chemoproteomic table.",
           call. = FALSE)
    }
    return(list(
      site_id = rec$site_id,
      gene_symbol = rec$gene_symbol,
      cysteine_position = rec$cysteine_position,
      uniprot_accession = rec$uniprot_accession,
      in_atlas = FALSE,
      engaged = rec$engaged,
      evidence_tier = rec$evidence_tier,
      evidence_tier_label = rec$evidence_tier_label,
      functional = FALSE,
      ligandable = FALSE,
      functional_ligandable = FALSE,
      n_smcls = rec$n_smcls,
      n_smcls_prioritised = rec$n_smcls_prioritised,
      smcls = rec$smcls,
      note = paste(
        "Site is not in the Cys_editing atlas (Tier 4, untested). Absence from the atlas",
        "is missing coverage, not evidence of non-function."
      ),
      external_resources = cpt_external_resources(rec$gene_symbol, rec$uniprot_accession),
      provenance = cpt_provenance()
    ))
  }

  row <- hit[1, ]
  bundle <- tryCatch(
    cpt_gene_site_records(row$gene_symbol, data_env, min_cr = 4, smcl_cap = 40L),
    error = function(e) NULL
  )
  rec <- NULL
  if (!is.null(bundle) && length(bundle$sites)) {
    matched <- Filter(function(s) {
      identical(as.integer(s$cysteine_position), as.integer(row$cysteine_position))
    }, bundle$sites)
    if (length(matched)) rec <- matched[[1]]
  }
  tier <- if (!is.null(rec)) rec$evidence_tier else {
    cpt_evidence_tier(TRUE, isTRUE(row$functional), isTRUE(row$ligandable))
  }

  list(
    site_id = row$site_id,
    gene_symbol = row$gene_symbol,
    protein_name = row$protein_name,
    cysteine_position = as.integer(row$cysteine_position),
    uniprot_accession = row$uniprot_accession,
    residue_mapping_status = row$residue_mapping_status,
    study_context = row$study_context,
    editor_support = row$editor_support,
    # Functional assessment
    functional = as.logical(row$functional),
    abe_functional = as.logical(row$abe_functional),
    cbe_functional = as.logical(row$cbe_functional),
    abe_mean_lfc = round(row$abe_mean_lfc, 3),
    cbe_mean_lfc = round(row$cbe_mean_lfc, 3),
    abe_guide_count = as.integer(row$abe_guide_count),
    cbe_guide_count = as.integer(row$cbe_guide_count),
    abe_neg_log10_p = round(row$abe_neg_log10_p, 3),
    cbe_neg_log10_p = round(row$cbe_neg_log10_p, 3),
    abe_neg_log10_fdr = round(row$abe_neg_log10_fdr, 3),
    cbe_neg_log10_fdr = round(row$cbe_neg_log10_fdr, 3),
    # Ligandability
    ligandability_score = if (is.na(row$ligandability_score)) NULL else round(row$ligandability_score, 1),
    ligandable = as.logical(row$ligandable),
    functional_ligandable = as.logical(row$functional_ligandable),
    # Structure
    proteomic_accessibility = row$proteomic_accessibility,
    relative_solvent_accessibility = if (is.na(row$relative_solvent_accessibility)) NULL else round(row$relative_solvent_accessibility, 3),
    alphafold_plddt = if (is.na(row$alphafold_plddt)) NULL else round(row$alphafold_plddt, 1),
    pdb_coverage = row$pdb_coverage,
    # Conservation
    conservation_score = if (is.na(row$conservation_score)) NULL else round(row$conservation_score, 2),
    ortholog_cys_count = if (is.na(row$ortholog_cys_count)) NULL else as.integer(row$ortholog_cys_count),
    ortholog_total = if (is.na(row$ortholog_total)) NULL else as.integer(row$ortholog_total),
    # Clinical
    clinvar_pathogenic = as.logical(row$clinvar_pathogenic),
    clinvar_annotation = if (is.na(row$clinvar_annotation)) NULL else row$clinvar_annotation,
    clinvar_phenotype = if (is.na(row$clinvar_phenotype)) NULL else row$clinvar_phenotype,
    in_atlas = TRUE,
    engaged = if (!is.null(rec)) isTRUE(rec$engaged) else FALSE,
    evidence_tier = as.integer(tier),
    evidence_tier_label = cpt_evidence_tier_label(tier),
    n_smcls = if (!is.null(rec)) rec$n_smcls else 0L,
    n_smcls_prioritised = if (!is.null(rec)) rec$n_smcls_prioritised else 0L,
    smcls = if (!is.null(rec)) rec$smcls else list(),
    external_resources = cpt_external_resources(row$gene_symbol, row$uniprot_accession),
    provenance = cpt_provenance()
  )
}

# ---- Tool: compare_subtypes -----------------------------------

#' Compare dependency of a gene across two subtypes.
#' @param gene Character: gene symbol
#' @param subtype1 Character: first OncotreeSubtype
#' @param subtype2 Character: second OncotreeSubtype
#' @param dataset Character: "CRISPR" or "RNAi"
#' @param data_env Environment/list containing loaded data
#' @return List with comparative statistics
api_compare_subtypes <- function(gene, subtype1, subtype2, dataset = "CRISPR", data_env) {
  dataset <- cpt_normalize_dataset(dataset)
  matrix <- if (dataset == "CRISPR") data_env$crispr_matrix else data_env$rnai_matrix
  meta <- data_env$cancer_model_data

  if (is.null(matrix)) stop("Dataset '", dataset, "' not loaded.", call. = FALSE)
  if (is.null(meta)) stop("Cancer model metadata not loaded.", call. = FALSE)

  col_idx <- find_gene_column(matrix, gene)
  if (is.null(col_idx)) {
    stop("Gene '", gene, "' not found in ", dataset, " dataset.", call. = FALSE)
  }
  actual_gene <- parse_gene_names(colnames(matrix)[col_idx])

  # Align metadata
  meta_aligned <- meta[match(rownames(matrix), meta$ModelID), ]
  idx1 <- which(meta_aligned$OncotreeSubtype %in% subtype1)
  idx2 <- which(meta_aligned$OncotreeSubtype %in% subtype2)

  if (length(idx1) < 3L) stop("Subtype '", subtype1, "' has <3 cell lines.", call. = FALSE)
  if (length(idx2) < 3L) stop("Subtype '", subtype2, "' has <3 cell lines.", call. = FALSE)

  vals1 <- matrix[idx1, col_idx]
  vals2 <- matrix[idx2, col_idx]
  vals1 <- vals1[!is.na(vals1)]
  vals2 <- vals2[!is.na(vals2)]

  mean1 <- mean(vals1)
  mean2 <- mean(vals2)

  t_result <- tryCatch(
    t.test(vals1, vals2),
    error = function(e) NULL
  )
  p_value <- if (!is.null(t_result)) t_result$p.value else NA_real_

  list(
    gene = actual_gene,
    dataset = dataset,
    subtype1 = list(
      name = subtype1,
      n_cell_lines = length(idx1),
      mean_effect = round(mean1, 4),
      sd = round(sd(vals1), 4)
    ),
    subtype2 = list(
      name = subtype2,
      n_cell_lines = length(idx2),
      mean_effect = round(mean2, 4),
      sd = round(sd(vals2), 4)
    ),
    difference = round(mean1 - mean2, 4),
    p_value = signif(p_value, 4),
    interpretation = paste0(
      actual_gene, " effect in ", subtype1, ": ", round(mean1, 3),
      " vs ", subtype2, ": ", round(mean2, 3),
      " (diff=", round(mean1 - mean2, 3), ", p=", signif(p_value, 3), ")"
    )
  )
}

# ---- Tool: platform_info --------------------------------------

#' Return platform metadata, version, and citation info.
#' @return List with platform information
api_platform_info <- function() {
  list(
    name = "CanProTarget",
    version = CPT_PLATFORM_VERSION,
    description = paste(
      "Cancer Protein Target prioritization platform integrating DepMap dependency data",
      "with chemoproteomic ligandability and functional cysteine annotations.",
      "Specialist tool for covalent oncology target questions; not a general science workbench."
    ),
    data_sources = list(
      list(name = "DepMap CRISPR", version = "23Q4", genes = 18443L, cell_lines = 1100L),
      list(name = "DepMap RNAi (DEMETER2)", version = "v6", genes = 17309L, cell_lines = 712L),
      list(name = "Cysteine Editing Atlas", reference = "Li et al. Nat Chem Biol 2023", sites = 13872L),
      list(
        name = "Chemoproteomic competition ratios",
        reference = "Six CysDB-indexed ligandability studies (not a live CysDB query)",
        records = "CR >= 4 probe-cysteine pairs"
      )
    ),
    citation = CPT_CITATION,
    repository = CPT_REPO,
    credits = list(
      authors = "Ong JP, Martins D, Bertoldo JB",
      lab = "Bertoldo Lab, Children's Cancer Institute / UNSW Sydney",
      original_app = "John-Paul Ong (cell_cpt), the Shiny foundation this platform builds on"
    ),
    mcp_version = "1.2.0",
    available_tools = c(
      "list_subtypes", "list_genes", "query_dependency",
      "top_dependencies", "gene_cysteines", "cysteine_detail",
      "compare_subtypes", "platform_info", "canprotarget_score",
      "pancancer_profile", "generate_report",
      "assess_target", "rank_targets", "rank_site_targets"
    ),
    provenance = cpt_provenance()
  )
}


# ---- Tool: canprotarget_score ----------------------------------

#' Compute CanProTarget Score for a gene in a specific cancer subtype.
#' @param gene Character: gene symbol
#' @param subtype Character: OncotreeSubtype
#' @param dataset Character: "CRISPR" or "RNAi"
#' @param data_env Environment/list containing loaded data
#' @return List with composite score and dimension breakdown
api_canprotarget_score <- function(gene, subtype, dataset = "CRISPR", data_env) {
  dataset <- cpt_normalize_dataset(dataset)
  result <- cpt_score_single(gene, subtype, dataset, data_env)
  # Enrich with shared links + provenance (score function stays pure-ish)
  uniprot <- NULL
  if (!is.null(data_env$cys_atlas)) {
    hit <- data_env$cys_atlas[
      toupper(data_env$cys_atlas$gene_symbol) == toupper(result$gene),
      "uniprot_accession",
      drop = TRUE
    ]
    if (length(hit) > 0 && !is.na(hit[1]) && nzchar(hit[1])) uniprot <- hit[1]
  }
  result$external_resources <- cpt_external_resources(result$gene, uniprot)
  result$provenance <- cpt_provenance(list(
    dataset = result$dataset %||% dataset,
    analysis = "canprotarget_score",
    ranking = "dependency/selectivity percentiles vs all genes in subtype"
  ))
  result$methods_note <- paste(
    "CPT Score: weighted composite of available dimensions (dependency strength,",
    "cancer selectivity, cysteine ligandability, conservation, clinical evidence,",
    "ADME when available). Dependency and selectivity use percentile ranks within",
    "the full subtype gene universe. Missing dimensions are listed explicitly."
  )
  result
}


# ---- Tool: rank_targets ----------------------------------------

#' Rank selective dependencies in a subtype by CPT Score.
#' @param subtype Character: OncotreeSubtype
#' @param dataset Character: "CRISPR" or "RNAi"
#' @param n Integer: how many genes to return
#' @param pool Integer: size of selective-dependency pool to re-score
#' @param require_ligandable Logical: keep only genes with ligandability data
#' @param data_env Environment/list containing loaded data
#' @return Ranked list with CPT breakdowns
api_rank_targets <- function(subtype,
                             dataset = "CRISPR",
                             n = 20L,
                             pool = 100L,
                             require_ligandable = FALSE,
                             data_env) {
  out <- cpt_rank_targets(
    subtype = subtype,
    dataset = dataset,
    n = n,
    pool = pool,
    require_ligandable = require_ligandable,
    data_env = data_env
  )
  if (is.null(out$genes) || !length(out$genes) || is.null(data_env$smcl_index)) {
    return(out)
  }
  out$genes <- lapply(out$genes, function(g) {
    bundle <- tryCatch(
      cpt_gene_site_records(g$gene, data_env, min_cr = 4, smcl_cap = 5L),
      error = function(e) NULL
    )
    engaged <- if (!is.null(bundle)) {
      Filter(function(s) isTRUE(s$engaged), bundle$sites)
    } else {
      list()
    }
    if (length(engaged)) {
      o <- order(vapply(engaged, function(s) s$evidence_tier, integer(1)),
                 -vapply(engaged, function(s) s$n_smcls, integer(1)))
      best <- engaged[[o[1]]]
      top <- if (length(best$smcls)) best$smcls[[1]] else NULL
      g$best_engaged_site <- list(
        site_id = best$site_id,
        cysteine_position = best$cysteine_position,
        evidence_tier = best$evidence_tier,
        evidence_tier_label = best$evidence_tier_label,
        n_smcls = best$n_smcls,
        top_probe = if (is.null(top)) NULL else top$probe_name,
        top_CR = if (is.null(top)) NULL else top$CR
      )
    } else {
      g$best_engaged_site <- NULL
    }
    g
  })
  out
}

#' Rank engaged cysteines in a subtype by evidence tier, then Site CPT.
api_rank_site_targets <- function(subtype,
                                  dataset = "CRISPR",
                                  n = 20L,
                                  pool = 100L,
                                  min_cr = 4,
                                  max_targets = 20,
                                  effect_size_max = NULL,
                                  p_max = NULL,
                                  exclude_common_essentials = FALSE,
                                  data_env) {
  cpt_rank_site_targets(
    subtype = subtype,
    dataset = dataset,
    n = n,
    pool = pool,
    min_cr = min_cr,
    max_targets = max_targets,
    effect_size_max = effect_size_max,
    p_max = p_max,
    exclude_common_essentials = isTRUE(exclude_common_essentials),
    data_env = data_env
  )
}


# ---- Tool: assess_target ---------------------------------------

#' One-shot target assessment for agent workflows.
#' Bundles dependency stats, CPT score, cysteine summary, and suggested next steps.
#' @param gene Character: gene symbol
#' @param subtype Character: OncotreeSubtype
#' @param dataset Character: "CRISPR" or "RNAi"
#' @param data_env Environment/list containing loaded data
#' @return Nested list suitable for JSON / agent chaining
api_assess_target <- function(gene, subtype, dataset = "CRISPR", data_env) {
  dataset <- cpt_normalize_dataset(dataset)

  dependency <- api_query_dependency(gene, subtype, dataset, data_env)
  score <- api_canprotarget_score(gene, subtype, dataset, data_env)

  # Cysteine summary is optional (gene may not be in atlas)
  cysteines <- tryCatch(
    api_gene_cysteines(gene, data_env),
    error = function(e) NULL
  )

  # Pan-cancer: keep it light (top 10) for agent context
  pancancer <- tryCatch(
    api_pancancer_profile(gene, dataset, data_env, top_n = 10L),
    error = function(e) NULL
  )

  next_actions <- list()
  actual_gene <- dependency$gene

  site_bundle <- tryCatch(
    cpt_gene_site_records(actual_gene, data_env, min_cr = 4, smcl_cap = 8L),
    error = function(e) NULL
  )
  if (!is.null(site_bundle) && length(site_bundle$sites)) {
    site_bundle <- cpt_attach_site_cpt(site_bundle, score)
  }
  engaged_sites <- if (!is.null(site_bundle)) {
    Filter(function(s) isTRUE(s$engaged), site_bundle$sites)
  } else {
    list()
  }

  if (length(engaged_sites) > 0) {
    best_eng <- engaged_sites[[1]]
    next_actions <- c(next_actions, list(list(
      tool = "cysteine_detail",
      reason = paste0(
        "Inspect engaged residue ", best_eng$site_id,
        " (Tier ", best_eng$evidence_tier, ") and its SMCL records"
      ),
      params = list(site_id = best_eng$site_id)
    )))
  } else if (!is.null(cysteines) && length(cysteines$sites) > 0) {
    fl_sites <- Filter(function(s) isTRUE(s$functional_ligandable), cysteines$sites)
    if (length(fl_sites) > 0) {
      next_actions <- c(next_actions, list(list(
        tool = "cysteine_detail",
        reason = "Inspect a functional+ligandable cysteine for covalent handle detail",
        params = list(site_id = fl_sites[[1]]$site_id)
      )))
    } else {
      next_actions <- c(next_actions, list(list(
        tool = "gene_cysteines",
        reason = "Review annotated cysteines for this gene",
        params = list(gene = actual_gene)
      )))
    }
  }

  next_actions <- c(next_actions, list(
    list(
      tool = "pancancer_profile",
      reason = "Check whether dependency is subtype-selective or pan-essential",
      params = list(gene = actual_gene, dataset = dataset)
    ),
    list(
      tool = "generate_report",
      reason = "Produce a shareable HTML gene dependency report",
      params = list(
        report_type = "gene_dependency",
        gene = actual_gene,
        subtype = subtype,
        dataset = dataset
      )
    )
  ))

  # Plain-language summary for agents to relay
  cys_bit <- if (is.null(cysteines)) {
    "No cysteine atlas entry for this gene."
  } else {
    paste0(
      cysteines$total_sites, " cysteine site(s) annotated (",
      cysteines$n_functional_ligandable, " functional+ligandable)."
    )
  }

  summary <- paste(
    dependency$interpretation,
    score$interpretation,
    cys_bit,
    sep = " "
  )

  list(
    gene = actual_gene,
    subtype = subtype,
    dataset = dataset,
    summary = summary,
    dependency = dependency,
    cpt_score = score,
    cysteines = if (is.null(cysteines)) {
      list(available = FALSE, note = "Gene not found in cysteine editing atlas.")
    } else {
      list(
        available = TRUE,
        total_sites = cysteines$total_sites,
        n_functional = cysteines$n_functional,
        n_ligandable = cysteines$n_ligandable,
        n_functional_ligandable = cysteines$n_functional_ligandable,
        n_engaged = cysteines$n_engaged,
        n_tier_1 = cysteines$n_tier_1,
        n_tier_2 = cysteines$n_tier_2,
        n_tier_3 = cysteines$n_tier_3,
        n_tier_4 = cysteines$n_tier_4,
        summary = cysteines$summary,
        # Keep payload small: top sites only
        top_sites = utils::head(cysteines$sites, 10)
      )
    },
    engaged_sites = lapply(utils::head(engaged_sites, 10), function(s) {
      list(
        site_id = s$site_id,
        cysteine_position = s$cysteine_position,
        evidence_tier = s$evidence_tier,
        evidence_tier_label = s$evidence_tier_label,
        site_cpt_score = s$site_cpt_score,
        n_smcls = s$n_smcls,
        n_smcls_prioritised = s$n_smcls_prioritised,
        top_smcl = if (length(s$smcls)) s$smcls[[1]] else NULL,
        smcls = s$smcls
      )
    }),
    pancancer_top = if (is.null(pancancer)) NULL else {
      n_dep <- pancancer$n_dependency_subtypes
      sel_class <- if (n_dep == 0) {
        "not_a_dependency"
      } else if (n_dep <= 3) {
        "selective"
      } else if (n_dep <= 10) {
        "moderate_breadth"
      } else {
        "broad_or_common_essential"
      }
      list(
        total_subtypes = pancancer$total_subtypes,
        n_dependency_subtypes = n_dep,
        selectivity_class = sel_class,
        interpretation = pancancer$interpretation,
        top_subtypes = pancancer$subtypes
      )
    },
    next_actions = next_actions,
    external_resources = if (!is.null(score$external_resources)) {
      score$external_resources
    } else {
      cpt_external_resources(actual_gene)
    },
    citation = CPT_CITATION,
    provenance = cpt_provenance(list(analysis = "assess_target", dataset = dataset)),
    caveats = c(
      score$caveats,
      "assess_target is a convenience wrapper; drill into individual tools for full detail."
    )
  )
}
