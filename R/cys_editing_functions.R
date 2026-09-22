# ============================================================
# Script:   cys_editing_functions.R
# Purpose:  Pure helpers for importing the Cravatt Lab
#           Cys_editing functional-cysteine atlas.
# ============================================================
# Sections:
#   cpt_cys_required_columns()
#   cpt_cys_within_gene_fdr()
#   cpt_cys_true()
#   cpt_cys_build_atlas()
#   cpt_cys_annotate_probe_table()
#   cpt_cys_link_functional_sites()
# ============================================================

# This file is also sourced standalone (tests/testthat/test-cys-editing-atlas.R,
# docs/scripts/import_cys_editing_atlas.R) without R/functions.R, so keep a
# fallback for the shared gene-key normaliser. Canonical definition lives in
# R/functions.R; keep the two in sync.
if (!exists("cpt_gene_match_key", mode = "function")) {
  cpt_gene_match_key <- function(x) {
    toupper(trimws(sub("\\s*\\(\\d+\\)$", "", as.character(x))))
  }
}

cpt_cys_required_columns <- function(data, columns, object_name) {
  missing_columns <- setdiff(columns, colnames(data))
  if (length(missing_columns)) {
    stop(
      object_name, " is missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# The publication code reports empirical significance as -log10(p) and applies
# Benjamini-Hochberg correction separately to the cysteines within each gene.
cpt_cys_within_gene_fdr <- function(neg_log10_p, gene_symbol) {
  result <- rep(0, length(neg_log10_p))
  gene_indices <- split(seq_along(gene_symbol), gene_symbol, drop = TRUE)

  for (indices in gene_indices) {
    raw_p <- 10^(-neg_log10_p[indices])
    adjusted_p <- stats::p.adjust(raw_p, method = "fdr")
    result[indices] <- round(-log10(adjusted_p), 2)
  }

  result[is.na(result)] <- 0
  result
}

cpt_cys_true <- function(x) {
  !is.na(x) & x
}

# Reproduce filter_dropout_data(..., sel_filter = TRUE) from Cys_editing,
# then reshape the result into stable, user-facing CanProTarget fields.
cpt_cys_build_atlas <- function(cys_dropout,
                                dep_score,
                                af_rsa,
                                conserve_score = NULL,
                                ortho_count = NULL,
                                clinvar = NULL,
                                source_commit = NA_character_,
                                apply_context_filter = TRUE,
                                apply_manual_annotations = TRUE) {
  cpt_cys_required_columns(
    cys_dropout,
    c(
      "unique_ID", "HUGO", "target_Cys", "full_name", "cell",
      "mean_LFC_A", "mean_LFC_C", "A_guide_num", "C_guide_num",
      "p_A", "p_C", "KB_engage", "IADTB_access", "RSA",
      "denature_LFC", "PC_mean_LFC_A", "PC_mean_LFC_C",
      "KMS_mean_LFC_A", "KMS_mean_LFC_C"
    ),
    "Cys_dropout"
  )
  cpt_cys_required_columns(
    af_rsa,
    c(
      "symbol", "uniprot_accession", "residue_number", "aa",
      "aa_same_as_portal", "pLDDT", "covered_by_PDB", "covered_by_PDB_30"
    ),
    "AF_RSA"
  )

  if (!is.matrix(dep_score) && !is.data.frame(dep_score)) {
    stop("ceres_used must be a matrix or data frame.", call. = FALSE)
  }
  if (!all(c("PC14", "KMS26") %in% rownames(dep_score))) {
    stop("ceres_used must contain PC14 and KMS26 rows.", call. = FALSE)
  }

  x <- as.data.frame(cys_dropout, stringsAsFactors = FALSE)
  x$abe_neg_log10_fdr <- cpt_cys_within_gene_fdr(x$p_A, x$HUGO)
  x$cbe_neg_log10_fdr <- cpt_cys_within_gene_fdr(x$p_C, x$HUGO)

  # Shared dependencies must reproduce in both cell lines.
  both_abe_inconsistent <-
    (x$PC_mean_LFC_A > -0.3 | x$KMS_mean_LFC_A > -0.3) & x$cell == "both"
  both_cbe_inconsistent <-
    (x$PC_mean_LFC_C > -0.3 | x$KMS_mean_LFC_C > -0.3) & x$cell == "both"
  x$abe_neg_log10_fdr[cpt_cys_true(both_abe_inconsistent)] <- 0
  x$cbe_neg_log10_fdr[cpt_cys_true(both_cbe_inconsistent)] <- 0

  if (isTRUE(apply_context_filter)) {
    pc_abe_inconsistent <- (x$PC_mean_LFC_A - x$KMS_mean_LFC_A) > -0.3
    pc_cbe_inconsistent <- (x$PC_mean_LFC_C - x$KMS_mean_LFC_C) > -0.3
    kms_abe_inconsistent <- (x$KMS_mean_LFC_A - x$PC_mean_LFC_A) > -0.3
    kms_cbe_inconsistent <- (x$KMS_mean_LFC_C - x$PC_mean_LFC_C) > -0.3

    dep_delta <- dep_score["PC14", ] - dep_score["KMS26", ]
    pc_selective_genes <- names(which(dep_delta <= -0.8))
    kms_selective_genes <- names(which(dep_delta >= 0.8))

    x$abe_neg_log10_fdr[cpt_cys_true(
      pc_abe_inconsistent & x$HUGO %in% pc_selective_genes & x$cell == "PC"
    )] <- 0
    x$cbe_neg_log10_fdr[cpt_cys_true(
      pc_cbe_inconsistent & x$HUGO %in% pc_selective_genes & x$cell == "PC"
    )] <- 0
    x$abe_neg_log10_fdr[cpt_cys_true(
      kms_abe_inconsistent & x$HUGO %in% kms_selective_genes & x$cell == "KMS"
    )] <- 0
    x$cbe_neg_log10_fdr[cpt_cys_true(
      kms_cbe_inconsistent & x$HUGO %in% kms_selective_genes & x$cell == "KMS"
    )] <- 0
  }

  # Exact thresholds from filter_dropout_data(): empirical p < 0.05,
  # within-gene FDR < 0.1, and mean dropout LFC <= -0.6.
  x$abe_functional <- cpt_cys_true(
    x$p_A > 1.3 & x$abe_neg_log10_fdr > 1 & x$mean_LFC_A <= -0.6
  )
  x$cbe_functional <- cpt_cys_true(
    x$p_C > 1.3 & x$cbe_neg_log10_fdr > 1 & x$mean_LFC_C <= -0.6
  )
  x$functional <- x$abe_functional | x$cbe_functional

  x$ligandability_score <- x$KB_engage
  x$ligandability_manual_annotation <- rep(FALSE, nrow(x))

  # Cys_plot.R explicitly adds this known covalent-drug site before plotting.
  if (isTRUE(apply_manual_annotations)) {
    egfr_c797 <- x$unique_ID == "EGFR_797"
    x$ligandability_score[egfr_c797] <- 100
    x$ligandability_manual_annotation[egfr_c797] <- TRUE
  }

  x$ligandable <- cpt_cys_true(x$ligandability_score > 50)
  x$functional_ligandable <- x$functional & x$ligandable

  af_key <- paste0(af_rsa$symbol, "_", af_rsa$residue_number)
  af_match <- match(x$unique_ID, af_key)
  sequence_match <- af_rsa$aa_same_as_portal[af_match]
  mapping_status <- ifelse(
    is.na(af_match),
    "Unmapped to source AlphaFold table",
    ifelse(sequence_match == "y", "Source sequence matched", "Source sequence mismatch")
  )

  study_context <- c(
    both = "PC14 and KMS26",
    PC = "PC14-selective",
    KMS = "KMS26-selective"
  )[x$cell]
  study_context[is.na(study_context)] <- x$cell[is.na(study_context)]

  editor_support <- ifelse(
    x$abe_functional & x$cbe_functional,
    "ABE and CBE",
    ifelse(x$abe_functional, "ABE", ifelse(x$cbe_functional, "CBE", "Not significant"))
  )

  atlas <- data.frame(
    site_id = x$unique_ID,
    gene_symbol = x$HUGO,
    protein_name = x$full_name,
    cysteine_position = suppressWarnings(as.integer(x$target_Cys)),
    uniprot_accession = af_rsa$uniprot_accession[af_match],
    residue_mapping_status = mapping_status,
    source_residue_is_cysteine = af_rsa$aa[af_match] == "C",
    study_context = unname(study_context),
    editor_support = editor_support,
    abe_mean_lfc = x$mean_LFC_A,
    cbe_mean_lfc = x$mean_LFC_C,
    abe_guide_count = x$A_guide_num,
    cbe_guide_count = x$C_guide_num,
    abe_neg_log10_p = x$p_A,
    cbe_neg_log10_p = x$p_C,
    abe_neg_log10_fdr = x$abe_neg_log10_fdr,
    cbe_neg_log10_fdr = x$cbe_neg_log10_fdr,
    abe_functional = x$abe_functional,
    cbe_functional = x$cbe_functional,
    functional = x$functional,
    ligandability_score = x$ligandability_score,
    ligandable = x$ligandable,
    functional_ligandable = x$functional_ligandable,
    ligandability_manual_annotation = x$ligandability_manual_annotation,
    proteomic_accessibility = x$IADTB_access,
    relative_solvent_accessibility = x$RSA,
    alphafold_plddt = af_rsa$pLDDT[af_match],
    pdb_coverage = af_rsa$covered_by_PDB[af_match],
    pdb_coverage_30 = af_rsa$covered_by_PDB_30[af_match],
    denaturation_lfc = x$denature_LFC,
    pc14_abe_mean_lfc = x$PC_mean_LFC_A,
    pc14_cbe_mean_lfc = x$PC_mean_LFC_C,
    kms26_abe_mean_lfc = x$KMS_mean_LFC_A,
    kms26_cbe_mean_lfc = x$KMS_mean_LFC_C,
    stringsAsFactors = FALSE
  )

  # --- Evolutionary conservation (Cys_conserve_score.Rdat) ---
  if (!is.null(conserve_score)) {
    cpt_cys_required_columns(conserve_score, c("score1", "score2"), "Cys_conserve_score")
    conserve_match <- match(atlas$site_id, rownames(conserve_score))
    atlas$conservation_score <- conserve_score$score2[conserve_match]
    atlas$conservation_raw <- conserve_score$score1[conserve_match]
  } else {
    atlas$conservation_score <- NA_real_
    atlas$conservation_raw <- NA_real_
  }

  # --- Ortholog cysteine counts (ortho_count_df_20230419.Rdat) ---
  if (!is.null(ortho_count)) {
    if (!"C" %in% colnames(ortho_count)) {
      stop("ortho_count must contain a 'C' column.", call. = FALSE)
    }
    ortho_match <- match(atlas$site_id, rownames(ortho_count))
    atlas$ortholog_cys_count <- ortho_count[["C"]][ortho_match]
    atlas$ortholog_total <- rowSums(ortho_count)[ortho_match]
  } else {
    atlas$ortholog_cys_count <- NA_integer_
    atlas$ortholog_total <- NA_integer_
  }

  # --- ClinVar pathogenic cysteine mutations (clinvar_data.Rdat) ---
  if (!is.null(clinvar)) {
    cpt_cys_required_columns(clinvar, c("HUGO", "pos", "from_aa", "annotation", "phenotype"), "clinvar_data")
    cys_clinvar <- clinvar[clinvar$from_aa == "C", , drop = FALSE]
    cys_clinvar$unique_ID <- paste0(cys_clinvar$HUGO, "_", cys_clinvar$pos)
    # Collapse multiple ClinVar entries per site (keep most severe annotation)
    severity_order <- c("Pathogenic", "Pathogenic/Likely pathogenic", "Likely pathogenic")
    cys_clinvar$severity_rank <- match(cys_clinvar$annotation, severity_order)
    cys_clinvar$severity_rank[is.na(cys_clinvar$severity_rank)] <- 99L
    cys_clinvar <- cys_clinvar[order(cys_clinvar$severity_rank), , drop = FALSE]
    cys_clinvar_dedup <- cys_clinvar[!duplicated(cys_clinvar$unique_ID), , drop = FALSE]
    clinvar_match <- match(atlas$site_id, cys_clinvar_dedup$unique_ID)
    atlas$clinvar_pathogenic <- !is.na(clinvar_match)
    atlas$clinvar_annotation <- cys_clinvar_dedup$annotation[clinvar_match]
    atlas$clinvar_phenotype <- cys_clinvar_dedup$phenotype[clinvar_match]
  } else {
    atlas$clinvar_pathogenic <- rep(FALSE, nrow(atlas))
    atlas$clinvar_annotation <- rep(NA_character_, nrow(atlas))
    atlas$clinvar_phenotype <- rep(NA_character_, nrow(atlas))
  }

  atlas <- atlas[order(!atlas$functional_ligandable, !atlas$functional,
                       atlas$gene_symbol, atlas$cysteine_position), , drop = FALSE]
  rownames(atlas) <- NULL

  attr(atlas, "source_repository") <- "https://github.com/cravattlab/Cys_editing"
  attr(atlas, "source_commit") <- source_commit
  attr(atlas, "source_object") <- "Part5_global_analysis/Rdat/Cys_dropout.Rdat"
  attr(atlas, "source_structure_object") <- "Part5_global_analysis/Rdat/AF_RSA.Rdat"
  attr(atlas, "source_conservation_object") <- "Part5_global_analysis/Rdat/Cys_conserve_score.Rdat"
  attr(atlas, "source_ortholog_object") <- "Part5_global_analysis/Rdat/ortho_count_df_20230419.Rdat"
  attr(atlas, "source_clinvar_object") <- "Part5_global_analysis/Rdat/clinvar_data.Rdat"
  attr(atlas, "source_license") <- "MIT (copyright 2023 Jason Li)"
  attr(atlas, "study_doi") <- "10.1038/s41589-023-01428-w"
  attr(atlas, "imported_at") <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  atlas
}

# Add Cys_editing evidence to a CanProTarget cysteine-probe table. The primary
# join is exact UniProt accession + cysteine position; gene + position is used
# only when the source AlphaFold mapping lacks an accession or an isoform label
# prevents an exact match.
cpt_cys_annotate_probe_table <- function(probe_data, atlas) {
  if (is.null(probe_data)) {
    return(NULL)
  }
  cpt_cys_required_columns(
    probe_data,
    c("gene_name", "cysteineid"),
    "probe_data"
  )
  cpt_cys_required_columns(
    atlas,
    c(
      "site_id", "gene_symbol", "cysteine_position", "uniprot_accession",
      "residue_mapping_status", "study_context", "editor_support",
      "abe_mean_lfc", "cbe_mean_lfc", "functional", "ligandable",
      "functional_ligandable"
    ),
    "cys_editing_atlas"
  )

  out <- as.data.frame(probe_data, stringsAsFactors = FALSE)
  cysteine_id <- as.character(out$cysteineid)
  position <- suppressWarnings(as.integer(sub(".*_C", "", cysteine_id)))
  accession <- sub("_C.*$", "", cysteine_id)
  accession <- sub("-.*$", "", accession)

  # Preserve the most informative record if a source table ever contains a
  # duplicate mapping key.
  atlas_order <- order(
    !cpt_cys_true(atlas$functional_ligandable),
    !cpt_cys_true(atlas$functional),
    atlas$site_id
  )
  atlas_use <- atlas[atlas_order, , drop = FALSE]
  exact_key <- paste(atlas_use$uniprot_accession, atlas_use$cysteine_position, sep = ":")
  gene_key <- paste(toupper(atlas_use$gene_symbol), atlas_use$cysteine_position, sep = ":")

  exact_match <- match(paste(accession, position, sep = ":"), exact_key)
  # Tolerate either gene-label convention. Probe tables use bare symbols, but a
  # table joined to dependency results carries DepMap's Entrez suffix
  # ("KRAS (3845)"); matching that raw against the atlas silently fails.
  gene_match <- match(
    paste(cpt_gene_match_key(out$gene_name), position, sep = ":"),
    gene_key
  )
  match_index <- exact_match
  fallback <- is.na(match_index) & !is.na(gene_match)
  match_index[fallback] <- gene_match[fallback]
  mapped <- !is.na(match_index)

  out$cys_site_id <- atlas_use$site_id[match_index]
  out$cys_site_in_atlas <- mapped
  out$cys_join_status <- ifelse(
    !mapped,
    "Not present in Cys_editing atlas",
    ifelse(fallback, "Gene + cysteine position", "UniProt + cysteine position")
  )
  out$cys_functional <- ifelse(mapped, cpt_cys_true(atlas_use$functional[match_index]), FALSE)
  out$cys_ligandable <- ifelse(mapped, cpt_cys_true(atlas_use$ligandable[match_index]), FALSE)
  out$cys_functional_ligandable <- ifelse(
    mapped,
    cpt_cys_true(atlas_use$functional_ligandable[match_index]),
    FALSE
  )
  out$cys_editor_support <- atlas_use$editor_support[match_index]
  out$cys_study_context <- atlas_use$study_context[match_index]
  out$cys_abe_mean_lfc <- atlas_use$abe_mean_lfc[match_index]
  out$cys_cbe_mean_lfc <- atlas_use$cbe_mean_lfc[match_index]
  out$cys_residue_mapping_status <- atlas_use$residue_mapping_status[match_index]
  out
}

# Expand dependency-linked probe rows to every functional Cys_editing site in
# the same gene. This supports target prioritisation even when the available
# chemoproteomic probe binds a different cysteine, while making that distinction
# explicit in cys_probe_site_match / cys_relationship.
cpt_cys_link_functional_sites <- function(probe_data, atlas) {
  if (is.null(probe_data)) {
    return(NULL)
  }
  cpt_cys_required_columns(probe_data, c("gene_name", "cysteineid"), "probe_data")
  cpt_cys_required_columns(
    atlas,
    c(
      "site_id", "gene_symbol", "cysteine_position", "uniprot_accession",
      "study_context", "editor_support", "abe_mean_lfc", "cbe_mean_lfc",
      "functional", "ligandable", "functional_ligandable",
      "residue_mapping_status"
    ),
    "cys_editing_atlas"
  )

  probe <- as.data.frame(probe_data, stringsAsFactors = FALSE)
  probe <- probe[, !startsWith(names(probe), "cys_"), drop = FALSE]
  functional_sites <- atlas[cpt_cys_true(atlas$functional), , drop = FALSE]
  functional_sites <- data.frame(
    gene_name = functional_sites$gene_symbol,
    cys_site_id = functional_sites$site_id,
    cys_cysteine_position = functional_sites$cysteine_position,
    cys_uniprot_accession = functional_sites$uniprot_accession,
    cys_editor_support = functional_sites$editor_support,
    cys_study_context = functional_sites$study_context,
    cys_abe_mean_lfc = functional_sites$abe_mean_lfc,
    cys_cbe_mean_lfc = functional_sites$cbe_mean_lfc,
    cys_ligandable = cpt_cys_true(functional_sites$ligandable),
    cys_functional_ligandable = cpt_cys_true(functional_sites$functional_ligandable),
    cys_residue_mapping_status = functional_sites$residue_mapping_status,
    stringsAsFactors = FALSE
  )

  # Join on the normalised key, not the raw label. `probe` may be a plain probe
  # table (bare symbols) or one already joined to dependency results, in which
  # case gene_name carries DepMap's Entrez suffix ("KRAS (3845)"). A raw join
  # silently returns zero rows for the latter. Normalising also makes the match
  # case-insensitive. The display gene_name from `probe` is preserved.
  probe$.cpt_gene_key <- cpt_gene_match_key(probe$gene_name)
  functional_sites$.cpt_gene_key <- cpt_gene_match_key(functional_sites$gene_name)
  functional_sites$gene_name <- NULL

  linked <- dplyr::inner_join(
    probe,
    functional_sites,
    by = ".cpt_gene_key",
    relationship = "many-to-many"
  )
  linked$.cpt_gene_key <- NULL
  probe_position <- suppressWarnings(as.integer(sub(".*_C", "", linked$cysteineid)))
  linked$cys_probe_site_match <- !is.na(probe_position) &
    probe_position == linked$cys_cysteine_position
  linked$cys_relationship <- ifelse(
    linked$cys_probe_site_match,
    "Probe binds this functional site",
    "Functional site in same dependency gene"
  )
  linked
}
