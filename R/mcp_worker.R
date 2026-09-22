#!/usr/bin/env Rscript --vanilla
# ============================================================
# Script:   mcp_worker.R
# Purpose:  Warm R process for MCP server. Loads data once at startup,
#           then reads JSON requests from stdin and writes JSON responses
#           to stdout. Communicates with Python MCP server via r_bridge.py.
# Usage:    Rscript --vanilla R/mcp_worker.R
# Protocol: One JSON object per line on stdin, one JSON response per line on stdout.
#           Stderr is used for logging (not part of the protocol).
# ============================================================
# Sections:
#   log_msg()
#   dispatch()
# ============================================================

# Only run as a worker process when executed via Rscript, not when sourced
# (e.g. accidental source, or older Shiny loadSupport without _disable_autoload).
if (sys.nframe() != 0L) {
  # Being sourced — export nothing; do not load data or block on stdin.
  invisible(NULL)
} else {

# --- Startup: load dependencies and data ------------------------

suppressPackageStartupMessages({
  library(jsonlite)
})

# Determine project root (this script is at R/mcp_worker.R)
script_dir <- if (interactive()) {
  "."
} else {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    dirname(dirname(sub("^--file=", "", file_arg[1])))
  } else {
    "."
  }
}

# Allow override via environment variable
project_root <- Sys.getenv("CPT_PROJECT_ROOT", unset = script_dir)
setwd(project_root)

# Source API functions
source(file.path(project_root, "R", "api_functions.R"))
source(file.path(project_root, "R", "canprotarget_score.R"))
source(file.path(project_root, "R", "pancancer_profile.R"))
source(file.path(project_root, "R", "report_generator.R"))

# Log to stderr (won't interfere with JSON protocol on stdout)
log_msg <- function(...) {
  msg <- paste0("[mcp_worker] ", paste(..., collapse = " "), "\n")
  cat(msg, file = stderr())
}

log_msg("Starting up. Project root:", project_root)
log_msg("Loading data files...")

# --- Load data into environment ---------------------------------

data_env <- list(
  data_dir = file.path(project_root, "data")
)

# CRISPR matrix
crispr_path <- file.path(data_env$data_dir, "CRISPRGeneEffect_23Q4_clean.rds")
if (file.exists(crispr_path)) {
  data_env$crispr_matrix <- readRDS(crispr_path)
  log_msg("Loaded CRISPR matrix:", nrow(data_env$crispr_matrix), "x", ncol(data_env$crispr_matrix))
} else {
  data_env$crispr_matrix <- NULL
  log_msg("WARNING: CRISPR matrix not found at", crispr_path)
}

# RNAi matrix
rnai_path <- file.path(data_env$data_dir, "d2_gene_effect_headers_refined.rds")
if (file.exists(rnai_path)) {
  data_env$rnai_matrix <- readRDS(rnai_path)
  log_msg("Loaded RNAi matrix:", nrow(data_env$rnai_matrix), "x", ncol(data_env$rnai_matrix))
} else {
  data_env$rnai_matrix <- NULL
  log_msg("WARNING: RNAi matrix not found at", rnai_path)
}

# Cancer model metadata
meta_path <- file.path(data_env$data_dir, "cancer_model_data.rds")
if (file.exists(meta_path)) {
  data_env$cancer_model_data <- readRDS(meta_path)
  log_msg("Loaded cancer model data:", nrow(data_env$cancer_model_data), "rows")
} else {
  data_env$cancer_model_data <- NULL
  log_msg("WARNING: Cancer model data not found at", meta_path)
}

# Cysteine editing atlas
cys_path <- file.path(data_env$data_dir, "cys_editing_atlas.rds")
if (file.exists(cys_path)) {
  data_env$cys_atlas <- readRDS(cys_path)
  log_msg("Loaded cysteine atlas:", nrow(data_env$cys_atlas), "sites")
} else {
  data_env$cys_atlas <- NULL
  log_msg("WARNING: Cysteine atlas not found at", cys_path)
}

# Per-gene ADME lookup + compact SMCL index from the ~10.6M-row binding table.
# Keep only CR >= 4 rows (~30k) for residue/SMCL queries; drop the full table.
adme_path <- file.path(data_env$data_dir, "swissadme_preprocessed.rds")
bind_path <- file.path(data_env$data_dir, "protein_binding_lookup_factored.rds")
if (!file.exists(bind_path)) {
  bind_path <- file.path(data_env$data_dir, "protein_binding_lookup_preprocessed.rds")
}
data_env$adme_gene_scores <- NULL
data_env$smcl_index <- NULL
if (file.exists(bind_path)) {
  ok <- tryCatch({
    .bind <- readRDS(bind_path)
    data_env$smcl_index <- cpt_build_smcl_index(.bind, min_cr = 4)
    if (file.exists(adme_path)) {
      .adme <- readRDS(adme_path)
      data_env$adme_gene_scores <- cpt_build_adme_gene_scores(.bind, .adme)
      rm(.adme)
    }
    rm(.bind)
    gc(verbose = FALSE)
    TRUE
  }, error = function(e) {
    log_msg("WARNING: chemoproteomics index build failed:", conditionMessage(e))
    FALSE
  })
  if (isTRUE(ok)) {
    if (!is.null(data_env$smcl_index)) {
      log_msg("Built SMCL index:", nrow(data_env$smcl_index), "CR>=4 rows")
    }
    if (!is.null(data_env$adme_gene_scores)) {
      log_msg("Built ADME gene lookup:", nrow(data_env$adme_gene_scores), "genes")
    }
  }
} else {
  log_msg("Chemoproteomics RDS absent; SMCL/ADME index skipped")
}

log_msg("Data loading complete. Ready for requests.")

# --- Dispatch table ---------------------------------------------

dispatch <- function(request) {
  tool <- request$tool
  params <- request$params
  if (is.null(params)) params <- list()

  result <- switch(tool,
    "list_subtypes" = {
      dataset <- params$dataset %||% "CRISPR"
      api_list_subtypes(dataset, data_env)
    },
    "list_genes" = {
      dataset <- params$dataset %||% "CRISPR"
      api_list_genes(dataset, data_env)
    },
    "query_dependency" = {
      api_query_dependency(
        gene = params$gene,
        subtype = params$subtype,
        dataset = params$dataset %||% "CRISPR",
        data_env = data_env
      )
    },
    "top_dependencies" = {
      api_top_dependencies(
        subtype = params$subtype,
        dataset = params$dataset %||% "CRISPR",
        n = params$n %||% 20L,
        data_env = data_env
      )
    },
    "gene_cysteines" = {
      api_gene_cysteines(gene = params$gene, data_env = data_env)
    },
    "cysteine_detail" = {
      api_cysteine_detail(site_id = params$site_id, data_env = data_env)
    },
    "compare_subtypes" = {
      api_compare_subtypes(
        gene = params$gene,
        subtype1 = params$subtype1,
        subtype2 = params$subtype2,
        dataset = params$dataset %||% "CRISPR",
        data_env = data_env
      )
    },
    "platform_info" = {
      api_platform_info()
    },
    "canprotarget_score" = {
      api_canprotarget_score(
        gene = params$gene,
        subtype = params$subtype,
        dataset = params$dataset %||% "CRISPR",
        data_env = data_env
      )
    },
    "pancancer_profile" = {
      api_pancancer_profile(
        gene = params$gene,
        dataset = params$dataset %||% "CRISPR",
        data_env = data_env,
        top_n = params$top_n %||% 30L
      )
    },
    "generate_report" = {
      api_generate_report(
        report_type = params$report_type,
        params = params,
        data_env = data_env
      )
    },
    "rank_targets" = {
      api_rank_targets(
        subtype = params$subtype,
        dataset = params$dataset %||% "CRISPR",
        n = params$n %||% 20L,
        pool = params$pool %||% 100L,
        require_ligandable = isTRUE(params$require_ligandable),
        data_env = data_env
      )
    },
    "rank_site_targets" = {
      api_rank_site_targets(
        subtype = params$subtype,
        dataset = params$dataset %||% "CRISPR",
        n = params$n %||% 20L,
        pool = params$pool %||% 100L,
        min_cr = params$min_cr %||% 4,
        max_targets = params$max_targets %||% 20,
        effect_size_max = params$effect_size_max,
        p_max = params$p_max,
        exclude_common_essentials = isTRUE(params$exclude_common_essentials),
        data_env = data_env
      )
    },
    "assess_target" = {
      api_assess_target(
        gene = params$gene,
        subtype = params$subtype,
        dataset = params$dataset %||% "CRISPR",
        data_env = data_env
      )
    },
    "ping" = {
      list(status = "ok", timestamp = format(Sys.time(), tz = "UTC"))
    },
    {
      stop("Unknown tool: ", tool, call. = FALSE)
    }
  )
  result
}

# Null-coalesce operator (R < 4.4 compat)
`%||%` <- function(x, y) if (is.null(x)) y else x

# --- Main loop: read JSON from stdin, write JSON to stdout ------

# Signal readiness
cat(toJSON(list(status = "ready", pid = Sys.getpid()), auto_unbox = TRUE), "\n", sep = "")
flush(stdout())

log_msg("Entering main loop (pid:", Sys.getpid(), ")")

con <- file("stdin", "r")

repeat {
  line <- readLines(con, n = 1L, warn = FALSE)

  # EOF — parent process closed stdin

  if (length(line) == 0L) {
    log_msg("stdin closed, shutting down.")
    break
  }

  # Skip empty lines
  line <- trimws(line)
  if (!nzchar(line)) next

  # Parse request
  request <- tryCatch(
    fromJSON(line, simplifyVector = TRUE, simplifyDataFrame = FALSE),
    error = function(e) {
      list(.parse_error = conditionMessage(e))
    }
  )

  if (!is.null(request$.parse_error)) {
    response <- list(
      error = TRUE,
      message = paste("JSON parse error:", request$.parse_error)
    )
  } else if (is.null(request$tool)) {
    response <- list(
      error = TRUE,
      message = "Request missing 'tool' field."
    )
  } else {
    # Execute
    response <- tryCatch(
      {
        result <- dispatch(request)
        list(error = FALSE, result = result)
      },
      error = function(e) {
        log_msg("Error in tool '", request$tool, "':", conditionMessage(e))
        list(error = TRUE, message = conditionMessage(e))
      }
    )
  }

  # Write response (single line JSON)
  json_out <- toJSON(response, auto_unbox = TRUE, null = "null", na = "null")
  cat(json_out, "\n", sep = "")
  flush(stdout())
}

log_msg("Exiting.")

} # end if (sys.nframe() == 0L) main worker
