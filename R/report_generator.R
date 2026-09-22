# ============================================================
# Script:   report_generator.R
# Purpose:  Report rendering engine for CanProTarget.
#           Generates branded HTML reports from Rmd templates.
#           Called by: Shiny download handlers, MCP generate_report tool.
#           No Shiny dependencies; pure R + rmarkdown.
#
# Report types:
#   - gene_dependency: Full gene analysis in a cancer subtype
#   - cysteine_target: Cysteine-level annotation for a gene
#
# Usage:
#   cpt_render_report("gene_dependency", list(gene = "KRAS", subtype = "..."))
#   cpt_render_report("cysteine_target", list(gene = "EGFR", site_id = "EGFR_797"))
# ============================================================
# Sections:
#   Template discovery
#     cpt_template_dir(), cpt_available_reports()
#   Rendering
#     cpt_render_report()
#   Convenience wrappers for each report type
#     cpt_report_gene_dependency(), cpt_report_cysteine_target()
#   MCP API wrapper
#     api_generate_report()
# ============================================================

# ---- Template discovery ----------------------------------------

#' Find the report templates directory.
#' Searches: inst/report_templates/ relative to project root.
#' @param project_root Character: project root directory
#' @return Path to templates directory
cpt_template_dir <- function(project_root = ".") {

  path <- file.path(project_root, "inst", "report_templates")
  if (!dir.exists(path)) {
    stop("Report templates directory not found: ", path, call. = FALSE)
  }
  normalizePath(path)
}

#' List available report types.
#' @param project_root Character: project root directory
#' @return Character vector of report type names
cpt_available_reports <- function(project_root = ".") {
  tpl_dir <- cpt_template_dir(project_root)
  rmds <- list.files(tpl_dir, pattern = "\\.Rmd$", full.names = FALSE)
  sub("_report\\.Rmd$", "", rmds)
}

# ---- Rendering -------------------------------------------------

#' Render a CanProTarget report to HTML.
#'
#' @param report_type Character: one of "gene_dependency", "cysteine_target"
#' @param params List: parameters for the report (gene, subtype, dataset, site_id, etc.)
#' @param output_file Character: output filename. If NULL, auto-generated in tempdir().
#' @param output_dir Character: directory for output. Defaults to tempdir().
#' @param project_root Character: project root (for finding templates and data).
#' @return Path to rendered HTML file (invisible).
cpt_render_report <- function(report_type,
                              params = list(),
                              output_file = NULL,
                              output_dir = NULL,
                              project_root = ".") {

  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    stop("Package 'rmarkdown' is required for report generation.", call. = FALSE)
  }

  # Validate report type
  tpl_dir <- cpt_template_dir(project_root)
  tpl_file <- file.path(tpl_dir, paste0(report_type, "_report.Rmd"))
  if (!file.exists(tpl_file)) {
    available <- cpt_available_reports(project_root)
    stop("Unknown report type '", report_type, "'. Available: ",
         paste(available, collapse = ", "), call. = FALSE)
  }

  # Set defaults
  if (is.null(output_dir)) output_dir <- tempdir()
  if (is.null(output_file)) {
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    safe_name <- gsub("[^A-Za-z0-9_]", "_", paste0(report_type, "_",
                      params$gene %||% "report"))
    output_file <- paste0("CanProTarget_", safe_name, "_", timestamp, ".html")
  }

  # Ensure output_dir exists
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # Inject project_root and data_dir into params
  params$project_root <- normalizePath(project_root)
  if (is.null(params$data_dir)) params$data_dir <- "data"

  # Copy CSS to temp render location (rmarkdown resolves relative to .Rmd)
  css_src <- file.path(tpl_dir, "report_style.css")
  css_dest <- file.path(output_dir, "report_style.css")
  if (file.exists(css_src) && !file.exists(css_dest)) {
    file.copy(css_src, css_dest, overwrite = TRUE)
  }

  # Copy template to output_dir so relative CSS path works
  tpl_copy <- file.path(output_dir, basename(tpl_file))
  file.copy(tpl_file, tpl_copy, overwrite = TRUE)

  # Render
  out_path <- rmarkdown::render(
    input = tpl_copy,
    output_file = output_file,
    output_dir = output_dir,
    params = params,
    envir = new.env(parent = globalenv()),
    quiet = TRUE
  )

  # Clean up copied template
  unlink(tpl_copy)

  invisible(normalizePath(out_path))
}

# ---- Convenience wrappers for each report type -----------------

#' Generate a Gene Dependency Report.
#' @param gene Character: gene symbol
#' @param subtype Character: OncotreeSubtype
#' @param dataset Character: "CRISPR" or "RNAi"
#' @param output_file Character: filename (optional)
#' @param output_dir Character: directory (optional)
#' @param project_root Character: project root
#' @return Path to rendered HTML
cpt_report_gene_dependency <- function(gene, subtype, dataset = "CRISPR",
                                       output_file = NULL, output_dir = NULL,
                                       project_root = ".") {
  cpt_render_report(
    report_type = "gene_dependency",
    params = list(gene = gene, subtype = subtype, dataset = dataset),
    output_file = output_file,
    output_dir = output_dir,
    project_root = project_root
  )
}

#' Generate a Cysteine Target Report.
#' @param gene Character: gene symbol
#' @param site_id Character: specific site ID (e.g. "EGFR_797"), or NULL for all gene sites
#' @param output_file Character: filename (optional)
#' @param output_dir Character: directory (optional)
#' @param project_root Character: project root
#' @return Path to rendered HTML
cpt_report_cysteine_target <- function(gene, site_id = NULL,
                                       output_file = NULL, output_dir = NULL,
                                       project_root = ".") {
  cpt_render_report(
    report_type = "cysteine_target",
    params = list(gene = gene, site_id = site_id %||% ""),
    output_file = output_file,
    output_dir = output_dir,
    project_root = project_root
  )
}

# ---- MCP API wrapper -------------------------------------------

#' Generate a report via MCP. Returns the file path to the rendered HTML.
#' @param report_type Character: "gene_dependency" or "cysteine_target"
#' @param params List: report parameters
#' @param data_env List: loaded data environment (from mcp_worker)
#' @return List with status, file path, and report metadata
api_generate_report <- function(report_type, params, data_env) {
  project_root <- dirname(data_env$data_dir)

  # Validate required params per report type
  if (report_type == "gene_dependency") {
    if (is.null(params$gene)) stop("Parameter 'gene' is required.", call. = FALSE)
    if (is.null(params$subtype)) stop("Parameter 'subtype' is required.", call. = FALSE)
    params$dataset <- params$dataset %||% "CRISPR"
  } else if (report_type == "cysteine_target") {
    if (is.null(params$gene)) stop("Parameter 'gene' is required.", call. = FALSE)
  } else {
    available <- cpt_available_reports(project_root)
    stop("Unknown report_type '", report_type, "'. Available: ",
         paste(available, collapse = ", "), call. = FALSE)
  }

  # Render to a persistent output directory
  output_dir <- file.path(project_root, "reports")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  out_path <- cpt_render_report(
    report_type = report_type,
    params = params,
    output_dir = output_dir,
    project_root = project_root
  )

  list(
    status = "success",
    report_type = report_type,
    file_path = out_path,
    file_name = basename(out_path),
    parameters = params,
    generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
    citation = CPT_CITATION,
    source = CPT_REPO
  )
}

# Null-coalesce (safe redefinition)
if (!exists("%||%", mode = "function")) {
  `%||%` <- function(x, y) if (is.null(x)) y else x
}
