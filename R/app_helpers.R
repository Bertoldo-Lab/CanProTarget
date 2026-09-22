# ============================================================
# Script:   app_helpers.R
# Purpose:  Small helpers shared across the app: path building, data
#           signatures for bindCache, the status banner and the
#           data-unavailable card, and the row indexes that keep probe
#           lookups off a 10.6M-row scan.
# Exports:  app_data_path(), cpt_is_true(), *_sig(), app_status_banner(),
#           cpt_data_unavailable_card(), cpt_index_protein_binding(),
#           cpt_pb_subset(), read_data_versions()
# ============================================================
# Sections:
#   Data versioning manifest
#     read_data_versions(), data_provenance_string(), data_versions_table(), cpt_data_versions_ui()
# ============================================================

# Null-coalesce operator (available natively in R >= 4.4, defined here for safety)
if (!exists("%||%", mode = "function")) {
  `%||%` <- function(x, y) if (is.null(x)) y else x
}

app_data_path <- function(data_dir, ...) {
  file.path(data_dir, ...)
}

#' Vectorised truthiness across the mixed conventions in these datasets.
#'
#' The cysteine atlas uses real logicals; the binding table uses "yes"/NA
#' strings. Comparing either with `== TRUE` gives NA-riddled results, so all
#' flag columns go through here.
#'
#' @param x logical or character vector
#' @return logical vector, never NA
cpt_is_true <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  v <- tolower(trimws(as.character(x)))
  !is.na(v) & v %in% c("true", "yes", "y", "1")
}

swissadme_data_sig <- function(data_dir) {
  p <- file.path(data_dir, "swissadme_preprocessed.rds")
  if (file.exists(p)) file.mtime(p) else 0
}

cancer_model_data_sig <- function(data_dir) {
  p <- app_data_path(data_dir, "cancer_model_data.rds")
  if (file.exists(p)) file.mtime(p) else 0
}

app_status_banner <- function(..., type = c("info", "warning")) {
  type <- match.arg(type)
  cls <- if (identical(type, "warning")) "cpt-status-banner cpt-status-warning" else "cpt-status-banner cpt-status-info"
  shiny::tags$div(class = cls, ...)
}

#' Placeholder card when an optional dataset is missing (Protein Lookup / SwissADME).
#' @param title Short heading
#' @param detail What the user can do next
#' @param code Optional CPT error code string
cpt_data_unavailable_card <- function(title, detail, code = NULL) {
  shiny::tags$div(
    class = "alert alert-warning",
    style = paste(
      "margin: 12px 0; padding: 16px 18px; border-left: 4px solid #e0a800;",
      "background: #fff8e6; color: #4a3b00;"
    ),
    shiny::tags$h4(style = "margin-top: 0;", title),
    shiny::tags$p(detail),
    if (!is.null(code) && nzchar(code)) {
      shiny::tags$p(
        style = "margin-bottom: 0; font-size: 0.9em; opacity: 0.85;",
        code
      )
    },
    shiny::tags$p(
      style = "margin: 10px 0 0 0; font-size: 0.9em;",
      "See ", shiny::tags$code("docs/DATA_PROVENANCE.md"), " for how these files are built."
    )
  )
}

#' Suggest close gene/protein symbols (for typos).
#' @param query Character vector of user-entered names
#' @param candidates Character vector of valid symbols
#' @param n Max suggestions per query
#' @return Named list: query -> character vector of suggestions (may be empty)
cpt_suggest_symbols <- function(query, candidates, n = 5) {
  if (!length(query) || !length(candidates)) {
    return(stats::setNames(vector("list", length(query)), query))
  }
  cand <- unique(as.character(candidates))
  cand_l <- tolower(cand)
  out <- lapply(query, function(q) {
    ql <- tolower(trimws(as.character(q)))
    if (!nzchar(ql)) return(character(0))
    exact <- cand[cand_l == ql]
    if (length(exact)) return(exact[1])
    # Prefix then fuzzy
    pref <- cand[startsWith(cand_l, ql) | startsWith(ql, cand_l)]
    if (length(pref) >= n) return(utils::head(pref, n))
    fuzzy_idx <- agrep(ql, cand_l, max.distance = 0.2, value = FALSE)
    fuzzy <- if (length(fuzzy_idx)) cand[fuzzy_idx] else character(0)
    utils::head(unique(c(pref, fuzzy)), n)
  })
  stats::setNames(out, query)
}

# Row-index maps for the ~10M-row protein-binding table. Built once on RDS load
# so the probe lookups do not scan every row.
cpt_index_protein_binding <- function(lookup) {
  if (is.null(lookup) || !is.data.frame(lookup) || !nrow(lookup)) {
    return(lookup)
  }
  n <- nrow(lookup)
  if ("gene_name_key" %in% names(lookup) && is.null(attr(lookup, "rows_by_gene"))) {
    attr(lookup, "rows_by_gene") <- split(seq_len(n), lookup$gene_name_key, drop = TRUE)
  }
  if ("probe_name" %in% names(lookup) && is.null(attr(lookup, "rows_by_probe"))) {
    attr(lookup, "rows_by_probe") <- split(
      seq_len(n),
      tolower(as.character(lookup$probe_name)),
      drop = TRUE
    )
  }
  lookup
}

cpt_pb_subset <- function(lookup, keys, index = c("rows_by_gene", "rows_by_probe")) {
  index <- match.arg(index)
  if (is.null(lookup) || !is.data.frame(lookup)) return(lookup)
  keys <- unique(as.character(keys))
  keys <- keys[!is.na(keys) & nzchar(keys)]
  if (!length(keys) || !nrow(lookup)) {
    return(lookup[0, , drop = FALSE])
  }
  idx <- attr(lookup, index)
  if (is.null(idx)) return(NULL)
  rows <- unlist(idx[keys], use.names = FALSE)
  if (!length(rows)) return(lookup[0, , drop = FALSE])
  lookup[rows, , drop = FALSE]
}


# ---- Data versioning manifest ----------------------------------

#' Read data_versions.yaml and return as a list.
#' Returns NULL if yaml package is unavailable or file missing.
read_data_versions <- function(data_dir) {
  path <- file.path(data_dir, "data_versions.yaml")
  if (!file.exists(path)) return(NULL)
  if (!requireNamespace("yaml", quietly = TRUE)) return(NULL)
  tryCatch(
    yaml::yaml.load_file(path),
    error = function(e) NULL
  )
}

#' Format a one-line data provenance string for reports.
#' e.g. "DepMap 23Q4, Cys Atlas v1.0, SwissADME accessed 2024-01-15"
data_provenance_string <- function(data_versions) {
  if (is.null(data_versions)) return("Data versions not available")
  datasets <- data_versions$datasets
  if (is.null(datasets)) return("Data versions not available")

  parts <- vapply(datasets, function(d) {
    if (identical(d$status, "missing")) return(NA_character_)
    paste0(d$name, " ", d$version)
  }, character(1), USE.NAMES = FALSE)
  parts <- parts[!is.na(parts) & nzchar(parts)]

  if (!length(parts)) return("No datasets loaded")
  paste(parts, collapse = ", ")
}

#' Flatten data_versions.yaml datasets into a display data.frame.
#' @param data_versions List from read_data_versions()
#' @return data.frame with name, version, source, status, description, url
data_versions_table <- function(data_versions) {
  if (is.null(data_versions) || is.null(data_versions$datasets)) {
    return(data.frame(
      name = character(0), version = character(0), source = character(0),
      status = character(0), description = character(0), url = character(0),
      stringsAsFactors = FALSE
    ))
  }
  rows <- lapply(data_versions$datasets, function(d) {
    st <- if (!is.null(d$status) && nzchar(as.character(d$status)[1])) {
      as.character(d$status)[1]
    } else {
      "available"
    }
    data.frame(
      name = as.character(d$name %||% ""),
      version = as.character(d$version %||% ""),
      source = as.character(d$source %||% ""),
      status = st,
      description = as.character(d$description %||% ""),
      url = as.character(d$url %||% ""),
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(rows)
}

#' Shiny UI for the About-tab data versions table.
#' @param data_versions List from read_data_versions()
#' @return shiny.tag
cpt_data_versions_ui <- function(data_versions) {
  tbl <- data_versions_table(data_versions)
  platform <- data_versions$platform
  header <- if (!is.null(platform)) {
    shiny::tags$p(
      style = "margin-bottom: 12px;",
      shiny::tags$strong(platform$name %||% "CanProTarget"),
      " v", platform$version %||% "?",
      " · data refresh ", platform$last_data_refresh %||% "unknown",
      " · ", data_provenance_string(data_versions)
    )
  } else {
    shiny::tags$p(data_provenance_string(data_versions))
  }

  if (!nrow(tbl)) {
    return(shiny::tagList(
      header,
      shiny::tags$p(
        class = "text-muted",
        "No data_versions.yaml found. Add data/data_versions.yaml for provenance."
      )
    ))
  }

  body_rows <- lapply(seq_len(nrow(tbl)), function(i) {
    src <- if (nzchar(tbl$url[i])) {
      shiny::tags$a(href = tbl$url[i], target = "_blank", tbl$source[i])
    } else {
      tbl$source[i]
    }
    st_style <- if (identical(tbl$status[i], "missing")) {
      "padding: 8px; color: #a67c00; font-weight: 600;"
    } else {
      "padding: 8px; color: #2d6a4f;"
    }
    shiny::tags$tr(
      style = "border-bottom: 1px solid #ddd;",
      shiny::tags$td(style = "padding: 8px;", tbl$name[i]),
      shiny::tags$td(style = "padding: 8px;", tbl$version[i]),
      shiny::tags$td(style = "padding: 8px;", src),
      shiny::tags$td(style = st_style, tbl$status[i]),
      shiny::tags$td(style = "padding: 8px;", tbl$description[i])
    )
  })

  shiny::tagList(
    header,
    shiny::tags$table(
      style = "width: 100%; border-collapse: collapse; margin-bottom: 20px;",
      shiny::tags$thead(
        shiny::tags$tr(
          style = "border-bottom: 2px solid #478EB8;",
          shiny::tags$th(style = "padding: 8px; text-align: left;", "Dataset"),
          shiny::tags$th(style = "padding: 8px; text-align: left;", "Version"),
          shiny::tags$th(style = "padding: 8px; text-align: left;", "Source"),
          shiny::tags$th(style = "padding: 8px; text-align: left;", "Status"),
          shiny::tags$th(style = "padding: 8px; text-align: left;", "Description")
        )
      ),
      shiny::tags$tbody(body_rows)
    )
  )
}
