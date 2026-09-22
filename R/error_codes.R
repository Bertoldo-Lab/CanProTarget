# ============================================================
# Script:   error_codes.R
# Purpose:  Centralized error handling with numbered codes.
#           Every user-facing error gets a CPT-XXXX code for
#           easy identification in bug reports and logs.
# ============================================================
# Sections:
#   cpt_error_msg() — Format a user-facing error notification with a CPT error code
#   cpt_notify_error() — Show a numbered error notification in Shiny
# ============================================================

# Error code registry:
#
# CPT-1000 series: Data loading errors
#   CPT-1001  Cancer model metadata file missing
#   CPT-1002  Cancer model metadata file corrupt/unreadable
#   CPT-1003  Gene effect matrix file missing
#   CPT-1004  Gene effect matrix file corrupt/unreadable
#   CPT-1005  Cysteine editing atlas file missing
#   CPT-1006  Cysteine editing atlas file corrupt/unreadable
#   CPT-1007  Protein binding lookup file missing
#   CPT-1008  Protein binding lookup file corrupt/unreadable
#   CPT-1009  SwissADME data file missing or corrupt
#   CPT-1010  ADME per-gene lookup could not be built from binding + SwissADME
#   CPT-1010  Data version manifest missing
#
# CPT-2000 series: Analysis / computation errors
#   CPT-2001  Gene not found in dataset
#   CPT-2002  Subtype not found in metadata
#   CPT-2003  Insufficient cell lines for analysis
#   CPT-2004  Analysis computation failed (limma/t-test)
#   CPT-2005  Precomputed file read error
#
# CPT-3000 series: Input validation errors
#   CPT-3001  Invalid gene symbol format
#   CPT-3002  Invalid subtype selection
#   CPT-3003  No genes match current filters
#   CPT-3004  Non-character input where text expected
#   CPT-3005  Empty or NULL input
#
# CPT-4000 series: Module/UI errors
#   CPT-4001  Cysteine atlas unavailable for filtering
#   CPT-4002  Probe data unavailable
#   CPT-4003  Render failed (plot/table)
#   CPT-4004  Download handler failed
#
# CPT-5000 series: MCP/API errors
#   CPT-5001  R worker not responding
#   CPT-5002  Tool dispatch failed
#   CPT-5003  Invalid tool parameters
#   CPT-5004  Query timeout

#' Format a user-facing error notification with a CPT error code.
#' @param code Character: error code (e.g. "CPT-3004")
#' @param message Character: human-readable error description
#' @param detail Character: optional technical detail (logged, not always shown)
#' @return Character string formatted for showNotification
cpt_error_msg <- function(code, message, detail = NULL) {
  msg <- paste0("[", code, "] ", message)
  if (!is.null(detail) && nzchar(detail)) {
    msg <- paste0(msg, " (", detail, ")")
  }
  msg
}

#' Show a numbered error notification in Shiny.
#' @param code Character: error code
#' @param message Character: user-facing message
#' @param detail Character: optional detail
#' @param session Shiny session (NULL for non-Shiny context)
#' @param type Character: "error", "warning", or "message"
#' @param duration Integer: seconds to show (NULL for sticky)
cpt_notify_error <- function(code, message, detail = NULL,
                             session = NULL, type = "error", duration = 10) {
  full_msg <- cpt_error_msg(code, message, detail)

  # Log to stderr regardless
  cat(paste0("[", Sys.time(), "] ", full_msg, "\n"), file = stderr())

  # Show Shiny notification if in reactive context
  if (!is.null(session) || shiny::isRunning()) {
    shiny::showNotification(full_msg, type = type, duration = duration)
  }

  invisible(full_msg)
}
