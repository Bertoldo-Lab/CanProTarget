# ============================================================
# Script:   usage_logger.R
# Purpose:  Optional append-only usage log, one JSON object per line,
#           under logs/usage/. Off unless CPT_LOG_USAGE=1.
# Exports:  cpt_log_usage()
# Note:     No-op when the variable is unset or jsonlite is missing, so
#           call sites never have to guard. What is recorded and what is
#           deliberately not: docs/LOGGING.md.
# ============================================================
# Sections:
#   cpt_log_usage()
# ============================================================


cpt_log_usage <- function(event, details = list()) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    return(invisible(NULL))
  }
  cfg <- if (exists("app_config", inherits = TRUE)) {
    app_config$logging
  } else {
    list()
  }
  env_on <- identical(Sys.getenv("CPT_LOG_USAGE", ""), "1")
  if (!env_on && !isTRUE(cfg$enabled)) {
    return(invisible(NULL))
  }

  log_dir <- if (!is.null(cfg$dir) && nzchar(as.character(cfg$dir)[1])) {
    cfg$dir
  } else {
    "logs/usage"
  }
  if (!dir.exists(log_dir)) {
    dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  }

  rec <- c(
    list(
      ts = format(Sys.time(), "%Y-%m-%dT%H:%M:%OSZ", tz = "UTC"),
      event = as.character(event)
    ),
    details
  )

  line <- tryCatch(
    jsonlite::toJSON(rec, auto_unbox = TRUE),
    error = function(e) {
      paste0('{"ts":"error","event":"log_encode_failed"}')
    }
  )

  fpath <- file.path(log_dir, paste0("usage_", Sys.Date(), ".jsonl"))
  tryCatch(
    cat(as.character(line), "\n", sep = "", file = fpath, append = TRUE),
    error = function(e) invisible(NULL)
  )
  invisible(NULL)
}
