# ============================================================
# Script:   ui_kit.R
# Purpose:  Shared UI primitives so explanatory text has ONE home
#           and one behaviour across every tab.
#
#           Guidance collapses behind a click, so the default view is
#           results and the explanation is one icon away.
#
#           Built on <details>/<summary>: no JS, no popover init, works
#           inside DT cells, box headers and conditionalPanels alike.
# ============================================================
# Sections:
#   cpt_help() — Inline help icon that expands to a short explanation
#   cpt_note() — Collapsible note block with a visible one-line summary
#   cpt_section() — Section heading with optional help icon attached
#   cpt_label()
#   cpt_empty_state() — Empty / pre-run state, identical everywhere
#   cpt_filter_levels() — Make DT's per-column filters usable on categorical columns
# ============================================================

#' Inline help icon that expands to a short explanation.
#'
#' Use in place of a permanent tags$small() under a control or heading.
#'
#' @param ... Help content (text or tags)
#' @param label Optional visible text next to the icon
#' @return shiny.tag
cpt_help <- function(..., label = NULL) {
  shiny::tags$details(
    class = "cpt-help",
    shiny::tags$summary(
      title = "What is this?",
      shiny::icon("circle-info"),
      if (!is.null(label)) shiny::tags$span(class = "cpt-help-label", label)
    ),
    shiny::tags$div(class = "cpt-help-body", ...)
  )
}

#' Collapsible note block with a visible one-line summary.
#'
#' For caveats that must stay discoverable but should not occupy the
#' control panel permanently (e.g. the provisional ADME formula).
#'
#' @param summary_text One-line summary shown when collapsed
#' @param ... Body content, revealed on click
#' @param tone "caution" (amber) or "info" (blue)
#' @param open Start expanded?
#' @return shiny.tag
cpt_note <- function(summary_text, ..., tone = c("caution", "info"), open = FALSE) {
  tone <- match.arg(tone)
  shiny::tags$details(
    class = paste0("cpt-note cpt-note-", tone),
    open = if (isTRUE(open)) NA else NULL,
    shiny::tags$summary(summary_text),
    shiny::tags$div(class = "cpt-note-body", ...)
  )
}

#' Section heading with optional help icon attached.
#'
#' @param title Heading text
#' @param ... Help content; when supplied, rendered behind an icon
#' @param level Heading level ("h4" default)
#' @return shiny.tag
cpt_section <- function(title, ..., level = "h4") {
  help <- list(...)
  shiny::tags[[level]](
    class = "cpt-section-title",
    title,
    if (length(help)) do.call(cpt_help, help) else NULL
  )
}

#' Input label with the help icon beside the term, not floating beneath it.
#'
#' Pair with an input whose own label is NULL, so the icon sits on the label
#' line instead of becoming a stray control below the field.
#'
#' @param text Label text
#' @param ... Help content; omitted means a plain label
#' @param for_id Input id the label belongs to
#' @return shiny.tag
cpt_label <- function(text, ..., for_id = NULL) {
  help <- list(...)
  shiny::tags$label(
    class = "cpt-field-label",
    `for` = for_id,
    text,
    if (length(help)) do.call(cpt_help, help) else NULL
  )
}

#' Empty / pre-run state, identical everywhere.
#'
#' @param ... Message content
#' @param action Optional name of the button the user should press
#' @return shiny.tag
cpt_empty_state <- function(..., action = NULL) {
  shiny::tags$div(
    class = "cpt-empty-state",
    shiny::icon("arrow-left"),
    shiny::tags$span(...),
    if (!is.null(action)) {
      shiny::tags$b(paste0(" ", action))
    }
  )
}

#' Make DT's per-column filters usable on categorical columns.
#'
#' DT renders a free-text box for a character column and a dropdown of the
#' levels for a factor one, so a column like Status offered a box with no
#' indication of what could be typed into it. Columns with few distinct values
#' become factors and therefore dropdowns; high-cardinality ones (gene symbols,
#' site ids) stay text, where a search box is the right control and a list of
#' thousands of options is not.
#'
#' @param df data frame about to be handed to DT::datatable().
#' @param max_levels above this many distinct values a column stays text.
cpt_filter_levels <- function(df, max_levels = 60) {
  if (!is.data.frame(df) || !nrow(df)) return(df)
  for (nm in names(df)) {
    x <- df[[nm]]
    if (!is.character(x) && !is.logical(x)) next
    u <- unique(x[!is.na(x)])
    if (length(u) <= max_levels) df[[nm]] <- factor(x, levels = sort(as.character(u)))
  }
  df
}
