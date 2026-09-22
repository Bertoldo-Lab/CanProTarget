# ============================================================
# Script:   targets_panel.R
# Purpose:  Discover's layered target explorer.
#
# Design
#   The three evidence layers are independent toggles, so ligandability and
#   residue functionality can be examined without first calling dependencies.
#   Tiers are properties of a residue rather than of a dependency.
#
#   The row grain is explicit. Collapsing to one row per gene is what allows a
#   gene to appear "ligandable and functional" when the ligandable and
#   functional cysteines are different residues, which is the inference the
#   tiers exist to prevent. The grain control names what a row represents.
#
#   Everything reads data/gene_index.rds, so no layer costs a matrix load.
# ============================================================
# Sections:
#   targets_panel_ligandability_box()
#   targets_panel_residue_box()
#   targets_panel_actions()
#   targets_panel_results() — The single results region
#   cpt_empty_msg() — Competition ratio against probe promiscuity
#   isTRUE_vec()
#   cpt_add_hist()
#   cpt_cr_scatter()
#   targets_panel_server()
# ============================================================

CPT_LAYERS <- c(
  "Dependency"    = "dependency",
  "Ligandability" = "ligandability",
  "Cysteine function" = "functionality"
)

CPT_GRAINS <- c(
  "One row per gene"              = "gene",
  "One row per cysteine"          = "cysteine",
  "One row per cysteine x probe"  = "probe"
)

# The explorer's controls are two of the three layer boxes; the dependency box
# is built in dependencies_module so its inputs keep their existing ids. Each
# box carries its own include-checkbox, so a layer and its parameters are one
# object rather than a toggle here and its filters over there.
targets_panel_ligandability_box <- function(ns) {
  column(
    width = 4, class = "cpt-layer-box", id = ns("box_lig"),
    box(
      width = 12, status = "primary", solidHeader = TRUE, title = "Ligandability",
      actionButton(ns("tp_browse_lig"), "Browse the full ligandability results",
                   icon = icon("chart-column"),
                   class = "btn-default btn-block btn-sm cpt-browse-btn"),
      div(class = "cpt-box-summary", textOutput(ns("tp_lig_summary"), inline = TRUE)),
      tags$details(
        class = "cpt-box-details",
        tags$summary(class = "cpt-box-details-summary", "Change"),
        cpt_label("Minimum competition ratio", for_id = ns("tp_cr")),
        # The floor stops at 1.5, not 0. A competition ratio near 1 is equal
        # signal with and without the competitor, which is no engagement at
        # all, so the positions below it describe nothing: CR >= 1 keeps 5.5
        # million of the 10.6 million records and CR >= 0 keeps every one.
        # They also cost 20 s and more to re-tier, which read as the app
        # hanging. CR >= 1.5 keeps 792,000 and takes under three seconds.
        sliderInput(ns("tp_cr"), label = NULL, min = 1.5, max = 20, value = 4, step = 0.5),
        cpt_label("Maximum probe targets", for_id = ns("tp_max_targets")),
        numericInput(ns("tp_max_targets"), label = NULL, value = 20, min = 1,
                     max = 1252, step = 1)
      )
    )
  )
}

targets_panel_residue_box <- function(ns) {
  column(
    width = 4, class = "cpt-layer-box", id = ns("box_res"),
    box(
      width = 12, status = "primary", solidHeader = TRUE, title = "Cysteine function",
      actionButton(ns("tp_browse_atlas"), "Browse the full cysteine atlas",
                   icon = icon("bullseye"),
                   class = "btn-default btn-block btn-sm cpt-browse-btn"),
      div(class = "cpt-box-summary", textOutput(ns("tp_res_summary"), inline = TRUE)),
      tags$details(
        class = "cpt-box-details",
        tags$summary(class = "cpt-box-details-summary", "Change"),
        cpt_label("Evidence tiers", for_id = ns("tp_tiers")),
        checkboxGroupInput(
          ns("tp_tiers"), label = NULL,
          choices = c("Tier 1 · functional + ligandable" = "1",
                      "Tier 2 · functional only" = "2",
                      "Tier 3 · tested, not functional" = "3",
                      "Tier 4 · untested" = "4"),
          selected = c("1", "2", "3", "4")
        ),
        checkboxInput(ns("tp_clinvar"), "Require ClinVar pathogenic", FALSE),
        # Residue filters with no equivalent in the action bar. They apply to
        # the cysteine atlas browse view as well as the candidate table.
        #
        # Restricting to engaged sites is not among them: that is what the
        # Ligandability layer means, and every view except the raw atlas table
        # is already built from engaged sites when it is on. A checkbox for it
        # here was a second control for one idea, and the two could disagree.
        checkboxInput(ns("tp_atlas_lig"),
                      "Require atlas ligandability annotation (score > 50)", FALSE),
        cpt_label("Study context", for_id = ns("tp_context")),
        selectizeInput(ns("tp_context"), label = NULL, multiple = TRUE,
                       choices = c("PC14 and KMS26", "PC14-selective", "KMS26-selective"),
                       selected = NULL,
                       options = list(placeholder = "any context")),
        cpt_label("Editor support", for_id = ns("tp_editor")),
        selectInput(ns("tp_editor"), label = NULL,
                    choices = c("Any editor", "ABE", "CBE", "ABE and CBE"),
                    selected = "Any editor")
      )
    )
  )
}

# Grain and the action sit under the three boxes: they apply to whatever
# combination of layers is switched on above.
targets_panel_actions <- function(ns) {
  div(
    id = ns("tp_action_bar"),
    class = "cpt-action-bar cpt-action-grid",
    div(
      class = "cpt-action-grain",
      cpt_label("Evidence layers", for_id = ns("tp_layers")),
      # Nothing is on at first open. The tab then explains itself instead of
      # presenting a table the reader did not ask for and cannot interpret yet.
      checkboxGroupInput(ns("tp_layers"), label = NULL, choices = CPT_LAYERS,
                         selected = character(0),
                         inline = TRUE)
    ),
    div(
      class = "cpt-action-genes",
      cpt_label("Restrict to genes", for_id = ns("tp_genes")),
      # The placeholder is the only instruction most people will read, so it
      # has to show both forms rather than just the comma-separated one.
      textAreaInput(ns("tp_genes"), label = NULL, rows = 2,
                    placeholder = "optional \u2014 one per line or comma-separated, e.g. CDK2, BRD4, EZH2")
    ),
    div(
      class = "cpt-action-grain",
      cpt_label("Row grain", for_id = ns("tp_grain")),
      radioButtons(ns("tp_grain"), label = NULL, choices = CPT_GRAINS,
                   selected = "probe", inline = TRUE)
    ),
    # The CPT Score combines six evidence axes and only two of them are
    # dependency, so the weights belong with the controls that cut across the
    # layers rather than inside the Dependency box.
    div(
      class = "cpt-action-weights",
      tags$details(
        class = "cpt-collapse-block",
        tags$summary(class = "cpt-collapse-summary",
                     icon("sliders"), " CPT Score weights"),
        tags$div(
          class = "cpt-collapse-body",
          cpt_section("CPT Score weights", level = "h5"),
          fluidRow(
            column(4, sliderInput(ns("w_dependency_strength"), "Dependency strength",
              min = 0, max = 10, value = 3, step = 0.5)),
            column(4, sliderInput(ns("w_cancer_selectivity"), "Cancer selectivity",
              min = 0, max = 10, value = 3, step = 0.5)),
            column(4, sliderInput(ns("w_cysteine_ligandability"), "Cysteine ligandability",
              min = 0, max = 10, value = 2, step = 0.5))
          ),
          fluidRow(
            column(4, sliderInput(ns("w_conservation"), "Conservation",
              min = 0, max = 10, value = 1.5, step = 0.5)),
            column(4, sliderInput(ns("w_clinical_evidence"), "Clinical evidence (ClinVar)",
              min = 0, max = 10, value = 1, step = 0.5)),
            # ADME is computed (see cpt_build_adme_gene_scores) but defaults to
            # weight 0 so composite scores match the pre-implementation baseline
            # exactly. The slider is live: raising it changes scores and rankings.
            column(4, sliderInput(ns("w_adme_druggability"),
              "ADME / developability (provisional, default 0)",
              min = 0, max = 10, value = 0, step = 0.5))
          ),
          cpt_note(
            "Provisional formula, weight 0 by default.",
            "Scored for ~4,470 genes as the mean fragment developability (MW, cLogP, ",
            "TPSA, synthetic accessibility) of covalent probes engaging the gene at ",
            "CR >= 4. Reactivity filters are deliberately excluded: the Brenk alert on ",
            "971 of 1000 probes is the covalent warhead itself, so rewarding fewer ",
            "alerts would penalise probes for being covalent. Spread is narrow by ",
            "nature (IQR 92.6-98.7), so this axis shifts scores more than it reorders ",
            "them. Raising the weight changes every score and report, so treat it as a ",
            "methods decision. See ", tags$code("docs/CPT_SCORE.md"), " 4.6.",
            tags$br(), tags$br(),
            tags$strong("Cost: "),
            "computing this axis loads the 10.6M-row probe table (~20 s, ~0.9 GB), so it ",
            "runs only while the weight is above 0. At weight 0 the ADME column reads NA.",
            tone = "caution"
          ),
          fluidRow(
            column(3,
              actionButton(ns("cpt_weights_default"), "Defaults",
                class = "btn-default btn-sm btn-block",
                icon = icon("sliders"))
            ),
            column(3,
              actionButton(ns("cpt_weights_equal"), "Equal",
                class = "btn-default btn-sm btn-block",
                icon = icon("balance-scale"))
            )
          )
        )
      )
    ),
    # Run analysis lives here, with the layers, because it is the one action on
    # this page and it applies to whichever layers are on. Everything else
    # follows the controls without being asked (see the debounced `shown`
    # reactive), so this is the only button.
    div(
      class = "cpt-action-run",
      div(
        actionButton(ns("run_analysis"), "Run analysis",
                     icon = icon("play"), class = "btn-success"),
        actionButton(ns("reset_deps"), "Reset",
                     icon = icon("rotate-left"), class = "btn-default")
      ),
      # When the button wakes up, say which views are showing the previous
      # settings, so it is clear what clicking it will change.
      uiOutput(ns("run_stale_note"))
    )
  )
}


#' The single results region.
#'
#' Discover shows results in exactly one place. `analysis_tabs` are extra
#' tabPanels (the dependency and ligandability analysis output) that used to
#' live in their own box below this one; they are siblings of Candidates now,
#' so only one table and one set of plots is ever on screen.
#' Discover's results: one table space and one graph space.
#'
#' Every table the active layers can produce shares a single box, and every
#' chart shares another. The reader picks which one is on screen; nothing is
#' stacked and nothing is hidden behind a tab whose label they have to guess.
#' The two switchers are built on the server from the active layers, so the
#' choices are exactly the views that have something to say.
targets_panel_results <- function(ns, dep_extra = NULL,
                                  atlas_table = NULL, atlas_charts = NULL) {
  no_layers  <- paste0("!input['", ns("tp_layers"), "'] || ",
                       "input['", ns("tp_layers"), "'].length === 0")
  any_layers <- paste0("input['", ns("tp_layers"), "'] && ",
                       "input['", ns("tp_layers"), "'].length > 0")

  pick <- function(input_id, value, ...) {
    conditionalPanel(
      condition = paste0("input['", ns(input_id), "'] == '", value, "'"), ...)
  }
  # hide.ui = FALSE: an output that starts inside a hidden conditionalPanel
  # keeps shinycssloaders' hide class and lays out at zero height otherwise.
  spin <- function(x) shinycssloaders::withSpinner(x, type = 6,
                                                   color = "#478EB8", hide.ui = FALSE)
  tbl <- function(id) div(class = "cpt-scroll-x", spin(DTOutput(ns(id))))
  cht <- function(id) spin(plotlyOutput(ns(id), height = "440px"))
  dl  <- function(id, label = "Download shown rows (.xlsx)") {
    downloadButton(ns(id), label, class = "btn-sm btn-default")
  }

  tagList(
    conditionalPanel(
      condition = no_layers,
      div(
        class = "cpt-guide",
        cpt_section("How this tab works", level = "h4"),
        tags$p("Discover opens with nothing switched on. Choose the evidence a ",
               "candidate has to carry, and the table builds itself from there."),

        cpt_section("Three evidence layers", level = "h5"),
        tags$ul(
          tags$li(tags$b("Dependency"), " \u2014 genes a cancer subtype depends on ",
                  "more than the remaining models. Pick a subtype in the Dependency ",
                  "box, then run the analysis."),
          tags$li(tags$b("Ligandability"), " \u2014 cysteines a covalent probe engages, ",
                  "with the competition ratio measuring how completely."),
          tags$li(tags$b("Cysteine function"), " \u2014 what the base-editing atlas says ",
                  "about the residue itself: functional, tested and not functional, ",
                  "or never tested.")
        ),
        tags$p("The layers are independent. Switch one on and you are browsing ",
               "that evidence on its own; switch on two or three and a candidate ",
               "has to satisfy all of them. Each layer you enable adds its own ",
               "analysis tab beside this one, and its own columns to the table."),

        cpt_section("Then set two things", level = "h5"),
        tags$p(tags$b("Row grain"), " decides what one row means. A gene-level row ",
               "can read as ligandable and functional when those are two different ",
               "cysteines."),
        tags$p(tags$b("Restrict to genes"), " takes a pasted list \u2014 one per line, ",
               "commas, or a column straight out of a spreadsheet. Leave it empty to ",
               "search every gene.")
      )
    ),
    # One delegated handler for both switchers. Bound here rather than in each
    # rendered switcher, so re-rendering a row does not stack listeners.
    tags$script(HTML(
      "$(document).off('click.cptswitch').on('click.cptswitch', '.cpt-switch-btn', function(){",
      "  var b = $(this), w = b.closest('.cpt-switch');",
      "  w.find('.cpt-switch-btn').removeClass('active');",
      "  b.addClass('active');",
      "  Shiny.setInputValue(w.attr('data-input'), b.attr('data-value'));",
      "});"
    )),

    conditionalPanel(
      condition = any_layers,
      uiOutput(ns("tp_funnel")),
      # What the two spaces below are currently showing, always on screen.
      div(class = "cpt-scope-note", uiOutput(ns("tp_context_line"))),

      # ---- one table space ------------------------------------------------
      box(
        width = 12, status = "info", solidHeader = FALSE,
        uiOutput(ns("tp_table_switch")),
        uiOutput(ns("tp_table_note")),

        pick("tp_table_pick", "candidates",
             uiOutput(ns("tp_scope_note")), uiOutput(ns("tp_prompt")),
             tbl("tp_table"), uiOutput(ns("tp_drawer")),
             uiOutput(ns("tp_gene_warning")), br(), dl("tp_download")),
        pick("tp_table_pick", "dep_ranked",
             dep_extra),
        pick("tp_table_pick", "dep_all",
             tbl("tp_dep_table"), br(), dl("dl_dep_browse")),
        pick("tp_table_pick", "lig_best",
             tbl("best_probe_table"), br(), dl("dl_best_probe", "Download (.xlsx)")),
        pick("tp_table_pick", "lig_combined",
             tbl("combined_table"), br(), dl("dl_combined", "Download (.xlsx)")),
        pick("tp_table_pick", "lig_engaged",
             tbl("tp_lig_table"), br(), dl("dl_lig_browse")),
        pick("tp_table_pick", "res_engaged",
             tbl("res_table"), br(), dl("dl_res_analysis")),
        pick("tp_table_pick", "res_atlas", atlas_table),
        pick("tp_table_pick", "dep_res",
             tbl("tp_dep_res_table"), br(),
             dl("dl_dep_res", "Download (.xlsx)")),
        pick("tp_table_pick", "all_layers",
             tbl("tp_all_table"), br(),
             dl("dl_all_layers", "Download (.xlsx)")),
        pick("tp_table_pick", "res_integrated",
             uiOutput(ns("cys_probe_status")),
             tbl("cys_probe_table"), br(),
             dl("dl_cys_probe_targets", "Download integrated targets (.xlsx)"))
      ),

      # ---- one graph space ------------------------------------------------
      box(
        width = 12, status = "info", solidHeader = FALSE,
        uiOutput(ns("tp_chart_switch")),
        uiOutput(ns("tp_chart_note")),

        pick("tp_chart_pick", "dep_volcano",  cht("tp_dep_volcano")),
        pick("tp_chart_pick", "dep_hist",     cht("tp_dep_hist")),
        pick("tp_chart_pick", "dep_vs",       cht("tp_dep_vs")),
        pick("tp_chart_pick", "lig_scatter",  cht("tp_lig_scatter")),
        pick("tp_chart_pick", "lig_cr",       cht("tp_lig_cr_hist")),
        pick("tp_chart_pick", "lig_probes",   cht("tp_lig_probe_hist")),
        pick("tp_chart_pick", "lig_tiers",    cht("tp_lig_tiers")),
        pick("tp_chart_pick", "lig_select",   cht("tp_selectivity")),
        pick("tp_chart_pick", "dep_tier",     cht("tp_dep_tier")),
        pick("tp_chart_pick", "lig_vol_best", cht("lig_volcano_best")),
        pick("tp_chart_pick", "res_tiers",    cht("res_tier_chart")),
        pick("tp_chart_pick", "res_join",     cht("res_join_chart")),
        pick("tp_chart_pick", "res_abe_cbe",  atlas_charts$abe_cbe),
        pick("tp_chart_pick", "res_overview", atlas_charts$overview)
      )
    )
  )
}

#' Competition ratio against probe promiscuity.
#'
#' Two places draw this chart — the candidate set and the full ligandability
#' browse — so they share one builder. Same palette, same marker size, same
#' cutoff lines, same out-of-range handling, whichever one you are looking at.
#'
#' @param d Rows with CR, n_targets, evidence_tier, gene_name (cysteineid and
#'   probe_name are used in the hover when present)
#' @param cr,mt The two cutoffs, drawn as dashed lines
#' Add a histogram trace only when it has rows.
#'
#' plotly raises "Must supply `x` and/or `y` attributes" on a zero-length
#' trace, and an empty side of a cutoff is an ordinary state here, not an
#' error: the index stores engagement at CR >= 4, so at the default cutoff
#' the "below" side has no rows at all.
#' Empty-state text that names the actual cause.
#'
#' Saying "no engaged cysteines in the index" when the index is fine and the
#' gene restriction emptied the view sends the user to the wrong control.
#'
#' @param base Message to use when nothing is restricting the view
#' @param restricted TRUE when a gene restriction is active
cpt_empty_msg <- function(base, restricted) {
  if (isTRUE(restricted)) {
    "No rows match \u201cRestrict to genes\u201d. Clear or widen that box."
  } else {
    base
  }
}

#' TRUE-only test that treats NA as FALSE, for counting evidence flags.
isTRUE_vec <- function(x) !is.na(x) & x

cpt_add_hist <- function(p, x, name, colour, size) {
  if (!length(x)) return(p)
  plotly::add_histogram(p, x = x, name = name,
                        marker = list(color = colour),
                        xbins = list(size = size))
}

cpt_cr_scatter <- function(d, cr, mt, by_tier = TRUE) {
  # Evidence tier is residue evidence. With that layer off the chart is purely
  # ligandability, so it is drawn in one colour rather than encoding something
  # the reader did not ask to see.
  d <- d[!is.na(d$CR) & !is.na(d$n_targets), , drop = FALSE]
  if (!nrow(d)) return(NULL)

  # CR saturates at 20 and n_targets is a per-probe constant, so many records
  # land on identical coordinates. Jitter is display-only; hover and table
  # report the stored values.
  #
  # The jitter is seeded so the chart is stable between renders, but seeding
  # must not leak: set.seed() at render time would reset the global RNG for
  # everything else in the session. Save and restore it.
  old_seed <- if (exists(".Random.seed", .GlobalEnv)) {
    get(".Random.seed", .GlobalEnv)
  } else {
    NULL
  }
  on.exit({
    if (is.null(old_seed)) {
      suppressWarnings(rm(".Random.seed", envir = .GlobalEnv))
    } else {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(1)
  d$x <- d$n_targets * exp(stats::runif(nrow(d), -0.045, 0.045))
  d$tier_lab <- if (by_tier) paste("Tier", d$evidence_tier) else "Engaged cysteine"
  d$hover <- paste0(
    "<b>", d$gene_name, "</b>",
    if ("cysteineid" %in% names(d)) paste0(" ", d$cysteineid) else "",
    if ("probe_name" %in% names(d)) paste0("<br>probe ", d$probe_name) else "",
    "<br>CR ", signif(d$CR, 3), " · ", d$n_targets, " proteins engaged",
    if (by_tier) paste0("<br>Tier ", d$evidence_tier, ": ",
                        CPT_TIER_LABELS[as.character(d$evidence_tier)]) else "")

  # Records above the assay's bound of 20 stretch a linear axis to 2,500 and
  # flatten everything else onto the baseline. Drawn on a marked out-of-range
  # row with the true value in the hover, rather than dropped.
  oor <- d$CR > 20
  inr <- d[!oor, , drop = FALSE]
  inr$y <- inr$CR + stats::runif(nrow(inr), -0.18, 0.18)

  pal <- stats::setNames(unname(CPT_TIER_PAL[c("1", "2", "3", "4")]),
                         paste("Tier", 1:4))
  levels_drawn <- if (by_tier) paste("Tier", 4:1) else "Engaged cysteine"
  p <- plotly::plot_ly(source = "cpt_cr_scatter")
  for (lvl in levels_drawn) {
    dd <- inr[inr$tier_lab == lvl, , drop = FALSE]
    if (!nrow(dd)) next
    tier_n <- if (by_tier) substr(lvl, 6, 6) else NA_character_
    p <- plotly::add_trace(
      p, x = dd$x, y = dd$y, customdata = dd$gene_name,
      name = if (by_tier) paste0(lvl, " \u00b7 ", CPT_TIER_LABELS[[tier_n]]) else lvl,
      type = "scattergl", mode = "markers",
      marker = list(size = CPT_MARKER_SIZE,
                    opacity = if (by_tier) unname(CPT_TIER_ALPHA[[tier_n]]) else 0.85,
                    color = if (by_tier) unname(pal[[lvl]]) else CPT_PAL$primary,
                    line = list(width = 0)),
      text = dd$hover, hoverinfo = "text")
  }
  if (any(oor)) {
    dd <- d[oor, , drop = FALSE]
    p <- plotly::add_trace(
      p, x = dd$x, y = rep(20.9, nrow(dd)),
      name = "CR above 20 (out of range)",
      type = "scattergl", mode = "markers",
      marker = list(size = CPT_MARKER_SIZE_LG, symbol = "triangle-up",
                    color = CPT_PAL$alert, line = list(width = 0)),
      text = paste0(dd$hover, "<br><i>plotted on the out-of-range row</i>"),
      hoverinfo = "text")
  }
  cpt_plotly(
    p,
    xaxis = cpt_axis("Proteins engaged by the same probe (log scale)",
                     type = "log",
                     # plotly's autorange misbehaves once a "below" rect is
                     # added to a log axis, so the range is explicit.
                     range = c(log10(0.8), log10(1400)),
                     tickmode = "array",
                     tickvals = c(1, 5, 10, 20, 50, 100, 250, 500, 1000),
                     ticktext = c("1", "5", "10", "20", "50", "100", "250",
                                  "500", "1,000")),
    yaxis = cpt_axis("Competition ratio", range = c(3.3, 21.6)),
    shapes = list(
      # The region both cutoffs keep, shaded so the trade-off is readable
      # without reading the two lines off the axes.
      list(type = "rect", layer = "below",
           x0 = log10(0.8), x1 = log10(mt), y0 = cr, y1 = 21.6,
           fillcolor = "rgba(71, 142, 184, 0.09)", line = list(width = 0)),
      list(type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = cr, y1 = cr,
           line = list(color = CPT_PAL$alert, dash = "dash", width = 1.5)),
      list(type = "line", x0 = log10(mt), x1 = log10(mt), y0 = 0, y1 = 1,
           yref = "paper",
           line = list(color = CPT_PAL$alert, dash = "dash", width = 1.5))
    ),
    legend = list(orientation = "h", y = -0.2),
    margin = list(t = 20)
  )
}

targets_panel_server <- function(input, output, session, shared_data,
                                 get_subtype, get_dataset, extra_panel = NULL,
                                 on_open_gene = NULL) {

  # Results follow the controls. Debounced so dragging the CR slider or typing a
  # gene list does not recompute on every intermediate value; the button stays as
  # an explicit refresh for anyone who wants one.
  shown <- reactive({
    list(input$tp_layers, input$tp_grain, input$tp_cr, input$tp_max_targets,
         input$tp_tiers, input$tp_clinvar, input$tp_genes, input$tp_context,
         input$tp_editor, get_subtype(), get_dataset())
    TRUE
  }) |> debounce(600)

  # One results region, one view at a time. The browse buttons switch it rather
  # than appending panels below the candidates, so the page stays a fixed size
  # and discovery happens in a single place. Pressing the active button returns
  # to the candidate table.
  layers <- reactive(input$tp_layers %||% character(0))

  # A single line naming everything that decides what is on screen, so a reader
  # arriving at a table knows what produced it without hunting through boxes.
  output$tp_context_line <- renderUI({
    on <- layers()
    if (!length(on)) return(NULL)
    # CRISPR vs RNAi selects a dependency screen, so it says nothing about a
    # run with that layer off; naming it there implied the atlas was filtered
    # by a screen it has no connection to.
    bits <- character(0)
    st <- get_subtype()
    if ("dependency" %in% on) {
      bits <- c(bits, get_dataset())
      if (length(st) && nzchar(st[[1]])) bits <- c(bits, st[[1]])
    }
    lay <- c(dependency = "dependency", ligandability = "ligandability",
             functionality = "cysteine function")
    bits <- c(bits, paste(unname(lay[on]), collapse = " + "))
    # Each layer states the settings it was run under, so the line accounts for
    # everything on screen rather than only the chemistry.
    if ("dependency" %in% on) {
      bits <- c(bits, paste0("effect \u2264 ", input$effect_min %||% -0.1))
      if (isTRUE(input$apply_pvalue)) bits <- c(bits, "p<0.05")
      if (isTRUE(input$excl_common)) bits <- c(bits, "no common essentials")
      if (isTRUE(input$req_nc_sig)) bits <- c(bits, "significant vs non-cancer")
    }
    if ("ligandability" %in% on)
      bits <- c(bits, paste0("CR \u2265 ", input$tp_cr %||% 4,
                             " \u00b7 \u2264 ", input$tp_max_targets %||% 20, " targets"))
    if ("functionality" %in% on) {
      tiers <- input$tp_tiers %||% character(0)
      bits <- c(bits, if (length(tiers) == 4) "tiers 1-4"
                      else paste("tiers", paste(sort(tiers), collapse = ",")))
    }
    bits <- c(bits, switch(input$tp_grain %||% "probe",
                           gene = "one row per gene",
                           cysteine = "one row per cysteine",
                           probe = "one row per cysteine x probe"))
    g <- wanted_genes()
    if (!is.null(g)) bits <- c(bits, paste(length(g), "gene restriction"))
    tags$span(tags$b("Showing: "), paste(bits, collapse = " \u00b7 "))
  })
  outputOptions(output, "tp_context_line", suspendWhenHidden = FALSE)

  # The two switchers list exactly the views the active layers can produce, so
  # the reader never picks a table that has nothing behind it. A view that joins
  # two layers needs both of them on.
  # Reading order: the summary of every active layer first, then each layer's
  # own tables, then the joins that put two layers together. A reader who wants
  # the answer stops at the first; a reader checking the evidence works down.
  table_groups <- reactive({
    on <- layers()
    g <- list()
    g[["Overview"]] <- c("Summary" = "candidates")
    if ("dependency" %in% on)
      g[["Dependency"]] <- c("Cancer-selective genes" = "dep_ranked",
                             "Every gene in this subtype" = "dep_all")
    if ("ligandability" %in% on)
      g[["Ligandability"]] <- c("Best probe per protein" = "lig_best",
                                "Every engaged cysteine" = "lig_engaged")
    if ("functionality" %in% on)
      g[["Cysteine function"]] <- c("Engaged cysteines" = "res_engaged",
                                    "Cysteine atlas" = "res_atlas")
    combined <- character(0)
    if (all(c("dependency", "ligandability") %in% on))
      combined <- c(combined, "Dependency + ligandability" = "lig_combined")
    if (all(c("dependency", "functionality") %in% on))
      combined <- c(combined, "Dependency + cysteine function" = "dep_res")
    if (all(c("functionality", "ligandability") %in% on))
      combined <- c(combined, "Ligandability + cysteine function" = "res_integrated")
    # Last, the table that holds nothing back, once there is more than one
    # layer for it to combine.
    if (length(on) > 1L) combined <- c(combined, "All combined" = "all_layers")
    if (length(combined)) g[["Combined"]] <- combined
    g
  })

  table_choices <- reactive(unlist(unname(table_groups()), use.names = TRUE))

  chart_groups <- reactive({
    on <- layers()
    g <- list()
    if ("dependency" %in% on)
      g[["Dependency"]] <- c("Effect vs significance" = "dep_volcano",
                             "Effect-size distribution" = "dep_hist",
                             "Subtype vs other lines" = "dep_vs")
    if ("ligandability" %in% on)
      g[["Ligandability"]] <- c("Engagement vs selectivity" = "lig_scatter",
                                "Competition ratio" = "lig_cr",
                                "Probe promiscuity" = "lig_probes")
    if ("functionality" %in% on)
      g[["Cysteine function"]] <- c("Tier composition" = "res_tiers",
                                    "Evidence by join quality" = "res_join",
                                    "ABE versus CBE dropout" = "res_abe_cbe",
                                    "Evidence overview" = "res_overview")
    # Joins last, one for each pair of layers that is on.
    combined <- character(0)
    if (length(on) > 1L)
      combined <- c(combined, "Engagement vs selectivity (candidates)" = "lig_select")
    if (all(c("dependency", "ligandability") %in% on))
      combined <- c(combined, "Ligandable volcano" = "lig_vol_best")
    if (all(c("dependency", "functionality") %in% on))
      combined <- c(combined, "Dependency by tier" = "dep_tier")
    if (all(c("ligandability", "functionality") %in% on))
      combined <- c(combined, "Tiers of engaged sites" = "lig_tiers")
    if (length(combined)) g[["Combined"]] <- combined
    g
  })

  chart_choices <- reactive(unlist(unname(chart_groups()), use.names = TRUE))

  # One row per group, rather than one strip of every view. These are plain
  # buttons driving the input directly: the conditionalPanels downstream read
  # it client-side, so switching stays instant, which a server-rendered
  # selection would not.
  render_switch <- function(out_id, input_id, groups_r) {
    output[[out_id]] <- renderUI({
      gs <- groups_r()
      flat <- unlist(unname(gs), use.names = TRUE)
      if (!length(flat)) return(NULL)
      cur <- isolate(input[[input_id]])
      sel <- if (!is.null(cur) && cur %in% flat) cur else unname(flat[[1]])
      full_id <- session$ns(input_id)
      # Only label the rows when there is more than one to tell apart.
      show_tags <- length(gs) > 1L
      div(
        class = "cpt-switch", `data-input` = full_id,
        lapply(names(gs), function(g) {
          ch <- gs[[g]]
          div(
            class = "cpt-switch-row",
            if (show_tags) tags$span(class = "cpt-switch-tag", g),
            lapply(seq_along(ch), function(i) {
              val <- unname(ch[[i]])
              tags$button(
                type = "button",
                class = paste0("cpt-switch-btn",
                               if (identical(val, sel)) " active" else ""),
                `data-value` = val, names(ch)[i]
              )
            })
          )
        }),
        tags$script(HTML(sprintf(
          "Shiny.setInputValue('%s', '%s');", full_id, sel)))
      )
    })
    outputOptions(output, out_id, suspendWhenHidden = FALSE)
  }

  render_switch("tp_table_switch", "tp_table_pick", table_groups)
  render_switch("tp_chart_switch", "tp_chart_pick", chart_groups)

  # The layer sections are all on the page, so these buttons navigate rather
  # than switch views: each one turns its layer on if it is off, then scrolls to
  # that section. Jumping between layers stays one click, without anything being
  # hidden to make it work.
  jump_to <- function(input_id, layer, anchor_id) {
    observeEvent(input[[input_id]], {
      # Isolate the layer: "browse the full X results" means X on its own, not X
      # added to whatever else was already on.
      on <- isolate(layers())
      if (!identical(on, layer)) {
        updateCheckboxGroupInput(session, "tp_layers", selected = layer)
      }
      shinyjs::runjs(sprintf(
        "setTimeout(function(){var e=document.getElementById('%s'); if(e) e.scrollIntoView({behavior:'smooth', block:'start'});}, %d);",
        session$ns(anchor_id), if (identical(on, layer)) 100 else 900))
    }, ignoreInit = TRUE)
  }
  jump_to("tp_browse_dep",   "dependency",    "sec_dep")
  jump_to("tp_browse_lig",   "ligandability", "sec_lig")
  jump_to("tp_browse_atlas", "functionality", "sec_res")

  # Every section a chosen layer can show is on the page at once. There is no
  # view to switch and nothing folded behind a tab: switching a layer on adds
  # its tables and its charts below the candidate table, switching it off takes
  # them away. Hiding them was costing more in "where is that table?" than the
  # page length saves.
  output$tp_view_candidates <- renderText(if (length(layers())) "1" else "0")
  output$tp_browse_dep_on   <- renderText(if ("dependency" %in% layers()) "1" else "0")
  output$tp_browse_lig_on   <- renderText(if ("ligandability" %in% layers()) "1" else "0")
  output$tp_browse_atlas_on <- renderText(if ("functionality" %in% layers()) "1" else "0")
  for (o in c("tp_view_candidates", "tp_browse_dep_on", "tp_browse_lig_on",
              "tp_browse_atlas_on")) {
    outputOptions(output, o, suspendWhenHidden = FALSE)
  }

  # One-line summaries keep the applied filters visible while the boxes are
  # collapsed — without them, collapsing would hide the thresholds that decide
  # the numbers on screen.
  output$tp_lig_summary <- renderText({
    genes <- wanted_genes()
    paste0("CR \u2265 ", input$tp_cr %||% 4,
           " \u00b7 \u2264 ", input$tp_max_targets %||% 20, " targets",
           if (is.null(genes)) "" else paste0(" \u00b7 ", length(genes), " gene",
                                              if (length(genes) == 1) "" else "s"))
  })

  output$tp_res_summary <- renderText({
    tiers <- input$tp_tiers %||% character(0)
    tier_txt <- if (!length(tiers)) "no tiers" else if (length(tiers) == 4) {
      "tiers 1-4"
    } else {
      paste0("tier", if (length(tiers) > 1) "s" else "", " ", paste(sort(tiers), collapse = ", "))
    }
    ctx <- input$tp_context
    paste0(tier_txt,
           " \u00b7 ", if (!length(ctx)) "any context" else paste(length(ctx), "contexts"),
           " \u00b7 ", tolower(input$tp_editor %||% "any editor"),
           if (isTRUE(input$tp_clinvar)) " \u00b7 ClinVar only" else "")
  })

  observe({
    on <- layers()
    any_on <- length(on) > 0
    shinyjs::toggleClass(id = session$ns("box_lig"), class = "cpt-box-off",
                         condition = any_on && !("ligandability" %in% on), asis = TRUE)
    shinyjs::toggleClass(id = session$ns("box_res"), class = "cpt-box-off",
                         condition = any_on && !("functionality" %in% on), asis = TRUE)
    shinyjs::toggleClass(id = session$ns("box_dep"), class = "cpt-box-off",
                         condition = any_on && !("dependency" %in% on), asis = TRUE)
  })

  # Pasted gene list: batch lookup across every layer.
  wanted_genes <- reactive({
    raw <- input$tp_genes
    if (is.null(raw) || !nzchar(trimws(raw))) return(NULL)
    g <- unlist(strsplit(raw, "[,;\n\t ]+"))
    g <- cpt_gene_match_key(trimws(g))
    g <- g[nzchar(g)]
    if (!length(g)) NULL else unique(g)
  })

  restrict_to_genes <- function(df) {
    want <- wanted_genes()
    if (is.null(want) || is.null(df) || !nrow(df)) return(df)
    df[df$gene_key %in% want, , drop = FALSE]
  }
  has_layer <- function(x) x %in% layers()

  # The engaged-site set for the current competition-ratio floor.
  #
  # The index is built at CR >= 4, the engagement threshold the manuscript
  # uses, so it cannot answer for a lower floor: asking for CR >= 2 against it
  # returned the CR >= 4 answer unchanged, while the binding table holds 25,593
  # distinct cysteines at that floor against 9,753. Below 4 the set is rebuilt
  # from the binding table by the same recipe the index builder uses
  # (docs/scripts/build_gene_index.R), so the two agree at any threshold.
  engaged_sites_for_cutoff <- reactive({
    idx <- shared_data$gene_index()
    cr <- suppressWarnings(as.numeric(input$tp_cr %||% 4))
    if (!length(cr) || !is.finite(cr)) cr <- 4
    if (cr >= 4) return(idx$engaged_sites)

    pb <- tryCatch({
      if (!is.null(shared_data$protein_binding_lookup)) shared_data$protein_binding_lookup()
      else shared_data$proteinbindinglookup()
    }, error = function(e) NULL)
    if (is.null(pb) || !nrow(pb)) return(idx$engaged_sites)

    keep <- !is.na(pb$CR) & pb$CR >= cr
    # Insurance against a floor the slider no longer offers, reached by a
    # bookmarked URL or a restored session: re-tiering millions of rows takes
    # tens of seconds and looks like a hang.
    shiny::validate(need(
      sum(keep) <= 1.2e6,
      paste0("A competition ratio of ", cr, " keeps ",
             format(sum(keep), big.mark = ","),
             " records, too many to tier here. Raise the floor to 1.5 or above.")))
    eng <- pb[keep, intersect(c("gene_name", "cysteineid", "probe_name", "CR", "n_targets"),
                              names(pb)), drop = FALSE]
    eng <- as.data.frame(eng, stringsAsFactors = FALSE)
    for (nm in c("gene_name", "cysteineid", "probe_name")) {
      if (nm %in% names(eng)) eng[[nm]] <- as.character(eng[[nm]])
    }
    atlas <- tryCatch(shared_data$cys_editing_atlas(), error = function(e) NULL)
    if (is.null(atlas)) return(idx$engaged_sites)
    out <- tryCatch(cpt_annotate_engaged_tiers(eng, atlas), error = function(e) NULL)
    if (is.null(out) || !nrow(out)) return(idx$engaged_sites)
    if ("cys_site_id" %in% names(out)) {
      out$cys_clinvar_pathogenic <- cpt_is_true(
        atlas$clinvar_pathogenic[match(out$cys_site_id, atlas$site_id)])
    }
    out
  })

  # Site-grain rows, before the dependency layer is applied.
  site_rows <- reactive({
    idx <- shared_data$gene_index()
    shiny::validate(need(!is.null(idx), "Gene index not available."))

    if (has_layer("ligandability")) {
      df <- engaged_sites_for_cutoff()
      shiny::validate(need(!is.null(df) && nrow(df) > 0, "No engaged sites in the index."))
      df <- df[!is.na(df$CR) & df$CR >= (input$tp_cr %||% 4), , drop = FALSE]
      if (!is.null(input$tp_max_targets) && "n_targets" %in% names(df)) {
        df <- df[is.na(df$n_targets) | df$n_targets <= input$tp_max_targets, , drop = FALSE]
      }
      df$gene_key <- cpt_gene_match_key(df$gene_name)
      return(restrict_to_genes(df))
    }

    # Functionality without ligandability: every atlas-tested site, engaged or
    # not, so residue evidence can be browsed on its own terms.
    atlas <- shared_data$cys_editing_atlas()
    shiny::validate(need(!is.null(atlas), "Cysteine atlas not available."))
    restrict_to_genes(data.frame(
      gene_key = cpt_gene_match_key(atlas$gene_symbol),
      gene_name = atlas$gene_symbol,
      cysteineid = paste0(atlas$uniprot_accession, "_C", atlas$cysteine_position),
      probe_name = NA_character_, CR = NA_real_, n_targets = NA_integer_,
      evidence_tier = cpt_evidence_tier(TRUE, cpt_is_true(atlas$functional),
                                        cpt_is_true(atlas$ligandable)),
      cys_functional = cpt_is_true(atlas$functional),
      cys_ligandable = cpt_is_true(atlas$ligandable),
      cys_clinvar_pathogenic = cpt_is_true(atlas$clinvar_pathogenic),
      cys_join_status = "UniProt + cysteine position",
      cys_study_context = atlas$study_context,
      cys_editor_support = atlas$editor_support,
      stringsAsFactors = FALSE
    ))
  })

  # The funnel: how many genes survive each layer, in the manuscript's order.
  funnel <- reactive({
    idx <- shared_data$gene_index()
    req(idx)
    steps <- list()
    rows <- NULL

    dep_genes <- NULL
    if (has_layer("dependency")) {
      st <- get_subtype()
      ds <- get_dataset()
      dep <- idx$dependency
      dep <- dep[dep$dataset == ds & dep$subtype %in% st, , drop = FALSE]
      dep_genes <- unique(dep$gene_key)
      steps[["Dependency"]] <- length(dep_genes)
    }

    if (has_layer("ligandability") || has_layer("functionality")) {
      rows <- site_rows()
      if (!is.null(dep_genes)) rows <- rows[rows$gene_key %in% dep_genes, , drop = FALSE]
      if (has_layer("ligandability")) {
        steps[["Ligandability"]] <- length(unique(rows$gene_key))
      }
      if (has_layer("functionality")) {
        keep <- as.character(rows$evidence_tier) %in% (input$tp_tiers %||% character(0))
        if (isTRUE(input$tp_clinvar)) keep <- keep & cpt_is_true(rows$cys_clinvar_pathogenic)
        # Atlas context and editor filters, carried over from the atlas tab so
        # residue evidence is filtered in one place rather than two.
        ctx <- input$tp_context
        if (length(ctx) && "cys_study_context" %in% names(rows)) {
          keep <- keep & (is.na(rows$cys_study_context) |
                            rows$cys_study_context %in% ctx)
        }
        ed <- input$tp_editor %||% "Any editor"
        if (!identical(ed, "Any editor") && "cys_editor_support" %in% names(rows)) {
          keep <- keep & (!is.na(rows$cys_editor_support) &
                            rows$cys_editor_support == ed)
        }
        rows <- rows[keep, , drop = FALSE]
        steps[["Cysteine function"]] <- length(unique(rows$gene_key))
      }
    } else if (!is.null(dep_genes)) {
      # Dependency on its own: the rows still need a display symbol, or the
      # table renders with no columns at all.
      gi <- idx$genes
      rows <- data.frame(
        gene_key = dep_genes,
        gene_name = gi$symbol[match(dep_genes, gi$gene_key)],
        stringsAsFactors = FALSE
      )
      rows$gene_name[is.na(rows$gene_name)] <- rows$gene_key[is.na(rows$gene_name)]
    }

    list(steps = steps, rows = rows, dep_genes = dep_genes)
  })

  candidates <- reactive({
    f <- funnel()
    rows <- f$rows
    # An empty result is a legitimate answer (a filter that excludes everything),
    # so it must reach the table as zero rows. Returning a validation error here
    # left DT showing the previous query's rows.
    shiny::validate(need(length(layers()) > 0,
                         "Enable at least one evidence layer to see candidates."))
    if (is.null(rows) || !nrow(rows)) return(rows)
    idx <- shared_data$gene_index()
    grain <- input$tp_grain %||% "probe"

    if (!all(c("cysteineid") %in% names(rows))) grain <- "gene"

    # Order by tier then engagement strength where those columns exist. With the
    # dependency layer alone they do not, and ordering on a missing column threw
    # "argument 1 is not a vector" — which DT then papered over with stale rows.
    order_rows <- function(d) {
      keys <- list()
      if ("evidence_tier" %in% names(d)) keys <- c(keys, list(d$evidence_tier))
      if ("CR" %in% names(d)) keys <- c(keys, list(-ifelse(is.na(d$CR), -Inf, d$CR)))
      if (!length(keys)) return(d)
      d[do.call(order, keys), , drop = FALSE]
    }

    out <- switch(
      grain,
      probe = order_rows(rows),
      cysteine = {
        d <- order_rows(rows)
        if ("cysteineid" %in% names(d)) d[!duplicated(d$cysteineid), , drop = FALSE] else d
      },
      gene = {
        d <- order_rows(rows)
        d[!duplicated(d$gene_key), , drop = FALSE]
      }
    )
    if (!nrow(out)) return(out)

    # How many distinct probes engage this row, counted before the dedup above
    # collapses a site to its best probe. At gene grain that is every probe
    # reaching the gene; otherwise every probe reaching the cysteine.
    if ("probe_name" %in% names(rows)) {
      key_col <- if (grain == "gene") "gene_key" else "cysteineid"
      if (key_col %in% names(rows) && key_col %in% names(out)) {
        n_by <- tapply(as.character(rows$probe_name),
                       as.character(rows[[key_col]]),
                       function(x) length(unique(x)))
        out$n_probes <- as.integer(n_by[as.character(out[[key_col]])])
      }
    }

    # Attach dependency columns when that layer is on.
    if (has_layer("dependency")) {
      dep <- idx$dependency
      dep <- dep[dep$dataset == get_dataset() & dep$subtype %in% get_subtype(), , drop = FALSE]
      m <- match(out$gene_key, dep$gene_key)
      out$effect_size <- dep$effect_size[m]
      out$p_value <- dep$p_value[m]
      out$subtype <- dep$subtype[m]
    }
    out
  })

  # The candidate rows with everything joined on: the identity columns, and
  # the corrected significance and cohort means that live in the per-subtype
  # effect-size file rather than in the index. Both the summary and the
  # all-layers table read this, so they cannot disagree about a value.
  enriched_rows <- reactive({
    out <- candidates()
    if (is.null(out)) out <- data.frame()
    if (!nrow(out)) return(out)

    # Identity first, in the vocabulary the rest of the app uses: the gene is
    # the dependency unit, the protein (UniProt accession) is what cysteine
    # numbering belongs to, and the site is a position within that protein.
    if ("cysteineid" %in% names(out)) {
      cid <- as.character(out$cysteineid)
      out$protein_id <- sub("_C.*$", "", cid)
      out$site <- sub("^.*_C", "C", cid)
    }
    if (has_layer("dependency") && "gene_name" %in% names(out)) {
      uni <- tryCatch(dep_universe(), error = function(e) NULL)
      if (!is.null(uni) && nrow(uni)) {
        m <- match(cpt_gene_match_key(out$gene_name), cpt_gene_match_key(uni$gene_name))
        take <- function(target, candidates, digits) {
          src <- intersect(candidates, names(uni))[1]
          if (is.na(src)) return(invisible(NULL))
          out[[target]] <<- round(suppressWarnings(as.numeric(uni[[src]][m])), digits)
        }
        take("q_value",           c("q.value", "q_value"), 4)
        take("neg_log10_p",       c("neg_log10_p_value", "neg_log10_p"), 2)
        take("cancer_avg",        c("Cancer_Avg", "cancer_avg"), 3)
        take("noncancer_avg",     c("NonCancer_Avg", "noncancer_avg"), 3)
        take("avg",               c("Avg", "avg"), 3)
        take("pval_vs_noncancer", c("pval_vs_NonCancer", "pval_vs_noncancer"), 4)
      }
    }
    out
  })

  # Everything the active layers know about each row, unabridged. The summary
  # above it is a selection -- the columns a reader judges a candidate on --
  # so the last table in the list is the one that holds nothing back: the
  # assay context, the structure, and the atlas detail the summary leaves out.
  all_layers_df <- reactive({
    d <- enriched_rows()
    if (!nrow(d)) return(data.frame(`No rows` = character(0), check.names = FALSE))

    # SMILES, dataset and cell line are measurements, so they come from the
    # binding table rather than the index.
    cr4 <- tryCatch(shared_data$protein_binding_cr4(), error = function(e) NULL)
    if (is.data.frame(cr4) && nrow(cr4) && "probe_name" %in% names(d)) {
      k <- paste(tolower(as.character(cr4$probe_name)), as.character(cr4$cysteineid))
      m <- match(paste(tolower(as.character(d$probe_name)), as.character(d$cysteineid)), k)
      for (nm in c("SMILES", "Dataset", "Cell_Line")) {
        if (nm %in% names(cr4)) d[[nm]] <- as.character(cr4[[nm]])[m]
      }
    }

    col <- function(nm) if (nm %in% names(d)) d[[nm]] else NULL
    yn  <- function(nm) { x <- col(nm); if (is.null(x)) NULL else ifelse(cpt_is_true(x), "yes", "no") }
    out <- data.frame(
      Gene = col("gene_name"), Protein = col("protein_id"), Site = col("site"),
      check.names = FALSE, stringsAsFactors = FALSE
    )
    add <- function(label, value) if (!is.null(value)) out[[label]] <<- value

    num <- function(nm, digits) {
      x <- col(nm)
      if (is.null(x)) NULL else round(suppressWarnings(as.numeric(x)), digits)
    }
    if (has_layer("dependency")) {
      add("Subtype", col("subtype"))
      # The index stores these at full precision; printing sixteen decimals of
      # a gene effect says nothing the fourth does not.
      add("Effect size", num("effect_size", 4))
      add("p", num("p_value", 4))
      add("q", col("q_value"))
      add("-log10 p", col("neg_log10_p"))
      add("Mean in subtype", col("cancer_avg"))
      add("Mean in non-cancer", col("noncancer_avg"))
      add("Mean across all lines", col("avg"))
      add("p vs non-cancer", col("pval_vs_noncancer"))
    }
    if (has_layer("ligandability")) {
      add("Probe", col("probe_name"))
      add("SMILES", col("SMILES"))
      add("CR", num("CR", 2))
      add("Probe targets", col("n_targets"))
      add("Probes available", col("n_probes"))
      add("Dataset", col("Dataset"))
      add("Cell line", col("Cell_Line"))
    }
    if (has_layer("functionality")) {
      add("Tier", col("evidence_tier"))
      add("Tier meaning", col("evidence_tier_label"))
      add("Functional (atlas)", yn("cys_functional"))
      add("Atlas ligandability", yn("cys_ligandable"))
      add("ClinVar", yn("cys_clinvar_pathogenic"))
      add("Editor support", col("cys_editor_support"))
      add("Study context", col("cys_study_context"))
      add("Residue mapping", col("cys_residue_mapping_status"))
      add("Atlas site id", col("cys_site_id"))
      add("Matched on", col("cys_join_status"))
    }

    # The composite belongs in the table that holds everything. It is a score
    # across layers, so it stays out of the per-layer views, where it would
    # imply evidence those views do not show.
    if (has_layer("dependency") && !is.null(shared_data$cpt_scores)) {
      cg <- tryCatch(shared_data$cpt_scores(), error = function(e) NULL)
      if (is.data.frame(cg) && nrow(cg) && "gene_name" %in% names(cg)) {
        m <- match(cpt_gene_match_key(d$gene_name), cpt_gene_match_key(cg$gene_name))
        if ("CPT_Score" %in% names(cg)) add("CPT Score", round(cg$CPT_Score[m], 2))
        if ("CPT_Rank" %in% names(cg))  add("CPT rank", cg$CPT_Rank[m])
      }
    }
    if (!nrow(out)) return(out)
    out[, !vapply(out, function(x) all(is.na(x)), logical(1)), drop = FALSE]
  })

  # Dependency + cysteine function. The pair that has no chemistry in it: how
  # strongly the subtype depends on the gene, beside what the atlas knows about
  # the residue. Probe columns are deliberately absent -- that is the other two
  # joins' business.
  dep_res_df <- reactive({
    d <- enriched_rows()
    if (!nrow(d)) return(data.frame(`No rows` = character(0), check.names = FALSE))
    col <- function(nm) if (nm %in% names(d)) d[[nm]] else NULL
    yn  <- function(nm) { x <- col(nm); if (is.null(x)) NULL else ifelse(cpt_is_true(x), "yes", "no") }
    out <- data.frame(
      Gene = col("gene_name"), Protein = col("protein_id"), Site = col("site"),
      check.names = FALSE, stringsAsFactors = FALSE
    )
    add <- function(label, value) if (!is.null(value)) out[[label]] <<- value
    add("Subtype", col("subtype"))
    add("Effect size", col("effect_size"))
    add("p", col("p_value"))
    add("q", col("q_value"))
    add("-log10 p", col("neg_log10_p"))
    add("Mean in subtype", col("cancer_avg"))
    add("Mean in non-cancer", col("noncancer_avg"))
    add("Tier", col("evidence_tier"))
    add("Tier meaning", col("evidence_tier_label"))
    add("Functional (atlas)", yn("cys_functional"))
    add("Atlas ligandability", yn("cys_ligandable"))
    add("ClinVar", yn("cys_clinvar_pathogenic"))
    add("Editor support", col("cys_editor_support"))
    add("Study context", col("cys_study_context"))
    add("Matched on", col("cys_join_status"))
    if (!nrow(out)) return(out)
    out <- out[, !vapply(out, function(x) all(is.na(x)), logical(1)), drop = FALSE]
    if ("Tier" %in% names(out)) out <- out[order(out$Tier, na.last = TRUE), , drop = FALSE]
    out
  })

  output$tp_dep_res_table <- DT::renderDT({
    DT::datatable(dep_res_df(), rownames = FALSE, filter = "top", selection = "none",
                  options = list(pageLength = 15, scrollX = TRUE))
  })

  output$dl_dep_res <- downloadHandler(
    filename = function() paste0("canprotarget_dependency_cysteine_", Sys.Date(), ".xlsx"),
    content = function(file) {
      df <- dep_res_df()
      rows <- input$tp_dep_res_table_rows_all
      if (!is.null(rows) && length(rows)) df <- df[rows, , drop = FALSE]
      write_xlsx(df, file)
    }
  )

  output$tp_all_table <- DT::renderDT({
    DT::datatable(all_layers_df(), rownames = FALSE, filter = "top",
                  selection = "none",
                  options = list(pageLength = 15, scrollX = TRUE))
  })

  output$dl_all_layers <- downloadHandler(
    filename = function() paste0("canprotarget_all_layers_", Sys.Date(), ".xlsx"),
    content = function(file) {
      df <- all_layers_df()
      rows <- input$tp_all_table_rows_all
      if (!is.null(rows) && length(rows)) df <- df[rows, , drop = FALSE]
      write_xlsx(df, file)
    }
  )

  display_df <- reactive({
    out <- enriched_rows()
    grain <- input$tp_grain %||% "probe"

    cols <- c(Gene = "gene_name", Protein = "protein_id")
    if (grain != "gene") cols <- c(cols, Site = "site")

    # Then one block per enabled layer, in the order of the funnel above.
    if (has_layer("dependency")) {
      # The dependency block carries the numbers a reader judges an association
      # on: the effect, its nominal and corrected significance, and the cohort
      # means it is a contrast against.
      cols <- c(cols, Subtype = "subtype", `Effect size` = "effect_size",
                p = "p_value", q = "q_value", `-log10 p` = "neg_log10_p",
                `Mean in subtype` = "cancer_avg",
                `Mean in non-cancer` = "noncancer_avg",
                `Mean across all lines` = "avg",
                `p vs non-cancer` = "pval_vs_noncancer")
    }
    if (has_layer("functionality")) {
      cols <- c(cols, Tier = "evidence_tier",
                `Functional (atlas)` = "cys_functional",
                `Atlas ligandability` = "cys_ligandable",
                ClinVar = "cys_clinvar_pathogenic",
                `Matched on` = "cys_join_status")
    }
    if (has_layer("ligandability")) {
      # Name the probe behind the competition ratio. At cysteine and gene
      # grain the retained row is the highest-CR probe for that site (see
      # order_rows + !duplicated), so it is reported as the best probe.
      cols <- c(cols,
                if (grain == "probe") c(Probe = "probe_name")
                else c(`Best probe` = "probe_name"))
      cols <- c(cols, CR = "CR", `Probe targets` = "n_targets",
                `Probes available` = "n_probes")
    }
    cols <- cols[cols %in% names(out)]
    if (!length(cols) || !nrow(out)) {
      empty <- as.data.frame(stats::setNames(
        rep(list(character(0)), max(length(cols), 1)),
        if (length(cols)) names(cols) else "No rows"
      ), check.names = FALSE)
      return(empty)
    }
    df <- out[, unname(cols), drop = FALSE]
    names(df) <- names(cols)
    num <- vapply(df, is.numeric, logical(1))
    df[num] <- lapply(df[num], function(x) round(x, 4))
    # "true"/"false" is how the data frame prints, not how a result should
    # read. An em dash for FALSE keeps the column scannable; NA stays blank,
    # because unknown and negative are different answers.
    lgl <- vapply(df, is.logical, logical(1))
    df[lgl] <- lapply(df[lgl], function(x) {
      ifelse(is.na(x), "", ifelse(x, "Yes", "\u2014"))
    })
    df
  })

  # Each step narrows the one before it, so only the first counts genes
  # outright and the rest say what they added. "genes after functionality" was
  # both awkward and a name the interface does not use for that layer.
  funnel_label <- function(step, n, i) {
    unit <- if (n == 1L) "gene" else "genes"
    if (i == 1L) {
      switch(step,
        "Dependency"      = paste(if (n == 1L) "dependency" else "dependencies",
                                  "in this subtype"),
        "Ligandability"   = paste(unit, "engaged by a covalent probe"),
        "Cysteine function" = paste(unit, "with a tested cysteine"),
        paste(unit, "after", tolower(step)))
    } else {
      switch(step,
        "Ligandability"   = "also engaged by a covalent probe",
        "Cysteine function" = "also with a tested cysteine",
        "Dependency"      = "also a dependency here",
        paste("also after", tolower(step)))
    }
  }

  output$tp_funnel <- renderUI({
    if (!shown() || !length(layers())) return(NULL)
    f <- tryCatch(funnel(), error = function(e) NULL)
    if (is.null(f) || !length(f$steps)) return(NULL)
    steps <- f$steps
    tags$div(
      class = "cpt-funnel",
      lapply(seq_along(steps), function(i) {
        tags$div(
          class = "cpt-funnel-step",
          tags$span(class = "cpt-funnel-n", format(steps[[i]], big.mark = ",")),
          tags$span(class = "cpt-funnel-label", funnel_label(names(steps)[i],
                                                              steps[[i]], i))
        )
      })
    )
  })

  # Say what the table is a table OF. Without the dependency layer these are
  # every engaged residue in the data, not results for a particular cancer, and
  # a reader arriving mid-page cannot tell the difference.
  output$tp_scope_note <- renderUI({
    if (!length(layers())) return(NULL)
    if ("dependency" %in% layers()) {
      st <- get_subtype()
      st <- if (length(st) && nzchar(st[1])) st[1] else "no subtype selected"
      return(tags$p(class = "cpt-scope-note",
                    "Candidates in ", tags$b(st), ", filtered by every layer below."))
    }
    # Nothing to say when no cancer is in scope: the Dependency checkbox is
    # directly above and shows its own state.
    NULL
  })

  output$tp_prompt <- renderUI(NULL)

  output$tp_table <- renderDT({
    if (!shown()) return(NULL)
    # With every layer off there is nothing to ask for. A validation error here
    # leaves the previous query's rows on screen, so return an empty table.
    if (!length(layers())) {
      return(datatable(
        data.frame(`No evidence layer selected` = character(0), check.names = FALSE),
        rownames = FALSE,
        options = list(dom = "t", language = list(
          emptyTable = "Switch on at least one evidence layer to see candidates."))
      ))
    }
    df <- display_df()
    # Row click opens the gene drawer; single selection keeps that unambiguous.
    # emptyTable rather than validate(): a validation error inside renderDT
    # leaves the previous query's rows on screen.
    datatable(df, rownames = FALSE, filter = "top", selection = "single",
              options = list(
                pageLength = 15, scrollX = TRUE,
                language = list(emptyTable = paste(
                  "No candidate survives the current layers and filters.",
                  "Loosen a threshold, widen the tier selection, or clear the",
                  "gene list."))
              ))
  })

  # Clicking a row opens the gene beside the table rather than sending the user
  # to another tab. Same card component the Gene tab's data feeds, so gene
  # information keeps one source.
  # Driven by cell clicks rather than row selection: any click in the row opens
  # the gene, with no separate "select the row first" step.
  drawer_row <- reactiveVal(NULL)

  observeEvent(input$tp_table_cell_clicked, {
    info <- input$tp_table_cell_clicked
    if (is.null(info) || is.null(info$row)) return()
    drawer_row(as.integer(info$row))
  }, ignoreInit = TRUE)

  # A new query closes the open drawer row.
  observeEvent(list(input$tp_layers, input$tp_grain, input$tp_tiers, input$tp_cr,
                    input$tp_max_targets, input$tp_clinvar),
               drawer_row(NULL), ignoreInit = TRUE)

  output$tp_drawer <- renderUI({
    if (!shown()) return(NULL)
    sel <- drawer_row()
    if (is.null(sel) || !length(sel) || is.na(sel)) return(NULL)
    df <- tryCatch(candidates(), error = function(e) NULL)
    if (is.null(df) || nrow(df) < sel[[1]]) return(NULL)
    key <- df$gene_key[sel[[1]]]
    idx <- shared_data$gene_index()
    tags$div(
      class = "cpt-drawer",
      tags$div(
        class = "cpt-drawer-head",
        actionLink(session$ns("tp_drawer_gene"), "Open in Target tab",
                   icon = icon("arrow-up-right-from-square"),
                   class = "cpt-drawer-open"),
        actionLink(session$ns("tp_drawer_close"), NULL, icon = icon("xmark"),
                   class = "cpt-drawer-close", title = "Close")
      ),
      tags$div(class = "cpt-drawer-body",
               cpt_gene_card(idx, key, dataset = get_dataset()),
               if (is.function(extra_panel)) extra_panel(key))
    )
  })

  observeEvent(input$tp_drawer_close, {
    drawer_row(NULL)
    DT::selectRows(DT::dataTableProxy("tp_table"), NULL)
  }, ignoreInit = TRUE)

  output$tp_gene_warning <- renderUI({
    if (!shown()) return(NULL)
    idx <- tryCatch(shared_data$gene_index(), error = function(e) NULL)
    if (is.null(idx)) return(NULL)
    want <- wanted_genes()

    # When nothing survives, name the filter responsible rather than leaving the
    # user to guess: the <= N target rule excludes promiscuous probes, and that
    # is usually what empties a gene list (CDK2 and EZH2 are engaged at CR 10
    # and 20, but only by probes hitting far more than 20 proteins).
    rows_now <- tryCatch(candidates(), error = function(e) NULL)
    blocked <- NULL
    if ((is.null(rows_now) || !nrow(rows_now)) && has_layer("ligandability")) {
      es <- engaged_sites_for_cutoff()
      if (!is.null(es) && nrow(es)) {
        keys <- cpt_gene_match_key(es$gene_name)
        pool <- if (is.null(want)) es else es[keys %in% want, , drop = FALSE]
        if (nrow(pool)) {
          blocked <- cpt_note(
            sprintf("%d engaged cysteine%s excluded by the probe-target limit.",
                    nrow(pool), if (nrow(pool) == 1) "" else "s"),
            "These residues are engaged at CR >= ", input$tp_cr %||% 4,
            ", but only by probes that also engage more than ",
            input$tp_max_targets %||% 20, " other proteins. Raise the maximum ",
            "probe targets to see them; the limit is a selectivity filter, not a ",
            "statement that the site is unengaged.",
            tone = "caution"
          )
        }
      }
    }

    missing <- if (is.null(want)) character(0) else setdiff(want, idx$genes$gene_key)
    if (!length(missing)) return(blocked)
    tagList(blocked, cpt_note(
      sprintf("%d of %d symbols are not in the index.", length(missing), length(want)),
      "Not found: ", paste(utils::head(missing, 25), collapse = ", "),
      if (length(missing) > 25) " …" else NULL,
      ". They may be aliases, or absent from the dependency, chemoproteomics and ",
      "base-editing sources entirely.",
      tone = "caution"
    ))
  })

  output$tp_selectivity <- renderPlotly({
    shiny::validate(need(shown(), ""))
    df <- candidates()
    shiny::validate(need(has_layer("ligandability"),
                         "Turn on the ligandability layer to see probe selectivity."))
    shiny::validate(need(
      "CR" %in% names(df) && any(!is.na(df$CR)),
      cpt_empty_msg("No probe records in the current selection.",
                    !is.null(wanted_genes()))))
    p <- cpt_cr_scatter(df, input$tp_cr %||% 4, input$tp_max_targets %||% 20,
                        by_tier = "functionality" %in% layers())
    shiny::validate(need(!is.null(p), "No probe records in the current selection."))
    p
  })

  # ---- Browse: full distributions with the cutoffs drawn on them ------------

  # Whole-subtype effect sizes, read straight from the precomputed cache so the
  # excluded side of the cutoff is visible too (the index keeps only passers).
  dep_universe <- reactive({
    st <- get_subtype()
    ds <- get_dataset()
    req(length(st) > 0, nzchar(st[[1]]))
    path <- cpt_effectsize_path(shared_data$data_dir, ds, st[[1]])
    shiny::validate(need(file.exists(path), paste(
      "No precomputed effect sizes for", st[[1]], "in", ds, "data."
    )))
    cpt_read_effectsizes(path)
  })

  # ---- Browse: dependency -----------------------------------------------

  # Both selection rules, applied to the whole subtype so the excluded side is
  # visible: effect size, p < 0.05, and the common-essential exclusion
  # (Avg < -0.5) that Discover applies on top of them.
  dep_browse <- reactive({
    df <- dep_universe()
    cut <- shared_data$dep_effect_min()
    # "Restrict to genes" governs the browse views too, so the tiles, the
    # charts and the table all describe the same rows.
    want <- wanted_genes()
    if (!is.null(want)) {
      df <- df[cpt_gene_match_key(df$gene_name) %in% want, , drop = FALSE]
    }
    df$p_value <- suppressWarnings(as.numeric(df$p_value))
    df$neglog_p <- -log10(pmax(df$p_value, .Machine$double.xmin))
    df$is_essential <- !is.na(df$Avg) & df$Avg < -0.5
    df$passes_effect <- !is.na(df$EffectSize) & df$EffectSize <= cut
    df$passes_p <- !is.na(df$p_value) & df$p_value < 0.05
    df$status <- ifelse(
      df$is_essential, "Common essential",
      ifelse(df$passes_effect & df$passes_p, "Selected",
             "Not selected"))
    list(df = df, cut = cut)
  })

  dep_vb <- function(expr, subtitle, icon_name, colour) {
    renderValueBox({
      b <- dep_browse()
      valueBox(format(expr(b), big.mark = ","), subtitle,
               icon = icon(icon_name), color = colour)
    })
  }
  # A funnel only reads as a funnel if each tile is a strict subset of the one
  # before it, so the row reads as a funnel: genes scored, then those
  # clearing the effect cutoff, then those also clearing p, then those left
  # after pan-essentials are removed. Counting every gene below -0.5 here
  # would not be a subset of the previous tile and would overstate how many
  # the pan-essential rule removes.
  output$dep_vb_total <- dep_vb(function(b) nrow(b$df),
                                "genes scored in this subtype", "dna", "light-blue")
  output$dep_vb_effect <- dep_vb(function(b) sum(b$df$passes_effect),
                                 "clear the effect-size cutoff", "arrow-down", "aqua")
  output$dep_vb_both <- dep_vb(
    function(b) sum(b$df$passes_effect & b$df$passes_p),
    "also clear p < 0.05", "filter", "blue")
  output$dep_vb_essential <- renderValueBox({
    b <- dep_browse()
    sel <- sum(b$df$status == "Selected")
    dropped <- sum(b$df$passes_effect & b$df$passes_p & b$df$is_essential)
    valueBox(
      format(sel, big.mark = ","),
      paste0("selected (", format(dropped, big.mark = ","),
             " pan-essential removed)"),
      icon = icon("circle-check"), color = "navy")
  })

  output$tp_dep_volcano <- plotly::renderPlotly({
    b <- dep_browse()
    d <- b$df[!is.na(b$df$EffectSize) & is.finite(b$df$neglog_p), , drop = FALSE]
    shiny::validate(need(nrow(d) > 0, cpt_empty_msg(
      "No scored genes for this subtype.", !is.null(wanted_genes()))))
    d$hover <- paste0(
      "<b>", d$gene_name, "</b><br>",
      "effect ", sprintf("%.3f", d$EffectSize),
      "<br>p ", signif(d$p_value, 3),
      "<br>mean across all lines ", sprintf("%.3f", d$Avg),
      "<br>", d$status)
    pal <- c("Selected" = CPT_PAL$ink, "Common essential" = CPT_PAL$light,
             "Not selected" = CPT_PAL$muted)
    # plot_ly() with type/mode set creates an empty placeholder trace before
    # the add_trace() calls; start bare, as the ligandability scatter does.
    p <- plotly::plot_ly(source = "cpt_dep_volcano")
    # Draw the greyed-out majority first so selected genes sit on top of it.
    for (lvl in c("Not selected", "Common essential", "Selected")) {
      dd <- d[d$status == lvl, , drop = FALSE]
      if (!nrow(dd)) next
      p <- plotly::add_trace(
        p, x = dd$EffectSize, y = dd$neglog_p, name = lvl,
        customdata = dd$gene_name,
        type = "scattergl", mode = "markers",
        marker = list(size = CPT_MARKER_SIZE,
                      opacity = if (lvl == "Not selected") 0.45 else 0.85,
                      color = unname(pal[[lvl]]),
                      line = list(width = 0)),
        text = dd$hover, hoverinfo = "text")
    }
    cpt_plotly(
      p,
      xaxis = cpt_axis("Gene effect size (more negative = stronger dependency)"),
      yaxis = cpt_axis("-log10 p"),
      shapes = list(
        list(type = "rect", layer = "below",
             x0 = min(d$EffectSize, na.rm = TRUE), x1 = b$cut,
             y0 = -log10(0.05), y1 = max(d$neglog_p, na.rm = TRUE),
             fillcolor = "rgba(71, 142, 184, 0.09)", line = list(width = 0)),
        list(type = "line", x0 = b$cut, x1 = b$cut, y0 = 0, y1 = 1,
             yref = "paper", line = list(color = CPT_PAL$alert, dash = "dash", width = 1.5)),
        list(type = "line", x0 = 0, x1 = 1, xref = "paper",
             y0 = -log10(0.05), y1 = -log10(0.05),
             line = list(color = CPT_PAL$alert, dash = "dash", width = 1.5))
      ),
      legend = list(orientation = "h", y = -0.18),
      margin = list(t = 20)
    )
  })

  # A histogram says how many genes sit where; the volcano cannot, because
  # 13,000 overlapping points hide their own density.
  output$tp_dep_hist <- plotly::renderPlotly({
    b <- dep_browse()
    d <- b$df[!is.na(b$df$EffectSize), , drop = FALSE]
    shiny::validate(need(nrow(d) > 0, cpt_empty_msg(
      "No scored genes for this subtype.", !is.null(wanted_genes()))))
    brk <- pretty(range(d$EffectSize), n = 60)
    w <- diff(brk)[1]
    kept <- d$EffectSize <= b$cut
    p <- plotly::plot_ly()
    p <- cpt_add_hist(p, d$EffectSize[!kept], "Excluded by the cutoff",
                      CPT_PAL$muted, w)
    p <- cpt_add_hist(p, d$EffectSize[kept], "Clears the cutoff",
                      CPT_PAL$ink, w)
    cpt_plotly(
      p, barmode = "stack",
      xaxis = cpt_axis("Gene effect size"),
      yaxis = cpt_axis("Genes"),
      shapes = list(list(type = "line", x0 = b$cut, x1 = b$cut, y0 = 0, y1 = 1,
                         yref = "paper",
                         line = list(color = CPT_PAL$alert, dash = "dash",
                                     width = 1.5))),
      legend = list(orientation = "h", y = -0.18), margin = list(t = 20))
  })

  # Selectivity directly: a gene on the diagonal is equally essential
  # everywhere, so distance below it is what "cancer-selective" means.
  output$tp_dep_vs <- plotly::renderPlotly({
    b <- dep_browse()
    d <- b$df
    shiny::validate(need(all(c("Cancer_Avg", "Other_Avg") %in% names(d)),
                         "This subtype's file has no per-group means."))
    d <- d[!is.na(d$Cancer_Avg) & !is.na(d$Other_Avg), , drop = FALSE]
    shiny::validate(need(nrow(d) > 0, cpt_empty_msg(
      "No scored genes for this subtype.", !is.null(wanted_genes()))))
    d$hover <- paste0("<b>", d$gene_name, "</b><br>this subtype ",
                      sprintf("%.3f", d$Cancer_Avg),
                      "<br>other lines ", sprintf("%.3f", d$Other_Avg),
                      "<br>difference ", sprintf("%.3f", d$Cancer_Avg - d$Other_Avg),
                      "<br>", d$status)
    pal <- c("Selected" = CPT_PAL$ink, "Common essential" = CPT_PAL$light,
             "Not selected" = CPT_PAL$muted)
    rng <- range(c(d$Cancer_Avg, d$Other_Avg), na.rm = TRUE)
    p <- plotly::plot_ly(source = "cpt_dep_vs")
    for (lvl in c("Not selected", "Common essential", "Selected")) {
      dd <- d[d$status == lvl, , drop = FALSE]
      if (!nrow(dd)) next
      p <- plotly::add_trace(
        p, x = dd$Other_Avg, y = dd$Cancer_Avg, name = lvl,
        customdata = dd$gene_name,
        type = "scattergl", mode = "markers",
        marker = list(size = CPT_MARKER_SIZE,
                      opacity = if (lvl == "Not selected") 0.35 else 0.85,
                      color = unname(pal[[lvl]]), line = list(width = 0)),
        text = dd$hover, hoverinfo = "text")
    }
    cpt_plotly(
      p,
      xaxis = cpt_axis("Mean effect in the other cell lines"),
      yaxis = cpt_axis("Mean effect in this subtype"),
      shapes = list(list(type = "line", x0 = rng[1], x1 = rng[2],
                         y0 = rng[1], y1 = rng[2],
                         line = list(color = "#b0bcc6", dash = "dot", width = 1))),
      legend = list(orientation = "h", y = -0.18), margin = list(t = 20))
  })

  dep_browse_table <- reactive({
    b <- dep_browse()
    d <- b$df
    out <- data.frame(
      Gene = d$gene_name,
      Status = d$status,
      `Effect size` = round(d$EffectSize, 3),
      `p` = signif(d$p_value, 3),
      `q` = if ("q.value" %in% names(d)) signif(d$q.value, 3) else NA_real_,
      `-log10 p` = if ("neg_log10_p_value" %in% names(d))
        round(d$neg_log10_p_value, 2) else NA_real_,
      # The summary beside this table reported the cohort means and the
      # non-cancer contrast; this one stopped short of them although the
      # per-subtype effect-size file carries both.
      `Mean in subtype` = round(d$Cancer_Avg, 3),
      `Mean in non-cancer` = if ("NonCancer_Avg" %in% names(d))
        round(d$NonCancer_Avg, 3) else NA_real_,
      `Mean in other lines` = round(d$Other_Avg, 3),
      `Mean across all lines` = round(d$Avg, 3),
      `p vs non-cancer` = if ("pval_vs_NonCancer" %in% names(d))
        signif(d$pval_vs_NonCancer, 3) else NA_real_,
      check.names = FALSE, stringsAsFactors = FALSE
    )
    out[order(out$`Effect size`, na.last = TRUE), , drop = FALSE]
  })

  output$tp_dep_table <- DT::renderDT({
    DT::datatable(
      # Status has three values; as a character column DT gave it a blank text
      # box, as a factor it gives a dropdown of the three.
      cpt_filter_levels(dep_browse_table()),
      rownames = FALSE, filter = "top", selection = "none",
      options = list(pageLength = 15, scrollX = TRUE,
                     order = list(list(2, "asc")),
                     language = list(emptyTable = cpt_empty_msg(
                       "No scored genes for this subtype.",
                       !is.null(wanted_genes()))))
    ) |>
      DT::formatStyle("Status", target = "row",
                      color = DT::styleEqual("Not selected", "#7b8a94",
                                             default = NULL))
  })

  output$dl_dep_browse <- downloadHandler(
    filename = function() paste0("canprotarget_dependency_", Sys.Date(), ".xlsx"),
    content = function(file) {
      df <- dep_browse_table()
      rows <- input$tp_dep_table_rows_all
      if (!is.null(rows) && length(rows)) df <- df[rows, , drop = FALSE]
      # The candidates export records the settings that produced it; these
      # did not, so a downloaded sheet could not be traced back to a run.
      writexl::write_xlsx(
        list(dependency = df,
             provenance = data.frame(
               setting = c("dataset", "subtype", "effect-size cutoff",
                           "p-value rule", "common-essential rule",
                           "rows exported", "exported"),
               value = c(get_dataset(),
                         paste(get_subtype(), collapse = "; "),
                         as.character(shared_data$dep_effect_min()),
                         "p < 0.05",
                         "mean effect < -0.5 excluded",
                         as.character(nrow(df)),
                         format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
               stringsAsFactors = FALSE)),
        path = file)
    }
  )

  # ---- Browse: ligandability ---------------------------------------------

  lig_browse <- reactive({
    # Same source as the candidate set, so the browse views and the summary
    # describe the same population at any competition-ratio floor.
    es <- engaged_sites_for_cutoff()
    req(!is.null(es), nrow(es) > 0)
    cr <- input$tp_cr %||% 4
    mt <- input$tp_max_targets %||% 20

    want <- wanted_genes()
    if (!is.null(want)) {
      es <- es[cpt_gene_match_key(es$gene_name) %in% want, , drop = FALSE]
    }

    # Row grain means the same thing here as on the candidates table: what one
    # row stands for. Records are stored per cysteine x probe, so the coarser
    # grains keep the strongest engagement for each cysteine or gene.
    grain <- input$tp_grain %||% "probe"
    # Counted before the dedup below, which keeps only the strongest probe.
    if (nrow(es)) {
      ckey <- if (grain == "gene") cpt_gene_match_key(es$gene_name) else es$cysteineid
      n_by <- tapply(as.character(es$probe_name), ckey, function(x) length(unique(x)))
      es$n_probes <- as.integer(n_by[as.character(ckey)])
    }
    if (nrow(es) && grain != "probe") {
      es <- es[order(-ifelse(is.na(es$CR), -Inf, es$CR)), , drop = FALSE]
      key <- if (grain == "gene") cpt_gene_match_key(es$gene_name) else es$cysteineid
      es <- es[!duplicated(key), , drop = FALSE]
    }

    # SMILES, dataset and cell line are measurements, so they live in the
    # binding table rather than the index. The CR >= 4 file holds exactly these
    # records, so the join is cheap.
    cr4 <- tryCatch(shared_data$protein_binding_cr4(), error = function(e) NULL)
    if (is.data.frame(cr4) && nrow(cr4) && nrow(es)) {
      k <- paste(tolower(as.character(cr4$probe_name)), as.character(cr4$cysteineid))
      m <- match(paste(tolower(as.character(es$probe_name)), as.character(es$cysteineid)), k)
      for (nm in c("SMILES", "Dataset", "Cell_Line")) {
        if (nm %in% names(cr4)) es[[nm]] <- as.character(cr4[[nm]])[m]
      }
    }

    es$passes <- !is.na(es$CR) & es$CR >= cr &
      (is.na(es$n_targets) | es$n_targets <= mt)
    list(es = es, cr = cr, mt = mt, grain = grain)
  })

  lig_vb <- function(expr, subtitle, icon_name, colour) {
    renderValueBox({
      b <- lig_browse()
      valueBox(format(expr(b), big.mark = ","), subtitle,
               icon = icon(icon_name), color = colour)
    })
  }
  output$lig_vb_sites <- lig_vb(function(b) nrow(b$es),
                                "engaged cysteine records", "atom", "light-blue")
  output$lig_vb_probes <- lig_vb(function(b) length(unique(b$es$probe_name)),
                                 "distinct covalent probes", "flask", "aqua")
  output$lig_vb_proteins <- lig_vb(
    function(b) length(unique(cpt_gene_match_key(b$es$gene_name))),
    "proteins engaged", "dna", "blue")
  # Unlike the dependency tiles this row is a summary, not a funnel: records,
  # probes and proteins are three different denominators. Only this tile is a
  # subset of the first, so it says which one it is a subset of.
  output$lig_vb_pass <- renderValueBox({
    b <- lig_browse()
    valueBox(
      format(sum(b$es$passes), big.mark = ","),
      paste0("of ", format(nrow(b$es), big.mark = ","),
             " records clear both cutoffs"),
      icon = icon("circle-check"), color = "navy")
  })

  output$tp_lig_scatter <- plotly::renderPlotly({
    b <- lig_browse()
    p <- cpt_cr_scatter(b$es, b$cr, b$mt,
                        by_tier = "functionality" %in% layers())
    shiny::validate(need(!is.null(p), cpt_empty_msg(
      "No engaged cysteines in the index.", !is.null(wanted_genes()))))
    p
  })

  # The competition ratio is a ratio bounded at 20 by the assay; a handful of
  # stored records exceed it. Say so rather than clipping them out of the plot.
  output$tp_lig_range_note <- renderUI({
    b <- lig_browse()
    n <- sum(!is.na(b$es$CR) & b$es$CR > 20)
    if (!n) return(NULL)
    cpt_note(
      paste0(n, " record", if (n == 1) "" else "s",
             " report a competition ratio above 20."),
      "A competition ratio is bounded at 20 by the assay, so these values are ",
      "out of range. They are plotted and tabulated as stored, not clipped, ",
      "so you can see them; treat them as saturated rather than as stronger ",
      "engagement than a CR of 20.",
      tone = "caution"
    )
  })

  # Where the competition-ratio cutoff actually bites.
  output$tp_lig_cr_hist <- plotly::renderPlotly({
    b <- lig_browse()
    es <- b$es[!is.na(b$es$CR), , drop = FALSE]
    shiny::validate(need(nrow(es) > 0, cpt_empty_msg(
      "No engaged cysteines in the index.", !is.null(wanted_genes()))))
    v <- pmin(es$CR, 20)
    p <- plotly::plot_ly()
    p <- cpt_add_hist(p, v[v < b$cr], "Below the cutoff", CPT_PAL$muted, 0.5)
    p <- cpt_add_hist(p, v[v >= b$cr], "Clears the cutoff", CPT_PAL$ink, 0.5)
    cpt_plotly(
      p, barmode = "stack",
      xaxis = cpt_axis("Competition ratio (values above 20 counted at 20)"),
      yaxis = cpt_axis("Engaged cysteine records"),
      shapes = list(list(type = "line", x0 = b$cr, x1 = b$cr, y0 = 0, y1 = 1,
                         yref = "paper",
                         line = list(color = CPT_PAL$alert, dash = "dash",
                                     width = 1.5))),
      legend = list(orientation = "h", y = -0.18), margin = list(t = 20))
  })

  # One bar per probe, not per record: how promiscuous the chemistry is, and
  # how few probes survive a 20-protein limit.
  output$tp_lig_probe_hist <- plotly::renderPlotly({
    b <- lig_browse()
    es <- b$es[!is.na(b$es$n_targets), , drop = FALSE]
    shiny::validate(need(nrow(es) > 0, cpt_empty_msg(
      "No engaged cysteines in the index.", !is.null(wanted_genes()))))
    per_probe <- tapply(es$n_targets, es$probe_name, function(x) x[1])
    v <- log10(pmax(as.numeric(per_probe), 1))
    keep <- as.numeric(per_probe) <= b$mt
    p <- plotly::plot_ly()
    p <- cpt_add_hist(p, v[!keep], "Too promiscuous", CPT_PAL$muted, 0.1)
    p <- cpt_add_hist(p, v[keep], "Within the limit", CPT_PAL$ink, 0.1)
    cpt_plotly(
      p, barmode = "stack",
      xaxis = cpt_axis("Proteins engaged by the probe",
                       tickmode = "array",
                       tickvals = log10(c(1, 5, 10, 20, 50, 100, 250, 500, 1000)),
                       ticktext = c("1", "5", "10", "20", "50", "100", "250",
                                    "500", "1,000")),
      yaxis = cpt_axis("Probes (one bar counts probes, not records)"),
      shapes = list(list(type = "line", x0 = log10(b$mt), x1 = log10(b$mt),
                         y0 = 0, y1 = 1, yref = "paper",
                         line = list(color = CPT_PAL$alert, dash = "dash",
                                     width = 1.5))),
      legend = list(orientation = "h", y = -0.18), margin = list(t = 20))
  })

  # How much of the engaged set carries functional evidence at all, split by
  # whether it clears the cutoffs. Tier 4 dominates; that is the point.
  output$tp_lig_tiers <- plotly::renderPlotly({
    b <- lig_browse()
    es <- b$es
    shiny::validate(need(nrow(es) > 0, cpt_empty_msg(
      "No engaged cysteines in the index.", !is.null(wanted_genes()))))
    tab <- table(factor(es$evidence_tier, levels = 1:4),
                 factor(ifelse(es$passes, "Clears both cutoffs", "Filtered out"),
                        levels = c("Clears both cutoffs", "Filtered out")))
    labs <- paste0("Tier ", 1:4)
    meaning <- unname(CPT_TIER_LABELS[as.character(1:4)])
    p <- plotly::plot_ly()
    p <- plotly::add_trace(
      p, x = as.numeric(tab[, "Clears both cutoffs"]), y = labs,
      name = "Clears both cutoffs", type = "bar", orientation = "h",
      marker = list(color = CPT_PAL$ink), customdata = meaning,
      hovertemplate = paste0("<b>%{y}</b>: %{customdata}",
                             "<br>%{x:,} records clear both cutoffs<extra></extra>"))
    p <- plotly::add_trace(
      p, x = as.numeric(tab[, "Filtered out"]), y = labs,
      name = "Filtered out", type = "bar", orientation = "h",
      marker = list(color = CPT_PAL$muted), customdata = meaning,
      hovertemplate = paste0("<b>%{y}</b>: %{customdata}",
                             "<br>%{x:,} records filtered out<extra></extra>"))
    cpt_plotly(
      p, barmode = "stack",
      xaxis = cpt_axis("Engaged cysteine records"),
      yaxis = cpt_axis(""),
      legend = list(orientation = "h", y = -0.18),
      margin = list(t = 20, l = 70))
  })

  # Dependency x cysteine function. The two layers meet on one question: among
  # the genes this subtype depends on, how strong is the dependency at each
  # grade of residue evidence? A box per tier answers it directly, where the
  # per-layer charts can only be read side by side.
  output$tp_dep_tier <- plotly::renderPlotly({
    d <- enriched_rows()
    shiny::validate(need(
      is.data.frame(d) && nrow(d) > 0 &&
        all(c("effect_size", "evidence_tier") %in% names(d)),
      cpt_empty_msg("No rows carry both a dependency and a tier.",
                    !is.null(wanted_genes()))))
    d <- d[!is.na(d$effect_size) & !is.na(d$evidence_tier), , drop = FALSE]
    shiny::validate(need(nrow(d) > 0,
      "No rows carry both a dependency and a tier."))

    tiers <- sort(unique(as.integer(d$evidence_tier)))
    ti <- match(as.integer(d$evidence_tier), tiers)
    site <- if ("site" %in% names(d)) d$site else rep("", nrow(d))

    # Points are spread across the box by rank within their tier rather than at
    # random: no RNG to seed, and the same data always draws the same chart.
    spread <- ave(seq_len(nrow(d)), ti, FUN = function(i) {
      n <- length(i)
      if (n == 1L) return(0)
      seq(-0.17, 0.17, length.out = n)[rank(i, ties.method = "first")]
    })

    # The box draws the distribution and the scatter draws the genes. A box
    # trace with boxpoints="all" hovers its own five-number summary over the
    # points, so the gene behind a dot could not be read at all.
    p <- plotly::plot_ly()
    for (k in seq_along(tiers)) {
      rows <- d[ti == k, , drop = FALSE]
      p <- plotly::add_trace(
        p, y = rows$effect_size, x = rep(k, nrow(rows)),
        type = "box", boxpoints = FALSE, hoverinfo = "skip", width = 0.5,
        line = list(color = CPT_PAL$ink), fillcolor = "rgba(71,142,184,0.18)",
        showlegend = FALSE)
    }
    p <- plotly::add_trace(
      p, type = "scatter", mode = "markers",
      x = ti + spread, y = d$effect_size,
      marker = list(size = 6, color = CPT_PAL$primary, opacity = 0.7,
                    line = list(width = 0.5, color = "#ffffff")),
      showlegend = FALSE,
      text = paste0("<b>", as.character(d$gene_name), "</b>",
                    ifelse(nzchar(site), paste0(" ", site), ""),
                    "<br>Tier ", as.integer(d$evidence_tier), " \u2014 ",
                    unname(CPT_TIER_LABELS[as.character(d$evidence_tier)])),
      hovertemplate = "%{text}<br>Gene effect %{y:.3f}<extra></extra>")

    cpt_plotly(
      p,
      xaxis = cpt_axis("", tickmode = "array", tickvals = seq_along(tiers),
                       ticktext = paste("Tier", tiers),
                       range = c(0.4, length(tiers) + 0.6)),
      yaxis = cpt_axis("Gene effect in this subtype"),
      showlegend = FALSE,
      shapes = list(list(type = "line", xref = "paper", x0 = 0, x1 = 1,
                         y0 = -0.5, y1 = -0.5,
                         line = list(dash = "dash", color = CPT_PAL$alert,
                                     width = 1))),
      margin = list(t = 20, l = 70))
  })

  lig_browse_table <- reactive({
    b <- lig_browse()
    es <- b$es
    cid <- as.character(es$cysteineid)
    # A missing column stands in as NA for every row, not as a single NA.
    # data.frame() recycles length one against a non-empty frame but not
    # against an empty one, so the scalar form threw the moment a filter
    # emptied the selection.
    col <- function(nm, fill = NA_character_) {
      if (nm %in% names(es)) es[[nm]] else rep(fill, nrow(es))
    }
    out <- data.frame(
      Gene = es$gene_name,
      Protein = sub("_C.*$", "", cid),
      Site = sub("^.*_C", "C", cid),
      Probe = es$probe_name,
      SMILES = col("SMILES"),
      CR = round(es$CR, 2),
      `Probe targets` = es$n_targets,
      `Probes available` = col("n_probes", NA_integer_),
      Dataset = col("Dataset"),
      `Cell line` = col("Cell_Line"),
      `Passes cutoffs` = ifelse(es$passes, "yes", "no"),
      check.names = FALSE, stringsAsFactors = FALSE
    )
    # No atlas columns here. This is the ligandability layer's own table, and
    # a per-layer view shows that layer only whatever else is switched on --
    # the joins are what the Combined row is for.
    #
    # Prune only when there is something to prune: all(is.na(x)) is TRUE for a
    # zero-length column, so on an empty selection this dropped every column.
    if (nrow(out)) {
      out <- out[, !vapply(out, function(x) all(is.na(x)), logical(1)), drop = FALSE]
    }
    if (!nrow(out) || !all(c("CR", "Probe targets") %in% names(out))) return(out)
    out[order(-pmin(out$CR, 20), out$`Probe targets`, na.last = TRUE), ,
        drop = FALSE]
  })

  output$tp_lig_table <- DT::renderDT({
    DT::datatable(
      lig_browse_table(),
      rownames = FALSE, filter = "top", selection = "none",
      options = list(pageLength = 15, scrollX = TRUE,
                     language = list(emptyTable = cpt_empty_msg(
                       "No engaged cysteines in the index.",
                       !is.null(wanted_genes()))))
    ) |>
      DT::formatStyle("Passes cutoffs", target = "row",
                      color = DT::styleEqual("no", "#7b8a94", default = NULL)) |>
      DT::formatStyle("CR",
                      color = DT::styleInterval(20, c(NA, "#cc3340")),
                      fontWeight = DT::styleInterval(20, c(NA, "bold")))
  })

  output$dl_lig_browse <- downloadHandler(
    filename = function() paste0("canprotarget_ligandability_", Sys.Date(), ".xlsx"),
    content = function(file) {
      df <- lig_browse_table()
      rows <- input$tp_lig_table_rows_all
      if (!is.null(rows) && length(rows)) df <- df[rows, , drop = FALSE]
      writexl::write_xlsx(
        list(ligandability = df,
             provenance = data.frame(
               setting = c("minimum competition ratio", "maximum probe targets",
                           "index floor", "rows exported", "exported"),
               value = c(as.character(input$tp_cr %||% 4),
                         as.character(input$tp_max_targets %||% 20),
                         "engagement records are stored at CR >= 4 only",
                         as.character(nrow(df)),
                         format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
               stringsAsFactors = FALSE)),
        path = file)
    }
  )

  output$tp_download <- downloadHandler(
    filename = function() paste0("canprotarget_candidates_", Sys.Date(), ".xlsx"),
    content = function(file) {
      df <- display_df()
      rows <- input$tp_table_rows_all
      # Export what is displayed, including any column filters the user typed.
      if (!is.null(rows) && length(rows)) df <- df[rows, , drop = FALSE]
      writexl::write_xlsx(
        list(candidates = df,
             filters = data.frame(
               setting = c("layers", "grain", "min CR", "max probe targets",
                           "tiers", "require ClinVar", "subtype", "dataset"),
               value = c(paste(input$tp_layers, collapse = ", "),
                         input$tp_grain, input$tp_cr, input$tp_max_targets,
                         paste(input$tp_tiers, collapse = ", "),
                         isTRUE(input$tp_clinvar),
                         paste(get_subtype(), collapse = "; "), get_dataset()),
               stringsAsFactors = FALSE)),
        file
      )
    }
  )

  # ---- Residue analysis ---------------------------------------------------

  # Same rows the ligandability browse works from, but read residue-first and
  # without the CR/target cutoffs, which are a chemistry filter: a tested,
  # non-functional residue is still residue evidence.
  res_rows <- reactive({
    b <- lig_browse()
    b$es
  })

  output$res_scope_note <- renderUI({
    b <- lig_browse()
    want <- wanted_genes()
    tags$div(
      class = "cpt-scope-note",
      paste0(
        "One row per ",
        switch(b$grain, gene = "gene", probe = "cysteine × probe", "cysteine"),
        if (is.null(want)) "" else paste0(", restricted to ", length(want), " gene(s)"),
        ". Tier is a property of the engaged residue, so the chemistry cutoffs ",
        "are reported here but not applied."
      )
    )
  })

  res_vb <- function(expr, subtitle, icon_name, colour) {
    renderValueBox({
      d <- res_rows()
      valueBox(format(expr(d), big.mark = ","), subtitle,
               icon = icon(icon_name), color = colour)
    })
  }
  output$res_vb_sites <- res_vb(function(d) nrow(d),
                                "engaged cysteines in scope", "atom", "light-blue")
  output$res_vb_functional <- res_vb(
    function(d) sum(isTRUE_vec(d$cys_functional)),
    "functional in the atlas", "bullseye", "aqua")
  output$res_vb_tier1 <- res_vb(
    function(d) sum(d$evidence_tier == 1, na.rm = TRUE),
    "Tier 1 (functional + ligandable)", "circle-check", "navy")
  output$res_vb_clinvar <- res_vb(
    function(d) sum(isTRUE_vec(d$cys_clinvar_pathogenic)),
    "ClinVar pathogenic", "triangle-exclamation", "blue")

  output$res_tier_chart <- plotly::renderPlotly({
    d <- res_rows()
    shiny::validate(need(nrow(d) > 0, cpt_empty_msg(
      "No engaged cysteines in scope.", !is.null(wanted_genes()))))
    n <- vapply(1:4, function(t) sum(d$evidence_tier == t, na.rm = TRUE), numeric(1))
    labs <- paste("Tier", 1:4)
    p <- plotly::plot_ly(
      x = labs, y = n, type = "bar",
      marker = list(color = unname(CPT_TIER_PAL[as.character(1:4)])),
      customdata = unname(CPT_TIER_LABELS[as.character(1:4)]),
      hovertemplate = "<b>%{x}</b>: %{customdata}<br>%{y:,} cysteines<extra></extra>")
    cpt_plotly(p, xaxis = cpt_axis(""), yaxis = cpt_axis("Engaged cysteines"),
               showlegend = FALSE, margin = list(t = 20))
  })

  # How the residue was matched to the atlas decides how much the tier is
  # worth, so it deserves its own view rather than a column buried in a table.
  output$res_join_chart <- plotly::renderPlotly({
    d <- res_rows()
    shiny::validate(need(nrow(d) > 0, cpt_empty_msg(
      "No engaged cysteines in scope.", !is.null(wanted_genes()))))
    st <- ifelse(is.na(d$cys_join_status), "unknown", d$cys_join_status)
    tb <- sort(table(st), decreasing = TRUE)
    p <- plotly::plot_ly(
      x = as.numeric(tb), y = names(tb), type = "bar", orientation = "h",
      marker = list(color = CPT_PAL$primary),
      hovertemplate = "%{y}<br>%{x:,} cysteines<extra></extra>")
    cpt_plotly(p, xaxis = cpt_axis("Engaged cysteines"), yaxis = cpt_axis(""),
               showlegend = FALSE, margin = list(t = 20, l = 220))
  })

  res_table_df <- reactive({
    d <- res_rows()
    if (!nrow(d)) return(data.frame(`No rows` = character(0), check.names = FALSE))
    cid <- as.character(d$cysteineid)
    out <- data.frame(
      # Residue first: this tab is about the site, not the gene.
      Site = sub("^.*_C", "C", cid),
      Protein = sub("_C.*$", "", cid),
      Gene = d$gene_name,
      Tier = d$evidence_tier,
      `Tier meaning` = d$evidence_tier_label,
      `Functional (atlas)` = d$cys_functional,
      `Atlas ligandability` = d$cys_ligandable,
      ClinVar = d$cys_clinvar_pathogenic,
      `Editor support` = d$cys_editor_support,
      `Study context` = d$cys_study_context,
      `Matched on` = d$cys_join_status,
      check.names = FALSE, stringsAsFactors = FALSE
    )
    # No probe columns here, for the same reason the ligandability table
    # carries no atlas columns: this is the cysteine function layer's own view.
    lgl <- vapply(out, is.logical, logical(1))
    out[lgl] <- lapply(out[lgl], function(x) {
      ifelse(is.na(x), "", ifelse(x, "Yes", "—"))
    })
    # CR is present only when the ligandability layer is on, so it can only be
    # a sort key then. Passing a NULL column to order() makes the argument
    # lengths differ, which threw and left the table blank on a residue-only run.
    ord <- if ("CR" %in% names(out)) {
      order(out$Tier, -ifelse(is.na(out$CR), -Inf, out$CR))
    } else {
      order(out$Tier)
    }
    out[ord, , drop = FALSE]
  })

  output$res_table <- DT::renderDT({
    DT::datatable(
      res_table_df(), rownames = FALSE, filter = "top", selection = "none",
      options = list(pageLength = 15, scrollX = TRUE,
                     language = list(emptyTable = cpt_empty_msg(
                       "No engaged cysteines in scope.",
                       !is.null(wanted_genes()))))
    )
  })

  output$dl_res_analysis <- downloadHandler(
    filename = function() paste0("canprotarget_residues_", Sys.Date(), ".xlsx"),
    content = function(file) {
      df <- res_table_df()
      rows <- input$res_table_rows_all
      if (!is.null(rows) && length(rows)) df <- df[rows, , drop = FALSE]
      writexl::write_xlsx(
        list(residues = df,
             provenance = data.frame(
               setting = c("row grain", "gene restriction", "tier source",
                           "rows exported", "exported"),
               value = c(input$tp_grain %||% "probe",
                         if (is.null(wanted_genes())) "none"
                         else paste(wanted_genes(), collapse = ", "),
                         "cpt_evidence_tier() on engaged cysteines (CR >= 4)",
                         as.character(nrow(df)),
                         format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
               stringsAsFactors = FALSE)),
        path = file)
    }
  )

  # A scientist who spots an interesting point wants that gene, not a hover.
  # Clicking any point routes to the Target tab, the same destination the
  # candidates drawer offers.
  chart_gene <- reactiveVal(NULL)
  for (src in c("cpt_dep_volcano", "cpt_dep_vs", "cpt_cr_scatter")) {
    local({
      this_src <- src
      observeEvent(plotly::event_data("plotly_click", source = this_src), {
        ev <- plotly::event_data("plotly_click", source = this_src)
        key <- tryCatch(as.character(ev$customdata)[[1]],
                        error = function(e) NULL)
        if (is.null(key) || !length(key) || is.na(key) || !nzchar(key)) return()
        key <- cpt_gene_match_key(key)
        chart_gene(key)
        if (is.function(on_open_gene)) on_open_gene(key)
      }, ignoreInit = TRUE)
    })
  }

  # Selected gene, so the host can route "Open in Gene tab".
  reactive({
    cg <- chart_gene()
    if (!is.null(cg) && nzchar(cg)) return(cg)
    sel <- drawer_row()
    if (is.null(sel) || !length(sel) || is.na(sel)) return(NULL)
    df <- tryCatch(candidates(), error = function(e) NULL)
    if (is.null(df) || nrow(df) < sel[[1]]) return(NULL)
    df$gene_key[sel[[1]]]
  })
}
