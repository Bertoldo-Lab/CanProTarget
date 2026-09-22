# ============================================================
# Script:   dependencies_module.R
# Purpose:  Shiny module for Cancer Dependencies tab:
#           - run gene effect analysis (DepMap-style)
#           - link dependencies to chemical probes
#           - provide plots + export tables
# Inputs:   shared_data from app.R (see below)
# Outputs:  UI + server functions for tabName = "gene_effect"
# ============================================================
# Sections:
#   1. dependencies_ui()
#   2. dependencies_server()
#      2.1 Initialization and dataset selection
#      2.2 Gene effect matrix loading
#      2.3 Cancer dependency analysis
#      2.4 Gene tables and exports
#      2.5 Volcano plot
#      2.6 Group comparison plot
#      2.7 Probe analysis and tables/plots
# ============================================================

library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(shinyjs)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)
library(DT)
library(writexl)
library(readr)

# ============================================================
# 1. dependencies_ui
# ============================================================

# sanitize_subtype() is defined in R/functions.R (sourced before this module)

# One gene-table panel: ranked table + volcano, optionally a per-gene group
# plot, shared by "Cancer-Selective Genes" and "All Genes";
# they differ only in which data frame they show and whether the group plot
# belongs with them. Input/output ids are unchanged so bookmarks still resolve.
gene_panel_ui <- function(ns, key, title, table_help, volcano_help) {
  # No box() here: this sits inside Discover's shared table space, and a box
  # within a box reads as a stray panel.
  panel <- list(
    tagList(
      do.call(cpt_section, c(list(title), table_help)),
      shinycssloaders::withSpinner(
        DTOutput(ns(paste0(key, "_genes_table"))),
        type = 6, color = "#478EB8"),
      br(),
      checkboxInput(ns(paste0("show_all_cols_", key)), "Show all columns", FALSE),
      downloadButton(ns(paste0("dl_", key, "_genes")),
        "Download (.xlsx)", class = "btn-sm btn-default")
    )
  )

  volcano_box <- box(
    width = 12, status = "info", solidHeader = FALSE,
    do.call(cpt_section, c(list("Volcano Plot"), volcano_help)),
    shinycssloaders::withSpinner(
      plotlyOutput(ns(paste0("volcano_plot_", key)), height = "400px"),
      type = 6, color = "#478EB8"),
    br(),
    downloadButton(ns(paste0("dl_volcano_", key)),
      "Download PNG", class = "btn-sm btn-default")
  )

  c(panel, list(volcano_box))
}

dependencies_ui <- function(id) {
  ns <- NS(id)

  # Three layer boxes of equal width, each owning its include-checkbox and its
  # own parameters, with row grain and the action beneath them. The layer
  # toggles sit below the three layer boxes rather than ahead of them; the dependency
  # parameters sat in a narrow sidebar, so the controls did not line up with
  # the layers they drive.
  tabItem(
    tabName = "gene_effect",
    div(
      class = "cpt-controls",
      fluidRow(
        column(
        width = 4, class = "cpt-layer-box", id = ns("box_dep"),
        box(
          title = "Dependency", width = 12, status = "primary", solidHeader = TRUE,
          actionButton(ns("tp_browse_dep"), "Browse the full dependency results",
                       icon = icon("chart-column"),
                       class = "btn-default btn-block btn-sm cpt-browse-btn"),
          div(class = "cpt-box-summary", textOutput(ns("dep_summary"), inline = TRUE)),
          tags$details(
            class = "cpt-box-details",
            tags$summary(class = "cpt-box-details-summary", "Change"),

        # Dataset selector (CRISPR vs RNAi)
        # RNAi is the primary discovery layer in the manuscript; CRISPR is the
        # orthogonal confirmatory layer. The default matches that, so the app
        # opens on the same footing as the published analysis.
        selectInput(ns("dataset"), "Gene Effect Dataset",
          choices  = c("RNAi"          = "RNAi",
                       "CRISPR (23Q4)" = "CRISPR (23Q4)"),
          selected = "RNAi"),

        # Cancer subtype picker (choices from data/cancer_subtypes_*.txt; see app_config)
        pickerInput(ns("cancer_subtypes"), "Cancer Subtype",
          choices = character(0),
          multiple = FALSE,
          options  = list(
            `live-search`  = TRUE,
            title          = "Select cancer"
          )),
        hr(),

        # One disclosure, not two: the thresholds are part of "Change", so
        # opening the box shows everything that decides the numbers rather than
        # hiding half of it behind a second collapse.
        cpt_section("Selection thresholds", level = "h5"),
        cpt_label("Effect size cutoff", for_id = ns("effect_min")),
        numericInput(ns("effect_min"), label = NULL, value = -0.1, step = 0.05),
        checkboxInput(ns("apply_pvalue"),
          "Apply linear-association p-value filter (< 0.05)", TRUE),
        checkboxInput(ns("excl_common"),
          "Exclude common essentials (Avg < -0.5)", TRUE),
        checkboxInput(ns("req_nc_sig"),
          "Require significance vs non-cancer (p < 0.05)", FALSE)
          ),
        )
        ),
        targets_panel_ligandability_box(ns),
        targets_panel_residue_box(ns)
      ),
      targets_panel_actions(ns)
    ),
    # Everything that produces a table or a chart goes through
    # targets_panel_results(): one table space, one chart space, the reader
    # chooses what is in them. The ranked dependency block and the atlas pieces
    # are the only outputs the parent cannot address from its own namespace, so
    # they are passed in.
    targets_panel_results(
      ns,
      dep_extra  = dep_ranked_ui(ns),
      atlas_table = cys_editing_table_ui("cys"),
      atlas_charts = list(
        abe_cbe  = cys_editing_chart_ui("cys", "abe_cbe"),
        overview = cys_editing_chart_ui("cys", "overview")
      )
    )
  )
}

#' The ranked cancer-selective table, shown once an analysis has been run.
dep_ranked_ui <- function(ns) {
  conditionalPanel(
    condition = paste0("output['", ns("deps_ready"), "'] == '1'"),
    uiOutput(ns("results_dataset_badge")),
    br(),
    gene_panel_ui(
      ns, "cancer", "Cancer-Selective Genes",
      table_help = list(
        "CPT_Score ranks dependency and selectivity against all genes in the ",
        "subtype. Cys_tier is the evidence tier of the gene's best engaged ",
        "cysteine (CR >= 4); Probes and Max_CR summarise covalent engagement."
      ),
      volcano_help = list()
    )
  )
}

# ============================================================
# 2. dependencies_server
# ============================================================

#' @param preset Optional reactive carrying list(subtype, dataset, gene) from
#'   the Target tab. When it fires, this tab configures itself and runs, so the
#'   user does not have to re-find the subtype they just clicked.
#' @param on_open_gene Optional callback(gene_key) used by the target
#'   explorer's drawer to route to the Target tab.
dependencies_server <- function(id, shared_data, preset = NULL, on_open_gene = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    displayed_dataset <- reactiveVal(NULL)
    displayed_subtype <- reactiveVal(NULL)   # display label (paste of subtypes)
    displayed_subtypes <- reactiveVal(NULL)  # character vector frozen at Run Analysis
    # Filters frozen at Run so volcano/badge match tables
    displayed_filters <- reactiveVal(NULL)
    deps_ready <- reactiveVal(FALSE)
    lig_ready <- reactiveVal(FALSE)
    analysis_generation <- reactiveVal(0L)  # bumps on each successful Run

    # Bookmark support: restore cancer_subtypes selection after choices load.
    # Shiny restores input values before dynamic choices are populated,
    # so we store the bookmarked subtype and apply it when the picker refreshes.
    bookmark_subtype <- reactiveVal(NULL)

    onRestored(function(state) {
      # state$input contains namespaced input values as restored from URL
      subtype_key <- ns("cancer_subtypes")
      if (!is.null(state$input[[subtype_key]])) {
        bookmark_subtype(state$input[[subtype_key]])
      }
    })

    # Surface missing metadata early (avoids silent spinner / empty picker)
    meta_missing_warned <- reactiveVal(FALSE)
    observe({
      meta <- shared_data$cancer_model_data()
      if (is.null(meta) && !isTRUE(meta_missing_warned())) {
        meta_missing_warned(TRUE)
        cpt_notify_error(
          "CPT-1001",
          "Cancer model metadata not available",
          detail = "cancer_model_data.rds missing or unreadable; analysis cannot run",
          type = "warning",
          duration = 15
        )
      }
    })

    output$deps_ready <- renderText({
      if (isTRUE(deps_ready())) "1" else "0"
    })
    outputOptions(output, "deps_ready", suspendWhenHidden = FALSE)

    output$lig_ready <- renderText({
      if (isTRUE(lig_ready())) "1" else "0"
    })
    outputOptions(output, "lig_ready", suspendWhenHidden = FALSE)

    refresh_cancer_subtype_picker <- function(dataset_label, preserve_selection) {
      smap <- shared_data$cancer_subtype_list_files
      fn <- smap[[dataset_label]]
      if (is.null(fn) || !nzchar(fn)) {
        updatePickerInput(session, "cancer_subtypes", choices = character(0), selected = NULL)
        return(invisible(NULL))
      }
      path <- shared_data$data_path(fn)
      choices <- read_cancer_subtype_list_file(path)
      if (!length(choices)) {
        showNotification(
          paste0(
            "No subtypes in ", basename(path), " for ", dataset_label, ". ",
            "Run docs/scripts/preprocess_data.R to generate lists (see docs/DATA_PROVENANCE.md), ",
            "or edit the file (one OncotreeSubtype per line; # starts a comment)."
          ),
          type = "warning",
          duration = 12
        )
      }
      # Check for bookmarked subtype first, then current selection
      bm_sub <- isolate(bookmark_subtype())
      cur_char <- if (!is.null(bm_sub) && nzchar(bm_sub)) {
        bookmark_subtype(NULL)  # consume it
        bm_sub
      } else if (isTRUE(preserve_selection)) {
        as.character(isolate(input$cancer_subtypes))
      } else {
        character(0)
      }
      sel <- if (length(cur_char) && cur_char[1] %in% choices) {
        cur_char[1]
      } else if (length(choices)) {
        choices[1]
      } else {
        NULL
      }
      updatePickerInput(session, "cancer_subtypes", choices = choices, selected = sel)
      invisible(NULL)
    }

    observeEvent(input$dataset, {
      refresh_cancer_subtype_picker(input$dataset, preserve_selection = TRUE)
    }, ignoreInit = FALSE)

    output$dep_summary <- renderText({
      st <- as.character(input$cancer_subtypes)
      paste0(
        input$dataset %||% "RNAi",
        " \u00b7 ", if (length(st) && nzchar(st[1])) st[1] else "no subtype",
        " \u00b7 effect \u2264 ", input$effect_min %||% -0.1,
        if (isTRUE(input$apply_pvalue)) " \u00b7 p<0.05" else "",
        if (isTRUE(input$excl_common)) " \u00b7 no common essentials" else ""
      )
    })

    # Layered target explorer. Shares this tab's subtype and dataset inputs so
    # one control panel drives both modes.
    # Radar, weight sensitivity and the HTML report follow the gene in the
    # drawer instead of occupying their own subtab. They are the only outputs
    # that need the run's subtype-and-weights context, so they appear only once
    # an analysis has been run.
    drawer_detail_ui <- function(gene_key) {
      if (!isTRUE(deps_ready())) {
        return(tags$p(class = "text-muted", style = "font-size: 12px;",
                      "Run the dependency analysis to score this gene against a subtype."))
      }
      tagList(
        tags$h5(class = "cpt-card-h", "CPT dimensions"),
        plotOutput(session$ns("cpt_radar_cancer"), height = "260px"),
        tags$details(
          class = "cpt-collapse-block",
          tags$summary(class = "cpt-collapse-summary",
                       icon("circle-info"), " Weight sensitivity"),
          tags$div(class = "cpt-collapse-body",
                   DTOutput(session$ns("cpt_sensitivity_table")))
        ),
        tags$div(
          style = "margin-top: 12px;",
          downloadButton(session$ns("dl_full_report"), "Download report (HTML)",
                         class = "btn-sm btn-primary", icon = icon("file-lines"))
        )
      )
    }

    # Pass the dependency cutoff through so the browse histogram can draw it.
    # cpt_scores travels with the shared data so the all-combined table can
    # carry the composite score. It is defined further down; a reactive is
    # lazy, so referring to it here is fine.
    shared_for_panel <- c(shared_data, list(
      dep_effect_min = reactive(input$effect_min %||% -0.1),
      cpt_scores = reactive(tryCatch(cancer_gene_cpt_df(), error = function(e) NULL))
    ))

    explorer_gene <- targets_panel_server(
      input, output, session, shared_for_panel,
      get_subtype = reactive(as.character(input$cancer_subtypes)),
      get_dataset = reactive(dataset_short_name(input$dataset %||% "RNAi")),
      extra_panel = drawer_detail_ui,
      on_open_gene = on_open_gene
    )

    observeEvent(input$tp_drawer_gene, {
      key <- explorer_gene()
      if (!is.null(key) && is.function(on_open_gene)) on_open_gene(key)
    }, ignoreInit = TRUE)

    # The run is a named function, not just an observer body, so the Gene tab
    # hand-off can invoke exactly the same path without faking a button click.
    # Run analysis is only meaningful when the settings differ from the ones
    # that produced what is on screen. Disabling it the rest of the time says
    # so, and stops people re-running an analysis they already have.
    # Normalised, because collapsing the "Change" disclosure makes Shiny report
    # some of these as NULL for a moment. Compared raw, that reads as a changed
    # setting and re-enabled the button when nothing had been touched.
    run_settings <- reactive({
      em <- suppressWarnings(as.numeric(input$effect_min))
      list(dataset = as.character(input$dataset %||% "RNAi"),
           subtypes = sort(as.character(input$cancer_subtypes %||% character(0))),
           effect_min = if (length(em) == 1 && is.finite(em)) em else -0.1,
           apply_pvalue = isTRUE(input$apply_pvalue),
           excl_common = isTRUE(input$excl_common),
           req_nc_sig = isTRUE(input$req_nc_sig))
    })
    ran_settings <- reactiveVal(NULL)

    observe({
      cur <- run_settings()
      has_subtype <- length(cur$subtypes) && any(nzchar(cur$subtypes))
      # Enabled when there is something to run and it is not what was run last.
      shinyjs::toggleState(
        "run_analysis",
        condition = isTRUE(has_subtype) &&
          (is.null(ran_settings()) || !identical(cur, ran_settings())))
    })

    # The precomputed path is a file read and a few filters, so the dependency
    # tables can simply appear when a subtype is chosen, exactly as Candidates
    # does. Run analysis stays for the case with no precomputed file, where the
    # gene effect matrix and a limma fit are needed.
    observe({
      if (!("dependency" %in% (input$tp_layers %||% character(0)))) return()
      st <- as.character(input$cancer_subtypes)
      if (!length(st) || !any(nzchar(st)) || length(st) != 1) return()
      if (isTRUE(isolate(deps_ready())) &&
          identical(isolate(ran_settings()), isolate(run_settings()))) return()
      idx <- isolate(precomputed_index())
      if (is.null(idx)) return()
      fpath <- cpt_effectsize_path(shared_data$data_dir,
                                   dataset_short_name(input$dataset), st[[1]])
      if (!file.exists(fpath)) return()   # needs the slow path; leave the button
      isolate(run_analysis_now())
    })

    output$run_stale_note <- renderUI({
      cur <- run_settings()
      if (is.null(ran_settings()) || identical(cur, ran_settings())) return(NULL)
      if (!length(cur$subtypes) || !any(nzchar(cur$subtypes))) return(NULL)
      stale <- c("Cancer-selective genes", "Every gene in this subtype",
                 "the dependency charts")
      if ("dependency" %in% (input$tp_layers %||% character(0)))
        stale <- c("Candidates", stale)
      tags$p(class = "cpt-run-stale",
             paste0("Showing previous settings: ",
                    paste(stale, collapse = ", "), "."))
    })
    outputOptions(output, "run_stale_note", suspendWhenHidden = FALSE)

    run_analysis_now <- function() {
      # Freeze analysis context for plots/radar; tables come from eventReactives.
      # deps_ready flips TRUE only after ge_results succeeds (see ge_results).
      #
      # IMPORTANT: Result widgets live inside conditionalPanel(deps_ready == 1), so
      # they are suspended while deps_ready is FALSE and will NOT pull ge_results().
      # We must force-evaluate the pipeline here or Run Analysis appears to do nothing.
      displayed_dataset(input$dataset)
      st <- as.character(input$cancer_subtypes)
      if (length(st) && any(nzchar(st))) {
        displayed_subtypes(st)
        displayed_subtype(paste(st, collapse = "; "))
      } else {
        displayed_subtypes(NULL)
        displayed_subtype(NULL)
        showNotification(
          "Select a cancer subtype before running analysis.",
          type = "warning", duration = 6
        )
        deps_ready(FALSE)
        return()
      }
      displayed_filters(list(
        effect_min = input$effect_min,
        apply_pvalue = isTRUE(input$apply_pvalue),
        excl_common = isTRUE(input$excl_common),
        req_nc_sig = isTRUE(input$req_nc_sig)
      ))
      ran_settings(isolate(run_settings()))
      deps_ready(FALSE)
      # Stale ligandability panels belong to the previous analysis
      lig_ready(FALSE)
      cpt_log_usage(
        "deps_run_analysis",
        list(
          dataset = input$dataset,
          subtype = if (length(st)) paste(st, collapse = "; ") else NA_character_,
          n_subtype = length(st)
        )
      )

      # compute_ge_results() loads a gene-effect matrix only if it has to fall
      # back to on-the-fly limma. The precomputed path reads one ~6 MB TSV, so
      # eagerly deserializing the 140 MB matrix here cost every run ~10 s of
      # blocked session for nothing.
      res <- tryCatch(
        compute_ge_results(),
        error = function(e) {
          msg <- conditionMessage(e)
          if (!inherits(e, "shiny.silent.error") || nzchar(msg)) {
            showNotification(
              paste0("Analysis error: ", if (nzchar(msg)) msg else class(e)[1]),
              type = "error", duration = 10
            )
          }
          NULL
        }
      )
      ge_results_val(res)
      removeNotification("cpt_deps_run")
      if (is.null(res) && !isTRUE(isolate(deps_ready()))) {
        showNotification(
          "Analysis did not complete. Check the subtype selection and that DepMap RDS files exist under data/.",
          type = "warning", duration = 10
        )
      } else if (!is.null(res)) {
        showNotification("Analysis complete.", type = "message", duration = 4)
      }
      invisible(NULL)
    }

    # Now that this tab has no cutoffs of its own, it has to say which ones it
    # is using and where they are set.
    output$lig_cutoff_note <- renderText({
      paste0("Using CR \u2265 ", input$tp_cr %||% 4,
             " and \u2264 ", input$tp_max_targets %||% 20,
             " probe targets, from the Ligandability box above.")
    })

    observeEvent(input$run_analysis, run_analysis_now(), ignoreInit = TRUE)
    # Same action, reachable from where its output is missing.
    # Target tab hand-off. Inputs set with update*Input() are not visible to the
    # server until the client echoes them back, so the preset is parked here and
    # the run fires on the echo rather than on a guessed delay.
    preset_pending <- reactiveVal(NULL)

    if (!is.null(preset)) {
      observeEvent(preset(), {
        p <- preset()
        req(p, p$subtype)
        ds_label <- if (identical(p$dataset, "RNAi")) "RNAi" else "CRISPR (23Q4)"
        updateSelectInput(session, "dataset", selected = ds_label)
        refresh_cancer_subtype_picker(ds_label, preserve_selection = FALSE)
        updatePickerInput(session, "cancer_subtypes", selected = p$subtype)
        preset_pending(p)
      }, ignoreInit = TRUE, ignoreNULL = TRUE)
    }

    observeEvent(input$cancer_subtypes, {
      p <- preset_pending()
      if (is.null(p)) return()
      if (!identical(as.character(input$cancer_subtypes)[1], as.character(p$subtype)[1])) return()
      preset_pending(NULL)
      run_analysis_now()
    }, ignoreInit = TRUE)

    # The full engagement table is loaded when something actually needs it,
    # rather than behind a button. Two things need it: opening the
    # Ligandability analysis tab, and asking for a competition ratio below 4,
    # which the precomputed index does not store — the index keeps CR >= 4
    # only, so anything lower can only come from the complete table.
    # Only a CR floor below 4 needs the full table now: everything at or above
    # the engagement threshold is served from the 3 MB engaged-only file.
    # Two different questions. The probe tables are needed whenever the
    # ligandability layer is on, whatever the cutoff; the full 10.6M-row table
    # is needed only for a CR floor under 4, since the engaged-only file starts
    # there. Answering the first with the second left every probe table empty
    # at the default cutoff.
    needs_full_table <- reactive({
      !is.null(input$tp_cr) && input$tp_cr < 4
    })

    needs_probe_tables <- reactive({
      "ligandability" %in% (input$tp_layers %||% character(0)) || isTRUE(needs_full_table())
    })

    observeEvent(needs_probe_tables(), {
      req(isTRUE(needs_probe_tables()))
      if (isTRUE(isolate(lig_ready()))) return()   # already in memory
      cpt_log_usage("deps_load_probes", list(
        trigger = if ("ligandability" %in% (input$tp_layers %||% character(0)))
          "layer" else "cr_below_4"))
      lig_ready(FALSE)
      res <- tryCatch(
        compute_probe_results(),
        error = function(e) {
          msg <- conditionMessage(e)
          if (!inherits(e, "shiny.silent.error") || nzchar(msg)) {
            showNotification(
              paste0("Probe load error: ", if (nzchar(msg)) msg else class(e)[1]),
              type = "error", duration = 10
            )
          }
          NULL
        }
      )
      probe_results_val(res)
      removeNotification("cpt_load_probes_run")
      if (is.null(res) && !isTRUE(isolate(lig_ready()))) {
        showNotification(
          paste0(
            "No probe tables loaded. Chemoproteomics data is missing or empty ",
            "(see docs/DATA_PROVENANCE.md). With the dependency layer on, run the analysis first."
          ),
          type = "warning", duration = 12
        )
      } else if (!is.null(res)) {
        n_best <- tryCatch(nrow(res$all_proteins_best_probe), error = function(e) NA_integer_)
        showNotification(
          paste0(
            "Probes loaded",
            if (!is.na(n_best)) paste0(" (", n_best, " best-probe rows)") else "",
            "."
          ),
          type = "message", duration = 5
        )
      }
    }, ignoreInit = TRUE)

    # CPT weights from sidebar (relative; scoring renormalizes over active dims).
    # Debounced: re-annotating the gene table costs ~2.3 s, and an undebounced
    # slider fires that on every intermediate value while the handle is dragged.
    cpt_user_weights <- shiny::debounce(
      reactive({
        cpt_coerce_weights(list(
          dependency_strength = input$w_dependency_strength,
          cancer_selectivity = input$w_cancer_selectivity,
          cysteine_ligandability = input$w_cysteine_ligandability,
          conservation = input$w_conservation,
          clinical_evidence = input$w_clinical_evidence,
          adme_druggability = input$w_adme_druggability
        ))
      }),
      millis = 500
    )

    apply_cpt_weight_values <- function(vals) {
      updateSliderInput(session, "w_dependency_strength", value = vals$dependency_strength)
      updateSliderInput(session, "w_cancer_selectivity", value = vals$cancer_selectivity)
      updateSliderInput(session, "w_cysteine_ligandability", value = vals$cysteine_ligandability)
      updateSliderInput(session, "w_conservation", value = vals$conservation)
      updateSliderInput(session, "w_clinical_evidence", value = vals$clinical_evidence)
      updateSliderInput(session, "w_adme_druggability", value = vals$adme_druggability)
    }

    observeEvent(input$cpt_weights_default, {
      apply_cpt_weight_values(CPT_DEFAULT_WEIGHTS)
      showNotification("CPT weights reset to defaults.", type = "message", duration = 3)
    }, ignoreInit = TRUE)

    observeEvent(input$cpt_weights_equal, {
      apply_cpt_weight_values(list(
        dependency_strength = 1,
        cancer_selectivity = 1,
        cysteine_ligandability = 1,
        conservation = 1,
        clinical_evidence = 1,
        adme_druggability = 1
      ))
      showNotification(
        paste("CPT weights set equal (1 each). This activates the provisional ADME",
              "axis, which is 0 by default \u2014 scores will shift."),
        type = "message", duration = 5
      )
    }, ignoreInit = TRUE)

    observeEvent(input$reset_deps, {
      updateSelectInput(session, "dataset", selected = "RNAi")
      refresh_cancer_subtype_picker("RNAi", preserve_selection = FALSE)
      updateNumericInput(session, "effect_min", value = -0.1)
      updateCheckboxInput(session, "apply_pvalue", value = TRUE)
      updateCheckboxInput(session, "excl_common", value = TRUE)
      updateCheckboxInput(session, "req_nc_sig", value = FALSE)
      apply_cpt_weight_values(CPT_DEFAULT_WEIGHTS)
      displayed_dataset(NULL)
      displayed_subtype(NULL)
      displayed_subtypes(NULL)
      displayed_filters(NULL)
      ran_settings(NULL)
      deps_ready(FALSE)
      lig_ready(FALSE)
      analysis_generation(0L)
    })

    output$results_dataset_badge <- renderUI({
      current_ds <- input$dataset
      shown_ds <- displayed_dataset()
      shown_st <- displayed_subtype()
      shown_f <- displayed_filters()
      current_st <- if (length(input$cancer_subtypes)) {
        paste(as.character(input$cancer_subtypes), collapse = "; ")
      } else {
        ""
      }
      req(current_ds, shown_ds)

      ds_ok <- identical(current_ds, shown_ds)
      st_ok <- identical(current_st, as.character(shown_st %||% ""))
      filt_ok <- TRUE
      if (!is.null(shown_f)) {
        filt_ok <- isTRUE(all.equal(
          list(
            effect_min = input$effect_min,
            apply_pvalue = isTRUE(input$apply_pvalue),
            excl_common = isTRUE(input$excl_common),
            req_nc_sig = isTRUE(input$req_nc_sig)
          ),
          shown_f,
          tolerance = 1e-9
        ))
      }
      st_part <- if (!is.null(shown_st) && nzchar(as.character(shown_st)[1])) {
        paste0(" | Subtype: ", as.character(shown_st)[1])
      } else {
        ""
      }

      # Nothing to say while the tables match the controls: the Showing line
      # above already names the dataset, the subtype and every threshold. This
      # output exists for the case below, where they have drifted apart.
      if (ds_ok && st_ok && filt_ok) return(NULL)

      changed_bits <- c(
        if (!ds_ok) paste0("dataset is now ", current_ds),
        if (!st_ok) paste0("subtype is now ", if (nzchar(current_st)) current_st else "(none)"),
        if (!filt_ok) "filters (effect size / p-value / essentials) changed"
      )
      app_status_banner(
        tags$b("Inputs changed. "),
        paste(changed_bits, collapse = "; "), ". ",
        "Tables/plots still show ", tags$b(shown_ds), st_part, ". Click ",
        tags$b("Run Analysis"), " to refresh.",
        type = "warning"
      )
    })

    # ---- 2.1 Precomputed effect index (for Run Analysis / TSV paths) -----------
    precomputed_index <- reactive({
      precomputed_index_dataframe(shared_data$data_dir)
    })

    ge_results_val <- reactiveVal(NULL)
    probe_results_val <- reactiveVal(NULL)
    ge_results <- reactive({ ge_results_val() })
    probe_results <- reactive({ probe_results_val() })

    compute_probe_results <- function() {
      showNotification(
        "Loading probe tables…",
        type = "message",
        duration = NULL,
        id = "load_probes_notif"
      )
      on.exit(removeNotification("load_probes_notif", session), add = TRUE)

      # Do not req(ge_results()) here: without the guard the probe load would
      # silently no-op when the analysis has not been run.
      cg_df <- NULL
      # Follows the Dependency evidence layer rather than a checkbox of its own:
      # two controls for one idea could disagree, and this one always lost.
      if (isTRUE("dependency" %in% (input$tp_layers %||% character(0)))) {
        gr <- tryCatch(ge_results(), error = function(e) NULL)
        if (is.null(gr) || is.null(gr$cancer_gene_df)) {
          showNotification(
            paste("Ligandability is restricted to the dependency results, which",
                  "have not been computed yet. Run the analysis in the Dependency",
                  "box, or switch the Dependency layer off to see every engaged",
                  "cysteine."),
            type = "warning",
            duration = 10
          )
          lig_ready(FALSE)
          return(NULL)
        }
        cg_df <- gr$cancer_gene_df
      }

      # At a CR floor of 4 or above the engaged-only file holds every row that
      # could be shown, and reads in 0.01 s against 18 s for the full table.
      pb <- if (!isTRUE(needs_full_table()) && !is.null(shared_data$protein_binding_cr4)) {
        cr4 <- tryCatch(shared_data$protein_binding_cr4(), error = function(e) NULL)
        if (is.data.frame(cr4) && nrow(cr4)) cr4 else NULL
      } else {
        NULL
      }
      if (is.null(pb)) {
        pb <- if (!is.null(shared_data$protein_binding_lookup)) {
          shared_data$protein_binding_lookup()
        } else if (!is.null(shared_data$proteinbindinglookup)) {
          shared_data$proteinbindinglookup()
        } else {
          NULL
        }
      }
      # Atlas is optional annotation; probe tables only need protein-binding RDS
      cys_atlas <- tryCatch(shared_data$cys_editing_atlas(), error = function(e) NULL)
      if (is.null(pb) || (is.data.frame(pb) && nrow(pb) == 0)) {
        showNotification(
          paste0(
            "Protein binding / chemoproteomics data is empty or missing. ",
            "Need real protein_binding_lookup_preprocessed.rds (or table-s2.xlsx + id_mapping.tsv). ",
            "See docs/DATA_PROVENANCE.md."
          ),
          type = "warning",
          duration = 12
        )
        lig_ready(FALSE)
        return(NULL)
      }

      result <- tryCatch({
        df <- pb

        # Restrict to cancer-dependency genes FIRST. The lookup is ~10.6M rows;
        # isoform/CR scans on the full table dominate the probe load on shinyapps.
        if (isTRUE("dependency" %in% (input$tp_layers %||% character(0))) &&
            !is.null(cg_df) && "gene_name" %in% colnames(df)) {
          keep_keys <- tolower(cpt_gene_match_key(cg_df$gene_name))
          indexed <- cpt_pb_subset(df, keep_keys, "rows_by_gene")
          if (!is.null(indexed)) {
            df <- indexed
          } else if ("gene_name_key" %in% colnames(df)) {
            df <- df %>% dplyr::filter(.data$gene_name_key %in% keep_keys)
          } else {
            df <- df %>%
              dplyr::filter(cpt_gene_match_key(gene_name) %in% toupper(keep_keys))
          }
        }

        # No isoform filter: none of the 10,588,541 accessions in this table
        # contains a hyphen, so the control it replaced removed nothing.

        # CR range filter
        if ("CR" %in% colnames(df)) {
          df <- df %>%
            dplyr::filter(CR >= (input$tp_cr %||% 4))
        }

        # Protein-level selectivity: distinct proteins with CR >= 4 (n_targets).
        # Do not use sum(CR >= 4); that counts cysteine sites and disagrees
        # with Protein Lookup for 548 of 998 probes.
        if ("n_targets" %in% colnames(df)) {
          df <- df %>% dplyr::mutate(targets = n_targets)
        } else if ("targets_total" %in% colnames(df)) {
          df <- df %>% dplyr::mutate(targets = targets_total)
        } else if ("probe_name" %in% colnames(df) && "CR" %in% colnames(df)) {
          df <- df %>%
            dplyr::group_by(probe_name) %>%
            dplyr::mutate(targets = dplyr::n_distinct(proteinid[CR >= 4 & !is.na(CR)])) %>%
            dplyr::ungroup()
        } else {
          # rep(), not 1L: assigning a scalar to a column of a zero-row frame
          # throws, and an empty probe frame is reachable from the filters.
          df$targets <- rep(1L, nrow(df))
        }

        # Keep probes with at least one target
        df <- df %>%
          dplyr::filter(is.na(targets) | targets >= 1)

        # Cancer-gene restriction already applied at the top of this block.
        filterable <- df

        # All probes table (already filtered to CR range / isoforms / n_targets >= 1)
        all_probes <- filterable

        # Ligandable flag: support "yes"/TRUE and fall back to all rows if column absent
        is_ligandable <- function(d) {
          if (!"ligandable" %in% colnames(d)) return(rep(TRUE, nrow(d)))
          lig <- d$ligandable
          if (is.logical(lig)) return(!is.na(lig) & lig)
          toupper(as.character(lig)) %in% c("YES", "TRUE", "1", "Y")
        }

        # Best probe per protein: fewest targets, then highest CR
        best_src <- filterable[is_ligandable(filterable), , drop = FALSE]
        if ("targets" %in% colnames(best_src)) {
          best_src <- best_src %>% dplyr::filter(targets <= (input$tp_max_targets %||% 20))
        }
        arrange_cols <- intersect(c("targets", "CR", "probe_name"), colnames(best_src))
        if (length(arrange_cols) && "gene_name" %in% colnames(best_src) && nrow(best_src)) {
          all_proteins_best_probe <- best_src %>%
            dplyr::arrange(targets, dplyr::desc(CR), probe_name) %>%
            dplyr::group_by(gene_name) %>%
            dplyr::slice(1) %>%
            dplyr::ungroup()
        } else {
          all_proteins_best_probe <- best_src
        }

        # Selective probes per protein (targets < maxTargets)
        sel_src <- filterable[is_ligandable(filterable), , drop = FALSE]
        if ("targets" %in% colnames(sel_src)) {
          sel_src <- sel_src %>% dplyr::filter(targets <= (input$tp_max_targets %||% 20))
        }
        if (nrow(sel_src) && all(c("targets", "CR", "probe_name") %in% colnames(sel_src))) {
          selective_probes_per_protein <- sel_src %>%
            dplyr::arrange(targets, dplyr::desc(CR), probe_name)
        } else {
          selective_probes_per_protein <- sel_src
        }

        # Combined dep + probe table.
        # Joins a dependency table (gene_name carries DepMap's Entrez suffix)
        # to a probe table (bare symbols), so it must go through the normalised
        # key. Joining the raw strings yields all-NA proteinid, which the filter
        # below then removes entirely -- silently emptying the combined table and
        # every panel downstream of it (combined table, selectivity scatter).
        ge_combined_final_data <- NULL
        if (isTRUE("dependency" %in% (input$tp_layers %||% character(0))) && !is.null(cg_df) &&
            "gene_name" %in% colnames(selective_probes_per_protein)) {
          lhs <- cg_df %>%
            dplyr::mutate(.cpt_gene_key = cpt_gene_match_key(gene_name))
          rhs <- selective_probes_per_protein %>%
            dplyr::mutate(.cpt_gene_key = cpt_gene_match_key(gene_name)) %>%
            dplyr::select(-gene_name)   # keep the display name from cg_df
          ge_combined_final_data <- lhs %>%
            dplyr::left_join(rhs, by = ".cpt_gene_key") %>%
            dplyr::select(-.cpt_gene_key)
          if ("proteinid" %in% colnames(ge_combined_final_data)) {
            ge_combined_final_data <- ge_combined_final_data %>%
              dplyr::filter(!is.na(proteinid) & proteinid != "")
          }
        }

        list(
          all_probes = cpt_cys_annotate_probe_table(all_probes, cys_atlas),
          all_proteins_best_probe = cpt_cys_annotate_probe_table(
            all_proteins_best_probe, cys_atlas
          ),
          selective_probes_per_protein = cpt_cys_annotate_probe_table(
            selective_probes_per_protein, cys_atlas
          ),
          ge_combined_final_data = cpt_cys_annotate_probe_table(
            ge_combined_final_data, cys_atlas
          )
        )
      }, error = function(e) {
        showNotification(paste("Probe analysis error:", e$message), type = "error")
        NULL
      })

      # Gate ligandability UI on successful load (not on button click alone)
      if (!is.null(result)) {
        lig_ready(TRUE)
      } else {
        lig_ready(FALSE)
      }
      result
    }


    # A gene-effect matrix is 140 MB on disk and ~10 s to deserialize, and it
    # blocks the whole R process while it loads. Only two things need one: the
    # on-the-fly limma fallback, and the per-gene group-comparison / radar
    # plots. So load on demand and keep it for the session, per dataset —
    # so a Run Analysis click does not re-read it from disk.
    ge_matrix_cache <- new.env(parent = emptyenv())

    load_ge_matrix <- function(ds_label = NULL) {
      ds_label <- ds_label %||% input$dataset
      if (is.null(ds_label) || !nzchar(ds_label)) return(NULL)
      fname <- shared_data$file_map[[ds_label]]
      if (is.null(fname) || !nzchar(fname)) return(NULL)
      fpath <- shared_data$data_path(fname)
      showNotification(paste0("Loading ", ds_label, " gene effect matrix..."), type = "message",
        duration = NULL, id = "load_ge")
      on.exit(removeNotification("load_ge"), add = TRUE)
      tryCatch({
        if (!file.exists(fpath)) {
          cpt_notify_error(
            "CPT-1003",
            paste0(ds_label, " gene effect matrix not found"),
            detail = paste("Expected", fpath)
          )
          return(NULL)
        }
        readRDS(fpath)
      }, error = function(e) {
        cpt_notify_error(
          "CPT-1004",
          paste0(ds_label, " gene effect matrix failed to load"),
          detail = conditionMessage(e)
        )
        NULL
      })
    }

    ge_matrix_for <- function(ds_label) {
      if (is.null(ds_label) || !nzchar(ds_label)) return(NULL)
      hit <- ge_matrix_cache[[ds_label]]
      if (!is.null(hit)) return(hit)
      mat <- load_ge_matrix(ds_label)
      if (!is.null(mat)) assign(ds_label, mat, envir = ge_matrix_cache)
      mat
    }

    ge_matrix <- reactive({ ge_matrix_for(displayed_dataset() %||% input$dataset) })

    # Prefer precomputed effect sizes when available; otherwise
    # fall back to on-the-fly lin_ass_pval + run_lm_ge + ge_analysis
    compute_ge_results <- function() {
      meta <- shared_data$cancer_model_data()
      if (is.null(meta)) {
        cpt_notify_error(
          "CPT-1001",
          "Cannot run analysis without cancer model metadata",
          detail = "cancer_model_data.rds"
        )
        deps_ready(FALSE)
        return(NULL)
      }
      req(length(input$cancer_subtypes) > 0, input$dataset)

      # --- Try precomputed path when only one subtype is selected ---
      idx <- precomputed_index()
      if (!is.null(idx) && length(input$cancer_subtypes) == 1) {
        short <- dataset_short_name(input$dataset)
        subtype <- input$cancer_subtypes[[1]]
        fpath <- cpt_effectsize_path(shared_data$data_dir, short, subtype)

        if (file.exists(fpath)) {
          df <- tryCatch(
            cpt_read_effectsizes(fpath),
            error = function(e) {
              showNotification(paste("Error reading precomputed file:", e$message),
                               type = "error")
              NULL
            }
          )
          if (!is.null(df)) {
            # Apply the same filters that construct cancer_gene_df
            cancer_df <- df %>%
              dplyr::filter(!is.na(EffectSize) & EffectSize <= input$effect_min)

            if (isTRUE(input$apply_pvalue)) {
              cancer_df <- cancer_df %>% dplyr::filter(p_value < 0.05)
            }

            if (isTRUE(input$excl_common)) {
              essential_genes <- cancer_df %>%
                dplyr::filter(Avg < -0.5) %>%
                dplyr::pull(gene_name)
              cancer_df <- cancer_df %>%
                dplyr::filter(!gene_name %in% essential_genes)
            }

            if (isTRUE(input$req_nc_sig)) {
              cancer_df <- cancer_df %>%
                dplyr::filter(!is.na(pval_vs_NonCancer) &
                                pval_vs_NonCancer < 0.05)
            }

            deps_ready(TRUE)
            analysis_generation(isolate(analysis_generation()) + 1L)
            return(list(
              all_gene_ge_df = df,
              cancer_gene_df = cancer_df
            ))
          }
        }
      }

      # --- Fallback: compute on the fly using cdsrmodels helpers ---
      # Only this path needs the full gene-effect matrix.
      progress <- shiny::Progress$new(session, min = 0, max = 100)
      on.exit(progress$close())
      progress$set(message = "Loading gene effect matrix...", value = 5)
      mat <- ge_matrix_for(input$dataset)
      if (is.null(mat)) {
        deps_ready(FALSE)
        return(NULL)
      }

      result <- tryCatch({
        validate_linear_ge_sample_size(mat, meta, input$cancer_subtypes)

        progress$set(message = "Computing linear associations...", value = 20)
        pval_res <- lin_ass_pval(mat, meta, input$cancer_subtypes)

        progress$set(message = "Running linear model...", value = 60)
        lm_res   <- run_lm_ge(mat,  meta, input$cancer_subtypes)

        progress$set(message = "Finalizing analysis...", value = 80)
        ge_analysis(
          pval_results             = pval_res,
          lm_results               = lm_res,
          apply_pvalue_filter      = isTRUE(input$apply_pvalue),
          min_effect_size          = input$effect_min,
          exclude_common_essentials = isTRUE(input$excl_common),
          filter_pval_vs_noncancer = isTRUE(input$req_nc_sig),
          effect_size_matrix       = mat,
          cancer_model_data        = meta,
          selected_subtype         = input$cancer_subtypes
        )
      }, error = function(e) {
        showNotification(paste("Analysis error:", e$message), type = "error")
        NULL
      })

      progress$set(value = 100)
      if (!is.null(result)) {
        deps_ready(TRUE)
        analysis_generation(isolate(analysis_generation()) + 1L)
      } else {
        deps_ready(FALSE)
      }
      result
    }

    # Convenience reactives for the two output data frames
    all_gene_ge_df <- reactive({
      req(ge_results())
      ge_results()$all_gene_ge_df
    })

    cancer_gene_df <- reactive({
      req(ge_results())
      ge_results()$cancer_gene_df
    })

    # ---- 2.4 Gene tables and exports ---------------------------
    # Priority columns shown by default (CPT first when present)
    # Default view: identity, rank, the dependency statistics that define the
    # set, then targetability. Everything else is behind "Show all columns" —
    # the full frame runs to about twenty columns and buries the first four.
    priority_cols <- c("gene_name", "CPT_Score", "CPT_Rank",
                       "EffectSize", "p_value",
                       "Cys_tier", "Cys_functional", "Probes", "Max_CR")

    # Helper: subset columns based on "show all" checkbox.
    #
    # A table answers the question the layers asked. Cys_tier and
    # Cys_functional are residue evidence and Probes and Max_CR are
    # ligandability, so on a dependency-only run they are columns the reader
    # never asked for and cannot act on. "Show all columns" still reveals
    # everything the frame holds.
    display_cols <- function(df, show_all) {
      if (show_all) return(df)
      on <- input$tp_layers %||% character(0)
      cols <- priority_cols
      if (!("functionality" %in% on))
        cols <- setdiff(cols, c("Cys_tier", "Cys_functional"))
      if (!("ligandability" %in% on))
        cols <- setdiff(cols, c("Probes", "Max_CR"))
      keep <- intersect(cols, colnames(df))
      df[, keep, drop = FALSE]
    }

    format_display_4dp <- function(df) {
      out <- df
      num_cols <- vapply(out, is.numeric, logical(1))
      out[num_cols] <- lapply(out[num_cols], function(x) round(x, 4))
      out
    }

    # The ADME axis is weighted 0 by default, but asking for its scores forces
    # the 10.6M-row binding table to load (~18 s deserialize + ~5 s index,
    # ~890 MB). At weight 0 it cannot move a single score, so only pay for it
    # once the axis is actually switched on.
    adme_scores_if_weighted <- function(w) {
      if (is.null(w) || !isTRUE(w$adme_druggability > 0)) return(NULL)
      if (!is.function(shared_data$adme_gene_scores)) return(NULL)
      tryCatch(shared_data$adme_gene_scores(), error = function(e) NULL)
    }

    # CPT annotation for cancer-specific table (percentiles vs full gene universe)
    # Ligandability and cysteine-tier columns come from the precomputed gene
    # index (4 MB), so the ranked table answers "is it targetable?" in the same
    # view as "is it a dependency?" without waiting for the full probe table.
    annotate_targetability <- function(df) {
      idx <- tryCatch(shared_data$gene_index(), error = function(e) NULL)
      if (is.null(idx) || is.null(df) || !nrow(df)) return(df)
      keys <- cpt_gene_match_key(df$gene_name)
      gi <- idx$genes
      m <- match(keys, gi$gene_key)
      df$Probes <- gi$n_probes[m]
      df$Max_CR <- round(gi$max_cr[m], 2)
      df$Cys_sites <- gi$n_sites[m]
      df$Cys_functional <- gi$n_functional[m]
      # Paper Tier of the gene's best ENGAGED cysteine (CR >= 4), per
      # cpt_evidence_tier(). Deliberately not the best atlas site: a functional
      # cysteine elsewhere in the protein is not evidence for the engaged one.
      df$Cys_tier <- gi$best_engaged_tier[m]
      df
    }

    cancer_gene_cpt_df <- reactive({
      cdf <- cancer_gene_df()
      adf <- all_gene_ge_df()
      req(cdf, adf)
      atlas <- tryCatch(shared_data$cys_editing_atlas(), error = function(e) NULL)
      w <- cpt_user_weights()
      annotate_targetability(
        cpt_annotate_gene_table(cdf, adf, cys_atlas = atlas,
                                weights = w, adme_data = adme_scores_if_weighted(w))
      )
    })

    # Shared data_env for CPT single-gene / compare (uses frozen analysis matrix)
    deps_data_env <- reactive({
      mat <- tryCatch(ge_matrix(), error = function(e) NULL)
      meta <- shared_data$cancer_model_data()
      atlas <- tryCatch(shared_data$cys_editing_atlas(), error = function(e) NULL)
      ds <- displayed_dataset() %||% input$dataset
      adme <- adme_scores_if_weighted(cpt_user_weights())
      env <- list(
        cancer_model_data = meta,
        cys_atlas = atlas,
        adme_gene_scores = adme,
        data_dir = shared_data$data_dir
      )
      if (!is.null(mat)) {
        if (grepl("CRISPR", ds %||% "", ignore.case = TRUE)) {
          env$crispr_matrix <- mat
        } else {
          env$rnai_matrix <- mat
        }
      }
      env
    })

    # Which gene the CPT detail outputs describe: the drawer selection first,
    # then whichever ranked table has a row selected.
    detail_gene <- reactive({
      explorer_gene() %||% input$selected_gene_cancer
    })

    selected_cpt_detail <- reactive({
      gene_id <- detail_gene()
      # Use frozen analysis context so radar/sensitivity match tables after Run
      st <- displayed_subtypes()
      ds_label <- displayed_dataset()
      req(gene_id, nzchar(gene_id), length(st) > 0, !is.null(ds_label))
      gene_clean <- sub(" \\(\\d+\\)$", "", gene_id)
      subtype <- st[[1]]
      ds <- if (grepl("CRISPR", ds_label %||% "", ignore.case = TRUE)) "CRISPR" else "RNAi"
      env <- deps_data_env()
      tryCatch(
        cpt_score_single(gene_clean, subtype, ds, env, weights = cpt_user_weights()),
        error = function(e) NULL
      )
    })

    # Table, xlsx export and volcano for one gene panel. Registered once per
    # panel instead of copied per panel; force() pins each argument so the two
    # registrations do not close over the same last value.
    wire_gene_panel <- function(key, data_reactive, volcano_reactive,
                                xlsx_stub, sort_by_score = FALSE) {
      force(key); force(data_reactive); force(volcano_reactive)
      force(xlsx_stub); force(sort_by_score)

      output[[paste0(key, "_genes_table")]] <- renderDT({
        df <- data_reactive()
        req(df)
        opts <- list(pageLength = 15, scrollX = TRUE)
        if (isTRUE(sort_by_score)) opts$order <- list(list(1, "desc"))
        display_cols(df, input[[paste0("show_all_cols_", key)]]) %>%
          format_display_4dp() %>%
          cpt_filter_levels() %>%
          datatable(options = opts, rownames = FALSE, filter = "top")
      })

      output[[paste0("dl_", key, "_genes")]] <- downloadHandler(
        filename = function() paste0(xlsx_stub, "_", Sys.Date(), ".xlsx"),
        content  = function(file) write_xlsx(data_reactive(), file)
      )

      output[[paste0("volcano_plot_", key)]] <- renderPlotly({
        ggplotly(volcano_reactive(), tooltip = c("text", "x", "y"))
      })

      output[[paste0("dl_volcano_", key)]] <- downloadHandler(
        filename = function() paste0("volcano_", key, "_", Sys.Date(), ".png"),
        content  = function(file) {
          ggplot2::ggsave(
            file,
            plot = make_volcano_static(data_reactive(), frozen_effect_min()),
            width = 10, height = 7, dpi = 300
          )
        }
      )
      invisible(NULL)
    }

    cpt_radar_gg_reactive <- reactive({
      det <- selected_cpt_detail()
      shiny::validate(need(!is.null(det), "Select a cancer-selective gene to score."))
      vec <- cpt_dimension_vector(det)
      shiny::validate(need(length(vec) > 0, "No CPT dimensions available."))
      cpt_dimension_radar_gg(
        vec,
        title = paste0(det$gene, " · CPT ",
                       if (is.na(det$cpt_score)) "NA" else round(det$cpt_score, 1))
      )
    })

    output$cpt_radar_cancer <- renderPlot({
      cpt_radar_gg_reactive()
    })

    output$dl_cpt_radar_cancer <- downloadHandler(
      filename = function() {
        g <- input$selected_gene_cancer %||% "gene"
        paste0("cpt_radar_", sub(" \\(\\d+\\)$", "", g), "_", Sys.Date(), ".png")
      },
      content = function(file) {
        ggplot2::ggsave(file, plot = cpt_radar_gg_reactive(),
                        width = 7, height = 6, dpi = 300, bg = "white")
      }
    )

    output$cpt_sensitivity_table <- renderDT({
      det <- selected_cpt_detail()
      shiny::validate(need(!is.null(det), "Select a gene under Group Comparison."))
      vec <- cpt_dimension_vector(det)
      sens <- cpt_weight_sensitivity(vec, weights = cpt_user_weights())
      shiny::validate(need(nrow(sens) > 0, "No dimensions available for sensitivity."))
      datatable(sens, options = list(pageLength = 12, scrollX = TRUE),
                rownames = FALSE)
    })

    # Gene pickers: refresh choices when analysis data is ready; auto-select top CPT
    # only when analysis_generation bumps (successful Run), so manual picks stick.
    pick_top_cancer_gene <- function(cancer_df) {
      if (is.null(cancer_df) || !nrow(cancer_df)) return(NULL)
      if ("CPT_Score" %in% names(cancer_df) && any(!is.na(cancer_df$CPT_Score))) {
        es <- cancer_df$EffectSize
        if (is.null(es)) es <- rep(0, nrow(cancer_df))
        es[is.na(es)] <- 0
        ord <- order(-cancer_df$CPT_Score, es)
        return(cancer_df$gene_name[ord[1]])
      }
      if ("EffectSize" %in% names(cancer_df)) {
        return(cancer_df$gene_name[order(cancer_df$EffectSize)[1]])
      }
      cancer_df$gene_name[1]
    }

    observeEvent(analysis_generation(), {
      gen <- analysis_generation()
      if (is.null(gen) || gen < 1L) return()
      all_df <- tryCatch(all_gene_ge_df(), error = function(e) NULL)
      cancer_df <- tryCatch(cancer_gene_cpt_df(), error = function(e) NULL)
      if (is.null(cancer_df) || !nrow(cancer_df)) {
        cancer_df <- tryCatch(cancer_gene_df(), error = function(e) NULL)
      }
      if (is.null(all_df) || !nrow(all_df)) return()
      cancer_choices <- if (!is.null(cancer_df) && nrow(cancer_df)) {
        sort(unique(cancer_df$gene_name))
      } else {
        character(0)
      }
      top_cancer <- pick_top_cancer_gene(cancer_df)
      updateSelectizeInput(session, "selected_gene_cancer",
        choices = cancer_choices, selected = top_cancer, server = TRUE)
    }, ignoreInit = TRUE)

    # ---- 2.5 Volcano plot --------------------------------------
    # Cancer-specific tab uses filtered genes; All Genes uses full universe.
    frozen_effect_min <- function() {
      f <- displayed_filters()
      if (!is.null(f) && !is.null(f$effect_min)) f$effect_min else input$effect_min
    }

    make_volcano_gg <- function(df, effect_min, subtitle = NULL, max_points = NULL) {
      req(df, "EffectSize" %in% colnames(df))

      if ("p_value" %in% colnames(df)) {
        pvec <- df$p_value
      } else if ("adj.P.Val" %in% colnames(df)) {
        pvec <- df$adj.P.Val
      } else {
        pvec <- rep(NA_real_, nrow(df))
      }

      pnum <- suppressWarnings(as.numeric(pvec))
      pnum[!is.finite(pnum) | pnum <= 0] <- NA_real_

      df$neg_log10_p <- -log10(pnum)
      df$significant <- !is.na(df$neg_log10_p) & df$EffectSize <= effect_min

      # ggplotly of 18k points stalls shinyapps; keep all hits, sample the rest.
      if (!is.null(max_points) && nrow(df) > max_points) {
        sig_idx <- which(df$significant)
        rest_idx <- which(!df$significant)
        n_keep <- max(0L, as.integer(max_points) - length(sig_idx))
        if (length(rest_idx) > n_keep) {
          rest_idx <- rest_idx[sample.int(length(rest_idx), n_keep)]
        }
        df <- df[sort(c(sig_idx, rest_idx)), , drop = FALSE]
      }

      p <- ggplot(df, aes(x = EffectSize, y = neg_log10_p,
                     colour = significant, text = gene_name)) +
        geom_point(alpha = 0.8, size = 2.6) +
        geom_vline(xintercept = effect_min,
                   linetype = "dashed", colour = CPT_PAL$alert) +
        scale_colour_manual(values = c("FALSE" = CPT_PAL$muted,
                                       "TRUE"  = CPT_PAL$ink)) +
        labs(x = "Effect Size", y = "-log10(p-value)", subtitle = subtitle) +
        cpt_theme() +
        theme(legend.position = "none")
      p
    }

    volcano_gg_cancer <- reactive({
      # Filtered cancer-hit set from last Run (not the full-universe volcano)
      make_volcano_gg(
        cancer_gene_df(), frozen_effect_min(),
        subtitle = "Cancer-selective gene set from last Run (see All Genes for full volcano)"
      )
    })


    # Separate download handlers so each tab can save correctly
    make_volcano_static <- function(df, effect_min) {
        req(df, "EffectSize" %in% colnames(df))

        if ("p_value" %in% colnames(df)) {
          pvec <- df$p_value
        } else if ("adj.P.Val" %in% colnames(df)) {
          pvec <- df$adj.P.Val
        } else {
          pvec <- rep(NA_real_, nrow(df))
        }

        pnum <- suppressWarnings(as.numeric(pvec))
        pnum[!is.finite(pnum) | pnum <= 0] <- NA_real_

        df$neg_log10_p <- -log10(pnum)
        df$significant <- !is.na(df$neg_log10_p) & df$EffectSize <= effect_min

        p_static <- ggplot(df, aes(x = EffectSize, y = neg_log10_p,
                                   colour = significant)) +
          geom_point(alpha = 0.8, size = 2.6) +
          geom_vline(xintercept = effect_min,
            linetype = "dashed", colour = CPT_PAL$alert) +
          scale_colour_manual(values = c("FALSE" = CPT_PAL$muted,
                                         "TRUE"  = CPT_PAL$ink)) +
          labs(x      = "Effect Size",
               y      = "-log10(p-value)") +
          cpt_theme() +
          theme(legend.position = "none")
      }

    # Both gene panels are registered here, after make_volcano_static() and the
    # volcano reactives they depend on exist.
    wire_gene_panel("cancer", cancer_gene_cpt_df, volcano_gg_cancer,
                    "cancer_specific_genes", sort_by_score = TRUE)

    # Group comparison moved to the Target tab, where box plots are drawn from
    # precomputed per-subtype summaries. Per-cell-line points are therefore no
    # longer plotted anywhere; restoring them needs either the matrix at runtime
    # or a per-line precompute.

    # ---- Presenting a joined frame ------------------------------
    # These tables are built by joining files that disagree about column names
    # and carry intermediates no reader wants. Take the first name that exists,
    # round it, drop columns that are empty for every row, and apply the row
    # grain, which is a property of the view rather than of one table.
    # A name that is not there stands in as NA once per row. Returning a bare
    # NA works while the frame has rows, because data.frame() recycles length
    # one, and throws the moment a filter empties it.
    col_first <- function(d, ..., fill = NA_character_) {
      nm <- intersect(c(...), names(d))
      if (!length(nm)) rep(fill, nrow(d)) else d[[nm[1]]]
    }
    col_num <- function(d, digits, ...) {
      x <- col_first(d, ...)
      if (is.null(x)) NULL else round(suppressWarnings(as.numeric(x)), digits)
    }
    col_yn <- function(d, ...) {
      x <- col_first(d, ...)
      # NA stays NA. Reading a missing flag as "no" would assert something the
      # data never said, and would stop drop_empty_cols() removing the column.
      ifelse(is.na(x), NA_character_, ifelse(cpt_cys_true(x), "yes", "no"))
    }
    drop_empty_cols <- function(out) {
      # all(is.na(x)) is TRUE for a zero-length column, so on an empty frame
      # this would drop every column and leave later code with nothing to
      # sort on. There is nothing to prune when there are no rows anyway.
      if (!nrow(out)) return(out)
      out[, !vapply(out, function(x) all(is.na(x)), logical(1)), drop = FALSE]
    }
    # Sort on the keys that are actually present. drop_empty_cols() removes a
    # column that is NA for every row, so naming one here unconditionally
    # passed NULL to order(), which fails with "argument is not a vector"
    # exactly when a filter leaves nothing to show.
    order_present <- function(out, ...) {
      specs <- list(...)
      keys <- list()
      for (sp in specs) {
        nm <- sp$col
        if (!nm %in% names(out)) next
        v <- out[[nm]]
        if (isTRUE(sp$desc)) {
          v <- suppressWarnings(as.numeric(v))
          v <- -ifelse(is.na(v), -Inf, v)
        }
        keys <- c(keys, list(v))
      }
      if (!length(keys) || !nrow(out)) return(out)
      out[do.call(order, c(keys, list(na.last = TRUE))), , drop = FALSE]
    }

    apply_grain <- function(out) {
      key <- switch(input$tp_grain %||% "probe",
                    gene = out$Gene,
                    cysteine = paste(out$Protein, out$Site),
                    NULL)
      if (is.null(key)) out else out[!duplicated(key), , drop = FALSE]
    }
    site_of <- function(cid) {
      cid <- as.character(cid)
      ifelse(is.na(cid), NA_character_, sub("^.*_C", "C", cid))
    }

    # Ligandability-layer table, so it shows ligandability columns. The cys_*
    # fields the probe annotator attaches are residue evidence and ride along
    # only when that layer is on; the internal join keys never show at all.
    best_probe_display <- reactive({
      req(probe_results())
      d <- probe_results()$all_proteins_best_probe
      req(is.data.frame(d), nrow(d) > 0)
      col <- function(nm, fill = NA_character_) {
        if (nm %in% names(d)) d[[nm]] else rep(fill, nrow(d))
      }
      out <- data.frame(
        Gene = col("gene_name"),
        Protein = col("proteinid"),
        Site = col("cysteineid"),
        `Best probe` = col("probe_name"),
        CR = round(suppressWarnings(as.numeric(col("CR"))), 2),
        `Probe targets` = if ("targets" %in% names(d)) d$targets else col("n_targets"),
        Dataset = col("Dataset"),
        `Cell line` = col("Cell_Line"),
        check.names = FALSE, stringsAsFactors = FALSE
      )
      # No atlas columns: Best probe per protein belongs to the ligandability
      # layer, and a per-layer view shows that layer only.
      order_present(out, list(col = "Probe targets"), list(col = "CR", desc = TRUE))
    })

    output$best_probe_table <- renderDT({
      datatable(best_probe_display(),
        rownames = FALSE, filter = "top", selection = "none",
        options = list(pageLength = 15, scrollX = TRUE))
    })

    # Dependency + ligandability. Identity first, then the dependency evidence,
    # then the chemistry, then the residue annotation when that layer is on.
    # Previously this rendered the joined frame as-is: every intermediate
    # column of the effect-size file at full float precision, and one row per
    # probe record whatever the row grain said.
    combined_display <- reactive({
      req(probe_results())
      d <- probe_results()$ge_combined_final_data
      req(is.data.frame(d), nrow(d) > 0)
      out <- data.frame(
        Gene = col_first(d, "gene_name"),
        Protein = col_first(d, "proteinid"),
        Site = site_of(col_first(d, "cysteineid")),
        `Effect size` = col_num(d, 4, "EffectSize", "effect_size"),
        `Mean in subtype` = col_num(d, 3, "Cancer_Avg", "cancer_avg"),
        `Mean across all lines` = col_num(d, 3, "Avg", "avg"),
        p = col_num(d, 4, "p.value", "p_value"),
        q = col_num(d, 4, "q.value", "q_value"),
        `-log10 p` = col_num(d, 2, "neg_log10_p_value", "neg_log10_p"),
        Probe = col_first(d, "probe_name"),
        SMILES = col_first(d, "SMILES"),
        CR = col_num(d, 2, "CR"),
        `Probe targets` = col_first(d, "targets", "n_targets"),
        Dataset = col_first(d, "Dataset"),
        `Cell line` = col_first(d, "Cell_Line"),
        check.names = FALSE, stringsAsFactors = FALSE
      )
      if ("functionality" %in% (input$tp_layers %||% character(0))) {
        out$`In atlas` <- col_yn(d, "cys_site_in_atlas")
        out$Functional <- col_yn(d, "cys_functional")
        out$`Atlas ligandable` <- col_yn(d, "cys_ligandable")
      }
      out <- drop_empty_cols(out)
      out <- order_present(out, list(col = "Probe targets"),
                                list(col = "CR", desc = TRUE))
      apply_grain(out)
    })

    output$combined_table <- renderDT({
      datatable(combined_display(),
        rownames = FALSE, filter = "top", selection = "none",
        options = list(pageLength = 15, scrollX = TRUE))
    })

    cys_probe_targets <- reactive({
      req(probe_results())
      df <- probe_results()$ge_combined_final_data
      if (is.null(df)) {
        df <- probe_results()$all_probes
      }
      req(df)
      df <- df %>%
        dplyr::filter(
          tolower(as.character(.data$ligandable)) == "yes"
        )
      atlas <- shared_data$cys_editing_atlas()
      df <- cpt_annotate_engaged_tiers(df, atlas)
      linked <- cpt_cys_link_functional_sites(df, atlas)
      if (!is.null(linked) && nrow(linked) > 0) {
        key_cols <- intersect(c("probe_name", "cysteineid", "gene_name"), names(df))
        if (length(key_cols)) {
          linked_key <- do.call(paste, c(linked[key_cols], sep = "\t"))
          df_key <- do.call(paste, c(df[key_cols], sep = "\t"))
          extra <- df[!df_key %in% linked_key, , drop = FALSE]
          if (nrow(extra) > 0) {
            extra$cys_probe_site_match <- extra$cys_site_in_atlas %in% TRUE &
              extra$cys_functional %in% TRUE
            extra$cys_relationship <- ifelse(
              extra$evidence_tier %in% c(1L, 2L),
              "Probe binds this functional site",
              "Engaged site; no functional Cys_editing match in this gene"
            )
            linked <- dplyr::bind_rows(linked, extra)
          }
        }
        df <- linked
      }
      if (isTRUE(input$cys_require_atlas_ligandable)) {
        df <- df %>% dplyr::filter(.data$cys_functional_ligandable %in% TRUE |
                                     .data$evidence_tier %in% 1L)
      }
      df <- df %>%
        dplyr::arrange(
          .data$evidence_tier,
          dplyr::desc(.data$cys_probe_site_match),
          dplyr::desc(.data$cys_functional_ligandable),
          .data$targets,
          dplyr::desc(.data$CR)
        )
      if ("EffectSize" %in% names(df)) {
        df <- df %>% dplyr::arrange(.data$evidence_tier, .data$EffectSize)
      }
      df
    })

    output$cys_probe_status <- renderUI({
      req(probe_results())
      combined <- probe_results()$ge_combined_final_data
      if (is.null(combined)) combined <- probe_results()$all_probes
      req(combined)
      functional <- cpt_cys_link_functional_sites(
        combined %>% dplyr::filter(tolower(as.character(.data$ligandable)) == "yes"),
        shared_data$cys_editing_atlas()
      )
      direct <- functional %>% dplyr::filter(.data$cys_probe_site_match)
      confirmed <- functional %>%
        dplyr::filter(.data$cys_functional_ligandable)
      app_status_banner(
        tags$b(dplyr::n_distinct(functional$gene_name)),
        " ligandable dependency gene(s) contain a functional cysteine; ",
        tags$b(dplyr::n_distinct(direct$gene_name)),
        " have a probe at that exact functional site; ",
        tags$b(dplyr::n_distinct(confirmed$gene_name)),
        " also have ligandability evidence in the Cys_editing atlas. Evidence tiers: 1 = exact functional + atlas ligandable; 2 = exact functional; 3 = tested non-functional; 4 = untested (absent from atlas). The relationship column still flags same-gene functional sites that are not the engaged residue.",
        type = "info"
      )
    })

    # Ligandability + residue evidence. Identity, then what is known about the
    # residue, then the chemistry that engages it. Previously this listed a
    # priority block and then appended every remaining column of the join.
    cys_probe_display <- reactive({
      d <- cys_probe_targets()
      req(is.data.frame(d), nrow(d) > 0)
      out <- data.frame(
        Gene = col_first(d, "gene_name"),
        Protein = col_first(d, "proteinid"),
        Site = site_of(col_first(d, "cysteineid")),
        Tier = col_first(d, "evidence_tier"),
        `Tier meaning` = col_first(d, "evidence_tier_label"),
        Functional = col_yn(d, "cys_functional"),
        `Atlas ligandable` = col_yn(d, "cys_ligandable"),
        `Functional + ligandable` = col_yn(d, "cys_functional_ligandable"),
        `Editor support` = col_first(d, "cys_editor_support"),
        `Study context` = col_first(d, "cys_study_context"),
        Relationship = col_first(d, "cys_relationship"),
        Probe = col_first(d, "probe_name"),
        SMILES = col_first(d, "SMILES"),
        CR = col_num(d, 2, "CR"),
        `Probe targets` = col_first(d, "targets", "n_targets"),
        check.names = FALSE, stringsAsFactors = FALSE
      )
      # Dependency numbers only when that layer is on.
      if ("dependency" %in% (input$tp_layers %||% character(0))) {
        out$`Effect size` <- col_num(d, 4, "EffectSize", "effect_size")
        out$p <- col_num(d, 4, "p_value", "p.value")
      }
      out <- drop_empty_cols(out)
      out <- order_present(out, list(col = "Tier"), list(col = "CR", desc = TRUE))
      apply_grain(out)
    })

    output$cys_probe_table <- renderDT({
      df <- cys_probe_display()
      shiny::validate(need(
        nrow(df) > 0,
        "No dependency-linked probes meet the current functional-cysteine filters."
      ))
      datatable(df, rownames = FALSE, filter = "top", selection = "none",
                options = list(pageLength = 15, scrollX = TRUE))
    })

    output$dl_best_probe <- downloadHandler(
      filename = function() paste0("best_probe_per_protein_", Sys.Date(), ".xlsx"),
      content  = function(file) write_xlsx(best_probe_display(), file)
    )

    output$dl_combined <- downloadHandler(
      filename = function() paste0("combined_dep_probe_", Sys.Date(), ".xlsx"),
      content  = function(file) write_xlsx(combined_display(), file)
    )

    output$dl_cys_probe_targets <- downloadHandler(
      filename = function() paste0("dependency_functional_cysteine_targets_", Sys.Date(), ".xlsx"),
      content = function(file) write_xlsx(cys_probe_display(), file)
    )

    # Probe selectivity scatter (targets vs CR)
    probe_scatter_gg <- reactive({
      req(probe_results())
      df <- probe_results()$ge_combined_final_data
      req(!is.null(df), nrow(df) > 0)

      # Use all selective probes from combined table with CR >= 4
      df <- df %>%
        dplyr::filter(!is.na(CR), CR >= 4,
                      !is.na(targets), targets <= (input$tp_max_targets %||% 20))

      req(nrow(df) > 0)

      ggplot(df, aes(x = targets, y = CR,
                     text = paste0(probe_name, " — ", gene_name))) +
        geom_point(alpha = 0.8, size = 2.6, colour = CPT_PAL$ink) +
        labs(x      = "Number of Targets",
             y      = "Competition Ratio (CR)") +
        cpt_theme() +
        theme(legend.position = "none")
    })

    output$probe_scatter_best <- renderPlotly({
      ggplotly(probe_scatter_gg(), tooltip = "text")
    })

    output$probe_scatter_combined <- renderPlotly({
      ggplotly(probe_scatter_gg(), tooltip = "text")
    })

    output$dl_probe_scatter_best <- downloadHandler(
      filename = function() paste0("selective_probes_best_", Sys.Date(), ".png"),
      content  = function(file) {
        ggplot2::ggsave(file, plot = probe_scatter_gg(),
                        width = 8, height = 6, dpi = 300)
      }
    )

    output$dl_probe_scatter_combined <- downloadHandler(
      filename = function() paste0("selective_probes_combined_", Sys.Date(), ".png"),
      content  = function(file) {
        ggplot2::ggsave(file, plot = probe_scatter_gg(),
                        width = 8, height = 6, dpi = 300)
      }
    )

    # Ligandable dependencies volcano (subset of all_gene_ge_df that has a CR ≥ 4 probe).
    lig_volcano_gg <- reactive({
      req(all_gene_ge_df())
      effect_cut <- frozen_effect_min()

      req(probe_results())
      # Probe gene names are bare symbols; all_gene_ge_df() carries the Entrez suffix.
      lig_keys <- cpt_gene_match_key(probe_results()$all_probes$gene_name)
      df <- all_gene_ge_df() %>%
        dplyr::filter(cpt_gene_match_key(gene_name) %in% lig_keys)

      req(nrow(df) > 0)

      # Build p-values as for main volcano
      if ("p_value" %in% colnames(df)) {
        pvec <- df$p_value
      } else if ("adj.P.Val" %in% colnames(df)) {
        pvec <- df$adj.P.Val
      } else {
        pvec <- rep(NA_real_, nrow(df))
      }

      pnum <- suppressWarnings(as.numeric(pvec))
      pnum[!is.finite(pnum) | pnum <= 0] <- NA_real_

      df$neg_log10_p <- -log10(pnum)
      df$significant <- !is.na(df$neg_log10_p) & df$EffectSize <= effect_cut

      ggplot(df, aes(x = EffectSize, y = neg_log10_p,
                     colour = significant, text = gene_name)) +
        geom_point(alpha = 0.8, size = 2.6) +
        geom_vline(xintercept = effect_cut,
          linetype = "dashed", colour = CPT_PAL$alert) +
        scale_colour_manual(values = c("FALSE" = CPT_PAL$muted,
                                       "TRUE"  = CPT_PAL$ink)) +
        labs(x      = "Effect Size",
             y      = "-log10(p-value)") +
        cpt_theme() +
        theme(legend.position = "none")
    })

    output$lig_volcano_best <- renderPlotly({
      ggplotly(lig_volcano_gg(), tooltip = c("text", "x", "y"))
    })

    output$dl_lig_volcano_best <- downloadHandler(
      filename = function() paste0("ligandable_volcano_best_", Sys.Date(), ".png"),
      content  = function(file) {
        ggplot2::ggsave(file, plot = lig_volcano_gg(),
                        width = 10, height = 7, dpi = 300)
      }
    )

    # ---- 2.8 Full report download --------------------------------
    # The HTML button sits above both Dependencies subtabs. Each subtab has its
    # own Group Comparison picker; Cancer-Selective Genes is auto-filled with the
    # top CPT gene, so reading only that picker ignores All Genes selections.
    resolve_report_gene <- function() {
      clean_gene <- function(x) {
        if (is.null(x) || !length(x)) return(NULL)
        x <- as.character(x)[[1L]]
        if (!nzchar(x)) return(NULL)
        sub(" \\(\\d+\\)$", "", x)
      }
      from_drawer <- clean_gene(tryCatch(explorer_gene(), error = function(e) NULL))
      if (!is.null(from_drawer)) return(from_drawer)
      picked <- clean_gene(input$selected_gene_cancer)
      if (!is.null(picked)) return(picked)
      cg <- tryCatch(cancer_gene_cpt_df(), error = function(e) NULL)
      if (is.null(cg) || !nrow(cg)) cg <- tryCatch(cancer_gene_df(), error = function(e) NULL)
      if (is.null(cg) || !nrow(cg)) return(NULL)
      g <- pick_top_cancer_gene(cg)
      clean_gene(g)
    }

    output$dl_full_report <- downloadHandler(
      filename = function() {
        g <- tryCatch(resolve_report_gene(), error = function(e) NULL)
        gene_part <- if (!is.null(g) && nzchar(g)) sanitize_subtype(g) else "gene"
        st_part <- if (!is.null(displayed_subtype()) && nzchar(displayed_subtype())) {
          sanitize_subtype(displayed_subtype())
        } else {
          "analysis"
        }
        paste0("CanProTarget_", gene_part, "_", st_part, "_", Sys.Date(), ".html")
      },
      content = function(file) {
        req(displayed_dataset(), displayed_subtype())
        ds <- displayed_dataset()
        st <- displayed_subtype()

        dataset_api <- if (exists("cpt_normalize_dataset", mode = "function")) {
          tryCatch(cpt_normalize_dataset(ds), error = function(e) {
            if (grepl("CRISPR", ds, ignore.case = TRUE)) "CRISPR" else "RNAi"
          })
        } else if (grepl("CRISPR", ds, ignore.case = TRUE)) {
          "CRISPR"
        } else {
          "RNAi"
        }

        withProgress(message = "Generating report...", value = 0.3, {
          report_gene <- resolve_report_gene()
          if (is.null(report_gene) || !nzchar(report_gene)) {
            showNotification("No genes found in current analysis.", type = "warning")
            # Write a minimal HTML so the browser download is not empty/broken
            writeLines(
              "<!DOCTYPE html><html><body><p>No genes found in current analysis.</p></body></html>",
              file
            )
            return(invisible(NULL))
          }

          # Strip Entrez ID if present: "GENE (12345)" -> "GENE"
          report_gene <- sub(" \\(\\d+\\)$", "", report_gene)

          # Prefer primary frozen subtype for report params
          report_subtype <- {
            st_vec <- displayed_subtypes()
            if (!is.null(st_vec) && length(st_vec)) st_vec[[1]] else st
          }

          incProgress(0.4, detail = paste("Rendering report for", report_gene))

          out_path <- tryCatch(
            cpt_report_gene_dependency(
              gene = report_gene,
              subtype = report_subtype,
              dataset = dataset_api,
              project_root = getwd()
            ),
            error = function(e) {
              showNotification(paste("Report error:", e$message), type = "error")
              NULL
            }
          )

          incProgress(0.3)
          if (!is.null(out_path) && file.exists(out_path)) {
            file.copy(out_path, file, overwrite = TRUE)
            unlink(out_path)
          } else {
            safe_gene <- gsub("[<>&\"]", "", report_gene)
            writeLines(
              paste0(
                "<!DOCTYPE html><html><body><p>Report generation failed for ",
                safe_gene, ".</p></body></html>"
              ),
              file
            )
          }
        })
      }
    )

    # Expose the current cancer-specific dependency symbols to other modules.
    # Before the first successful analysis this is intentionally empty.
    list(
      cancer_genes = reactive({
        # Keyed on a completed analysis, not on the button. The precomputed
        # path runs on its own without the button ever being pressed, and
        # guarding on the click count reported no dependency genes while 500
        # of them were on screen -- which emptied the atlas table, since the
        # dependency layer restricts it to exactly this list.
        if (analysis_generation() < 1L) {
          return(character(0))
        }
        result <- ge_results()
        if (is.null(result) || is.null(result$cancer_gene_df) ||
            !"gene_name" %in% colnames(result$cancer_gene_df)) {
          return(character(0))
        }
        unique(as.character(result$cancer_gene_df$gene_name))
      }),
      # The Target tab's per-cell-line table needs every line, not just the
      # subtype's, and the index deliberately stores only the subtype's points.
      # This hands over the same session-cached matrix Discover uses, so the
      # cost is paid once per dataset rather than once per tab.
      ge_matrix_for = ge_matrix_for,
      # The cysteine atlas browse view has no filter box of its own; it reads
      # the Cysteine function box and the action bar through these, so the two
      # cannot drift apart.
      residue_filters = list(
        genes = reactive({
          raw <- input$tp_genes
          if (is.null(raw) || !nzchar(trimws(raw))) return(NULL)
          g <- cpt_gene_match_key(trimws(unlist(strsplit(raw, "[,;\n\t ]+"))))
          g <- g[nzchar(g)]
          if (!length(g)) NULL else unique(g)
        }),
        contexts  = reactive(input$tp_context),
        editor    = reactive(input$tp_editor %||% "Any editor"),
        tiers     = reactive(input$tp_tiers %||% c("1", "2", "3", "4")),
        clinvar   = reactive(isTRUE(input$tp_clinvar)),
        atlas_lig = reactive(isTRUE(input$tp_atlas_lig)),
        # Follows the Ligandability layer rather than a checkbox of its own, so
        # the atlas table shows the same population as every other view.
        engaged   = reactive("ligandability" %in% (input$tp_layers %||% character(0))),
        # Restricting the atlas to the dependency results is what the Dependency
        # layer means, so it follows the layer rather than a checkbox that could
        # disagree with it.
        dep_only  = reactive("dependency" %in% (input$tp_layers %||% character(0)))
      ),
      # The Target tab's dependency table honours the same cutoff, so the two
      # tabs cannot disagree about what counts as a dependency.
      effect_min = reactive(input$effect_min %||% -0.1),
      # Full Discover context, so the Target tab reads the same dataset and
      # subtype rather than keeping a second set of controls that can disagree.
      dataset = reactive(dataset_short_name(input$dataset %||% "RNAi")),
      subtype = reactive({
        st <- as.character(input$cancer_subtypes)
        if (length(st) && nzchar(st[1])) st[1] else NULL
      }),
      # Setter so the Target tab's dataset dropdown drives the same input
      # rather than being a second, divergent control.
      set_dataset = function(value) {
        if (!identical(value, isolate(input$dataset))) {
          updateSelectInput(session, "dataset", selected = value)
        }
      },
      context_label = reactive({
        st <- as.character(input$cancer_subtypes)
        paste0(
          input$dataset %||% "RNAi",
          " \u00b7 ", if (length(st) && nzchar(st[1])) st[1] else "no subtype",
          " \u00b7 effect \u2264 ", input$effect_min %||% -0.1,
          if (isTRUE(input$apply_pvalue)) " \u00b7 p<0.05" else "",
          if (isTRUE(input$excl_common)) " \u00b7 no common essentials" else ""
        )
      })
    )

  }) # end moduleServer
} # end dependencies_server
