# ============================================================
# Script:   app.R
# Purpose:  CanProTarget — entry point for the Shiny app.
#           - Loads all packages
#           - Sets data directory (change data_dir for deployment)
#           - Sources all R/ module scripts
#           - Defines shared data reactives (loaded once, passed to all modules)
#           - Wires UI and server via Shiny modules
#
# Usage:
#   Local:        Open in RStudio → Run App  (data_dir = "data")
#   shinyapps.io: Publish button in RStudio; include data/ folder in bundle
#
# Static images (logo, graphical abstract) must live in www/ — Shiny serves that
# folder at the app URL root (independent of getwd() when using runApp("path")).
#
# Required packages (install once):
#   install.packages(c(
#     "shiny", "shinydashboard", "dplyr", "tidyr", "ggplot2",
#     "plotly", "DT", "tibble", "writexl", "magrittr",
#     "shinyjs", "shinycssloaders", "shinyWidgets", "colourpicker", "ggrepel",
#     "gtools", "jsonlite", "readr", "rmarkdown", "yaml", "remotes"
#   ))
#   remotes::install_github("broadinstitute/cdsrmodels")
#
#   Versions the release was built against: docs/R_ENVIRONMENT.md
# ============================================================
# Sections:
#   1. Package loading
#   2. Data directory configuration
#   3. File map (dataset name → filename)
#   4. Source helpers (functions.R, ggplot theme)
#   5. Source modules
#   6. UI definition
#   7. Server definition (shared data + module wiring)
# ============================================================

# ---- 1. Package loading ----------------------------------------
library(shiny)
library(shinydashboard)
library(DT)
library(dplyr)
library(magrittr)
library(tibble)
library(cdsrmodels)
library(ggplot2)
library(plotly)
library(colourpicker)
library(shinyjs)
library(tidyr)
library(writexl)
library(shinyWidgets)
library(ggrepel)
library(gtools)
library(shinycssloaders)

# Resolve namespace conflict: plotly imports jsonlite which masks shiny::validate.
# Force shiny's validate to win in the global search path.
validate <- shiny::validate

# ---- 2. Configuration + helpers --------------------------------
source("R/app_config.R")
source("R/app_helpers.R")
source("R/ui_kit.R")
source("R/error_codes.R")
source("R/usage_logger.R")
source("R/theme_canprotarget.R")
source("R/plot_palette.R")           # shared chart palette, marker sizes, plotly defaults
source("R/canprotarget_score.R")
source("R/pancancer_profile.R")
source("R/functions.R")
source("R/api_functions.R")
source("R/cys_editing_functions.R")
source("R/report_generator.R")

data_dir <- app_config$data_dir
file_map <- app_config$file_map
cancer_subtype_list_files <- app_config$cancer_subtype_list_files
data_path <- function(...) app_data_path(data_dir, ...)

# ---- 3. Source modules -----------------------------------------
# Load all Shiny module files from R/ directory
source("R/gene_module.R")                 # Target tab (target-first entry point)
source("R/targets_panel.R")               # Discover tab - layered explorer
source("R/dependencies_module.R")         # Discover tab (subtype-first workflow)
source("R/swissadme_module.R")            # Chemistry tab (SwissADME)
source("R/cys_editing_module.R")          # Cysteine Atlas tab

# Sidebar width (px) — configured in R/app_config.R
.sidebar_width <- app_config$sidebar_width

# ---- 4. UI definition ------------------------------------------
ui <- dashboardPage(
  title = "CanProTarget",
  dashboardHeader(
    title = span(
      style = "display: inline-flex; align-items: center; column-gap: 10px;",
      tags$img(
        src = "cpt_logo.png",
        height = "36px",
        alt = "CanProTarget",
        style = "vertical-align: middle;"
      ),
      "CanProTarget"
    ),
    titleWidth = .sidebar_width,
    tags$li(
      class = "dropdown",
      actionButton("share_view", label = "Share This View",
        icon = icon("link"),
        style = "background:transparent;border:none;color:#fff;margin-top:10px;")
    )
  ),
  dashboardSidebar(
    width = .sidebar_width,
    # Tabs are ordered by the question the user arrives with:
    #   Gene      — "what do we know about X?"        (no run required)
    #   Discover  — "what should we target in Y?"     (subtype workflow)
    #   Cysteine Atlas — "which sites are functional?"
    #   Chemistry — "what can we put on it?"          (probes + ADME)
    # tabName values are unchanged so existing bookmark URLs still resolve.
    sidebarMenu(
      id = "tabs",
      menuItem("Home", tabName = "home_tab", icon = icon("house")),
      menuItem("Discover", tabName = "gene_effect", icon = icon("dna")),
      menuItem("Target", tabName = "gene_tab", icon = icon("magnifying-glass")),
      menuItem("Chemistry", tabName = "chemistry_tab", icon = icon("flask")),
      menuItem("About", tabName = "about_tab", icon = icon("circle-info"))
    )
  ),
  dashboardBody(
    useShinyjs(),
    tags$head(tags$style(HTML(build_app_css()))),
    tabItems(
      # Home tab
      tabItem(
        tabName = "home_tab",
        fluidRow(
          box(
            width = 12, status = NULL, solidHeader = FALSE,
            style = "background: transparent; box-shadow: none;",
            div(
              class = "cpt-home-hero",
              style = "text-align: center; max-width: 960px; margin: 0 auto;",
              h2("Welcome to CanProTarget"),
              p(
                style = "color: #444; line-height: 1.6;",
                "Cancer chemoproteomics-based target prioritization pipeline."
              ),
              tags$img(
                src = "graphical_abstract.png",
                alt = "CanProTarget graphical abstract — pipeline overview",
                style = "max-width: 100%; height: auto; border: 1px solid #ddd; border-radius: 4px; background: #fff;"
              )
            )
          )
        )
      ),
      # Discover tab — layered target explorer; hosts the cysteine atlas browse
      dependencies_ui("deps"),
      # Target tab — one target across gene, protein and residue evidence
      gene_ui("gene"),
      # Chemistry tab — SwissADME only. Single-protein
      # probes are on the Target tab, and multi-protein lookup is the Discover
      # target explorer, which takes a pasted gene list and shows engaged
      # cysteines, competition ratios, target counts and probe selectivity.
      tabItem(
        tabName = "chemistry_tab",
        swissadme_ui("swiss", wrap = FALSE)
      ),
      # About tab (data versions filled from data/data_versions.yaml in server)
      tabItem(
        tabName = "about_tab",
        fluidRow(
          box(
            width = 12, status = "primary", solidHeader = TRUE,
            title = "About CanProTarget",
            div(
              style = "max-width: 960px; line-height: 1.7;",
              h3("Overview"),
              p("CanProTarget (cancer chemoproteomics-based target prioritization
                pipeline) allows for the combination of three kinds of evidence from
                three different platforms:"),
              tags$ol(
                tags$li("Cancer dependency map (DepMap) functional genomics data to
                        reveal dependencies within a cancer context."),
                tags$li("Cysteinome database (CysDB) chemoproteomic data to reveal
                        covalent ligands."),
                tags$li("Cysteine editing (Cys_editing) to rank potential targets by
                        residue-level functionality.")
              ),

              h3("Tabs"),
              tags$ul(
                tags$li(tags$strong("Home:"), " provides a graphical abstract
                  summarising the basis of each platform CanProTarget leverages."),
                tags$li(tags$strong("Discover:"), " allows for the analysis of the
                  evidence layers (dependency, ligandability, residue evidence) either
                  alone or concurrently."),
                tags$li(tags$strong("Target:"), " information on a target including
                  where it is a dependency, distribution of gene effect scores across
                  cell lines, the tiers of its engaged cysteines and the probes that
                  engage them."),
                tags$li(tags$strong("Chemistry:"), " SwissADME drug-likeness and
                  physicochemistry of covalent ligands with their protein-binding
                  profile.")
              ),

              h3("Data versions"),
              uiOutput("about_data_versions"),

              h3("Citation"),
              tags$dl(
                style = "margin-left: 20px;",
                tags$dt(style = "font-weight: 600; color: #1a3a5c;", "Software"),
                tags$dd(style = "margin-bottom: 12px;",
                        "Ong JP, Martins D, Bertoldo JB. ",
                        tags$em("CanProTarget"), ". Zenodo. ",
                        tags$span(style = "color: #7b8a94;", "DOI pending release.")),
                tags$dt(style = "font-weight: 600; color: #1a3a5c;", "Article"),
                tags$dd("Ong JP, Bell J, Martins D, Zhu J, Rodrigues T, Bertoldo JB. ",
                        tags$em("CanProTarget: a residue-resolved prioritization
                        platform for covalent cancer drug discovery"), ". ",
                        tags$span(style = "color: #7b8a94;", "In preparation."))
              ),
              h3("Licence"),
              p("Released under the GNU Affero General Public License v3.0."),

              h3("Source Code"),
              p(tags$a(href = "https://github.com/Bertoldo-Lab/CanProTarget", target = "_blank",
                       icon("github"), " GitHub Repository"))
            )
          )
        )
      )
    )
  )
)

# ---- 5. Server definition (shared data + module wiring) -------
server <- function(input, output, session) {

  # ---- Bookmarking: custom "Share This View" button ---------------
  # Instead of Shiny's default modal (shows ugly long URL), we copy to clipboard
  # and show a brief notification.
  onBookmark(function(state) {
    state$values$active_tab <- input$tabs
  })

  onBookmarked(function(url) {
    # Send URL to clipboard via JS and show notification
    js <- sprintf(
      "navigator.clipboard.writeText('%s').then(function() {}, function() {});",
      gsub("'", "\\\\'", url)
    )
    shinyjs::runjs(js)
    showNotification(
      "Link copied to clipboard. Share it to restore this exact view.",
      type = "message", duration = 4
    )
  })

  onRestored(function(state) {
    if (!is.null(state$values$active_tab)) {
      updateTabItems(session, "tabs", selected = state$values$active_tab)
    }
  })

  # Trigger bookmark when custom button is clicked

  observeEvent(input$share_view, {
    session$doBookmark()
  })

  # ---- Shared data reactives (loaded once, passed to modules) ----

  # SwissADME table (offline-built RDS only; see docs/DATA_PROVENANCE.md)
  swissadme_table <- reactive({
    tryCatch(
      process_swissadme_data(data_dir),
      error = function(e) {
        cpt_notify_error("CPT-1009", "SwissADME data load failed",
                         detail = conditionMessage(e))
        tibble::tibble()
      }
    )
  }) %>%
    bindCache(swissadme_data_sig(data_dir))

  # Engagement-only slice of the binding table: the 30,219 rows at CR >= 4, out
  # of 10.6 million. 3 MB against 555 MB, so it is read at startup and serves
  # every view that uses the CR >= 4 threshold. The full table below is read
  # only when a CR floor under 4 asks for measurements this file does not hold.
  protein_binding_cr4 <- reactive({
    p <- data_path("protein_binding_cr4.rds")
    if (!file.exists(p)) return(NULL)
    tryCatch({
      x <- readRDS(p)
      # Same normalisation and indexing the full table gets below. The stored
      # file keeps the source naming (ACRYL_0), while every probe name in the
      # app comes from the SwissADME table, which is simplified on load (AC0).
      # Skipping this leaves the two vocabularies unable to join.
      x <- x %>%
        dplyr::mutate(probe_name = as.factor(simplify_probe_name(
          as.character(.data$probe_name))))
      if ("Dataset" %in% names(x)) x$Dataset <- cpt_fix_dataset_names(x$Dataset)
      cpt_index_protein_binding(x)
    }, error = function(e) NULL)
  })

  # Protein binding lookup (~10M+ rows). Loaded on first use (Find Probes / Load Binding /
  # Ligandability probes), not at session start. SMILES aligned offline when building the RDS.
  proteinbindinglookup_raw <- reactive({
    # Either form will do: the factored copy is preferred, the original is
    # converted on the fly. Only the absence of both is a problem.
    lookup_path <- data_path("protein_binding_lookup_preprocessed.rds")
    fast_path <- data_path("protein_binding_lookup_factored.rds")
    if (!file.exists(fast_path) && !file.exists(lookup_path)) {
      cpt_notify_error("CPT-1007",
        "Protein binding data not available",
        detail = "See docs/DATA_PROVENANCE.md for setup instructions")
      return(NULL)
    }
    tryCatch(
      {
        # 10.6 million rows: every character column costs ~85 MB regardless of
        # what is in it, because R stores one pointer per element. Several of
        # these columns hold a handful of distinct values — Cell_Line has 8,
        # ligandable has 2 — so as factors they become integer vectors plus a
        # short level table: 936 MB -> 555 MB, with identical behaviour
        # (comparisons, %in%, match() and split() are unaffected, and the two
        # places that do string work already call as.character() first).
        #
        # Storing the file that way is better still, because the conversion
        # then happens once offline instead of on every session: the factored
        # RDS loads in 1.8 s against 22.8 s for the original. Use it when it
        # exists and fall back to converting on the fly when it does not, so
        # the app works either way.
        lookup <- if (file.exists(fast_path)) {
          readRDS(fast_path)
        } else {
          x <- readRDS(lookup_path)
          chr <- vapply(x, is.character, logical(1))
          if (any(chr)) x[chr] <- lapply(x[chr], as.factor)
          x
        }
        lookup <- lookup %>%
          dplyr::mutate(probe_name = as.factor(simplify_probe_name(
            as.character(.data$probe_name))))
        if ("Dataset" %in% names(lookup)) {
          lookup$Dataset <- cpt_fix_dataset_names(lookup$Dataset)
        }

        cpt_index_protein_binding(lookup)
      },
      error = function(e) {
        cpt_notify_error("CPT-1008", "Protein binding data corrupt",
                         detail = conditionMessage(e))
        NULL
      }
    )
  }) %>%
    bindCache({
      # Whichever copy is present; the factored one is what gets read.
      vapply(c("protein_binding_lookup_factored.rds",
               "protein_binding_lookup_preprocessed.rds"),
             function(f) { p <- data_path(f)
                           if (file.exists(p)) as.numeric(file.mtime(p)) else 0 },
             numeric(1))
    })

  # Cancer model metadata (cancer_model_data.rds only). Build offline -- docs/DATA_PROVENANCE.md.
  cancer_model_data <- reactive({
    p <- data_path("cancer_model_data.rds")
    if (!file.exists(p)) {
      cpt_notify_error("CPT-1001",
        "Cancer model metadata not found",
        detail = "Run docs/scripts/preprocess_data.R to generate cancer_model_data.rds")
      return(NULL)
    }
    tryCatch(
      readRDS(p),
      error = function(e) {
        cpt_notify_error("CPT-1002", "Cancer model metadata corrupt",
                         detail = conditionMessage(e))
        NULL
      }
    )
  }) %>%
    bindCache(cancer_model_data_sig(data_dir))

  # Publication-derived, cysteine-level atlas imported offline from
  # cravattlab/Cys_editing. This compact RDS is safe to load on first use.
  cys_editing_atlas <- reactive({
    path <- data_path("cys_editing_atlas.rds")
    if (!file.exists(path)) {
      cpt_notify_error("CPT-1005",
        "Cysteine editing atlas not found",
        detail = "Run: Rscript docs/scripts/import_cys_editing_atlas.R")
      return(NULL)
    }
    tryCatch(
      readRDS(path),
      error = function(e) {
        cpt_notify_error("CPT-1006", "Cysteine editing atlas corrupt",
                         detail = conditionMessage(e))
        NULL
      }
    )
  }) %>%
    bindCache({
      p <- data_path("cys_editing_atlas.rds")
      if (file.exists(p)) file.mtime(p) else 0
    })

  # Per-gene index behind the Gene tab (docs/scripts/build_gene_index.R).
  # Small on purpose: it exists so gene-level questions never touch the
  # 140 MB matrix or the 10.6M-row binding table.
  gene_index <- reactive({
    p <- data_path("gene_index.rds")
    if (!file.exists(p)) {
      cpt_notify_error(
        "CPT-1011",
        "Gene index not available",
        detail = "Run: Rscript docs/scripts/build_gene_index.R",
        type = "warning",
        duration = 12
      )
      return(NULL)
    }
    tryCatch(
      readRDS(p),
      error = function(e) {
        cpt_notify_error("CPT-1012", "Gene index corrupt",
                         detail = conditionMessage(e))
        NULL
      }
    )
  }) %>%
    bindCache({
      p <- data_path("gene_index.rds")
      if (file.exists(p)) file.mtime(p) else 0
    })

  # Data version manifest (About tab + modules / reports provenance)
  data_versions <- reactive({
    read_data_versions(data_dir)
  })

  output$about_data_versions <- renderUI({
    cpt_data_versions_ui(data_versions())
  })

  # Comprehensive shared_data list passed to all modules.
  # Canonical key: protein_binding_lookup. proteinbindinglookup kept as alias.
  # Per-gene ADME lookup collapsed from the binding table + SwissADME.
  # Lazy: only built when a CPT score with ADME is first requested, because the
  # binding table is ~10.6M rows. Cached on the file mtimes of both inputs.
  adme_gene_scores <- reactive({
    bind <- proteinbindinglookup_raw()
    sw <- swissadme_table()
    if (is.null(bind) || is.null(sw) || !nrow(bind) || !nrow(sw)) return(NULL)
    tryCatch(
      cpt_build_adme_gene_scores(bind, sw),
      error = function(e) {
        cpt_notify_error("CPT-1010", "ADME lookup build failed",
                         detail = conditionMessage(e))
        NULL
      }
    )
  }) %>%
    bindCache({
      p1 <- data_path("protein_binding_lookup_factored.rds")
      if (!file.exists(p1)) p1 <- data_path("protein_binding_lookup_preprocessed.rds")
      p2 <- data_path("swissadme_preprocessed.rds")
      c(if (file.exists(p1)) file.mtime(p1) else 0,
        if (file.exists(p2)) file.mtime(p2) else 0)
    })

  shared_data <- list(
    protein_binding_lookup    = proteinbindinglookup_raw,
    protein_binding_cr4       = protein_binding_cr4,
    proteinbindinglookup      = proteinbindinglookup_raw,
    swissadme_table           = swissadme_table,
    adme_gene_scores          = adme_gene_scores,
    cancer_model_data         = cancer_model_data,
    cys_editing_atlas         = cys_editing_atlas,
    gene_index                = gene_index,
    data_versions             = data_versions,
    file_map                  = file_map,
    cancer_subtype_list_files = cancer_subtype_list_files,
    data_path                 = data_path,
    data_dir                  = data_dir
  )

  # ---- Wire up module servers ------------------------------------
  # Gene -> Discover hand-off: selecting a subtype row on the Gene tab sets up
  # and runs that analysis on Discover, instead of making the user re-find it.
  discover_preset <- reactiveVal(NULL)
  gene_preset <- reactiveVal(NULL)

  dependencies_outputs <- dependencies_server(
    "deps", shared_data,
    preset = discover_preset,
    on_open_gene = function(gene_key) {
      gene_preset(gene_key)
      updateTabItems(session, "tabs", selected = "gene_tab")
    }
  )

  gene_server("gene", shared_data, on_jump = function(preset) {
    discover_preset(preset)
    updateTabItems(session, "tabs", selected = "gene_effect")
  }, preset = gene_preset,
     dep_effect_min = dependencies_outputs$effect_min,
     dep_context = dependencies_outputs)
  cys_editing_server(
    "cys", shared_data,
    dependency_genes = dependencies_outputs$cancer_genes,
    filters = dependencies_outputs$residue_filters
  )
  swissadme_server("swiss", shared_data)
}

# ---- Run app ---------------------------------------------------
shinyApp(ui, server, enableBookmarking = "url")
