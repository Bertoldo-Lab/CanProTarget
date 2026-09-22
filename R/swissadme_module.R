# ============================================================
# Script:   swissadme_module.R
# Purpose:  Shiny module for SwissADME tab (CanProTarget):
#           - display ADME/drug-likeness properties per probe
#           - show 2D structure (PubChem/ChemSpider fallback)
#           - protein binding tables by dataset and cell line
#           - BOILED-Egg and radar plots
# Inputs:   shared_data from app.R (swissadme_table, protein_binding_lookup)
# Outputs:  UI + server functions for tabName = "swissadme_tab"
# ============================================================
# Sections:
#   1. swissadme_ui()
#   2. swissadme_server()
#      2.1 Probe selector
#      2.2 2D structure pane
#      2.3 Molecular property tables
#      2.4 BOILED-Egg plot
#      2.5 Radar / drug-likeness plot
#      2.6 Protein binding tables and downloads
#      2.7 Full SwissADME table download
# ============================================================

library(shiny)
library(shinydashboard)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)
library(DT)
library(writexl)

# Plotly polar default: counterclockwise, rotation 0 → first category at East (3 o'clock).
# Cartesian replica for PNG export (ggplot coord_polar maps discrete axes differently).
radar_ccw_from_east_xy <- function(r_per_spoke) {
  n <- length(r_per_spoke)
  r <- c(r_per_spoke, r_per_spoke[[1L]])
  m <- length(r)
  theta <- (seq_len(m) - 1L) * (2 * pi / n)
  data.frame(x = r * cos(theta), y = r * sin(theta))
}

radar_axis_label_xy <- function(labels, r = 1.12) {
  n <- length(labels)
  theta <- (seq_len(n) - 1L) * (2 * pi / n)
  data.frame(
    x = r * cos(theta),
    y = r * sin(theta),
    label = labels,
    stringsAsFactors = FALSE
  )
}

# ============================================================
# 1. swissadme_ui
# ============================================================

#' @param wrap TRUE returns a standalone tabItem (legacy top-level tab);
#'   FALSE returns just the body so it can be embedded in the Chemistry tab.
swissadme_ui <- function(id, wrap = TRUE) {
  ns <- NS(id)

  body <- tagList(
    fluidRow(
      # ---- Left control panel --------------------------------
      box(
        title = "Probe & Filter Controls", width = 3, status = "primary",
        solidHeader = TRUE,

        # Probe selector (populated from SwissADME table)
        cpt_label("Select Probe", for_id = ns("probe_select")),
        selectInput(ns("probe_select"), label = NULL,
          choices  = character(0),
          selectize = TRUE),
        hr(),

        # CR cutoff for counting "targets" in the protein binding tables
        cpt_label("Minimum Competition Ratio (CR)", for_id = ns("cr_cutoff")),
        numericInput(ns("cr_cutoff"), label = NULL,
          value = 4, min = 0, step = 0.5),
        hr(),

        # Explicit load button for protein binding + full table
        actionButton(
          ns("load_binding"),
          "Load Protein Binding",
          class = "btn-sm btn-primary btn-block"
        ),
        br(),
        actionButton(
          ns("reset_swiss"),
          "Reset Controls",
          icon = icon("rotate-left"),
          class = "btn-sm btn-default btn-block"
        )
      ),

      # ---- Main display area ---------------------------------
      box(
        title = "Chemistry Explorer", width = 9, status = "info",
        solidHeader = TRUE,
        uiOutput(ns("data_unavailable")),
        tabsetPanel(
          id = ns("swiss_tabs"),

          # -- Combined Explorer: properties, structure, plots --
          tabPanel(
            "Explorer",
            div(
              class = "swissadme-explorer-pane",
              br(),
              uiOutput(ns("probe_loaded_title")),
              br(),
              # Top: 2D structure (left) + properties (right)
              fluidRow(
                column(
                  width = 6,
                  h4("2D Structure"),
                  br(),
                  uiOutput(ns("structure_img"))
                ),
                column(
                  width = 6,
                  h4("Molecular Properties"),
                  br(),
                  h5("Molecular Formula & Size"),
                  tableOutput(ns("tbl_formula")),
                  br(),
                  h5("Lipophilicity"),
                  tableOutput(ns("tbl_logp")),
                  br(),
                  h5("Solubility & ADME"),
                  tableOutput(ns("tbl_solubility")),
                  br(),
                  h5("Permeability & Alerts"),
                  tableOutput(ns("tbl_permeability"))
                )
              ),
              hr(),
              # Bottom: plots side-by-side
              fluidRow(
                column(
                  width = 6,
                  cpt_section("BOILED-Egg", level = "h4"),
                  radioButtons(
                    ns("boiled_show"),
                    label   = NULL,
                    choices = c(
                      "All probes (highlight selected)" = "all",
                      "Selected probe only"             = "selected_only"
                    ),
                    selected = "all",
                    inline = FALSE
                  ),
                  br(),
                  plotlyOutput(ns("boiled_egg"), height = "450px")
                ),
                column(
                  width = 6,
                  cpt_section("Drug-likeness Radar", level = "h4"),
                  br(),
                  plotlyOutput(ns("radar_plot"), height = "550px")
                )
              ),
              hr(),
              fluidRow(
                column(
                  width = 6,
                  downloadButton(ns("dl_boiled"),
                    "BOILED-Egg (PNG)",
                    class = "btn-sm btn-default btn-block")
                ),
                column(
                  width = 6,
                  downloadButton(ns("dl_radar"),
                    "Radar (PNG)",
                    class = "btn-sm btn-default btn-block")
                )
              )
            )
          ),

          # -- Protein binding tables --------------------------
          tabPanel(
            "Protein Binding",
            br(),
            conditionalPanel(
              condition = paste0("output['", ns("binding_ready"), "'] != '1'"),
              div(style = "min-height: 80px;")
            ),
            conditionalPanel(
              condition = paste0("output['", ns("binding_ready"), "'] == '1'"),
              tagList(
                uiOutput(ns("swiss_status_line")),
                br(),
                h5("All datasets combined"),
                shinycssloaders::withSpinner(
                  DTOutput(ns("pb_all")),
                  type = 6, color = "#2C6A94"),
                br(),
                downloadButton(ns("dl_protein_binding_pb_tab"),
                  "Download Protein Binding (.xlsx)",
                  class = "btn-sm btn-default"),
                br(), br(),
                shinycssloaders::withSpinner(
                  plotlyOutput(ns("pb_binding_scatter"), height = "420px"),
                  type = 6, color = "#2C6A94"),
                br(),
                downloadButton(ns("dl_pb_binding_plot"),
                  "Download Probe Protein Binding (PNG)",
                  class = "btn-sm btn-default")
              )
            )
          ),

          # -- Full SwissADME table ----------------------------
          tabPanel(
            "Full SwissADME Table",
            br(),
            DTOutput(ns("full_adme_table")),
            br(),
            downloadButton(ns("dl_swissadme_full_tab"),
              "Download Full SwissADME (.xlsx)",
              class = "btn-sm btn-default")
          )
        )
      )
    )
  )

  if (isTRUE(wrap)) tabItem(tabName = "swissadme_tab", body) else body
}

# ============================================================
# 2. swissadme_server
# ============================================================

swissadme_server <- function(id, shared_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    swiss_last_updated <- reactiveVal(NULL)
    binding_ready <- reactiveVal(FALSE)
    pb_loaded_val <- reactiveVal(NULL)
    # Lowest CR the currently loaded table can answer for. The CR >= 4 subset
    # is a strict row subset of the full table, so it is only valid at or above
    # that floor; the full table is valid everywhere (-Inf).
    pb_loaded_floor <- reactiveVal(Inf)

    pb_full_lookup <- function() {
      if (!is.null(shared_data$protein_binding_lookup)) {
        shared_data$protein_binding_lookup()
      } else {
        shared_data$proteinbindinglookup()
      }
    }

    # At the default cutoff the 30k-row subset answers the same question as the
    # 10.6M-row table, so read it instead and skip the slow load entirely.
    pb_lookup <- function(cr_floor) {
      if (is.finite(cr_floor) && cr_floor >= 4 &&
          !is.null(shared_data$protein_binding_cr4)) {
        cr4 <- shared_data$protein_binding_cr4()
        if (is.data.frame(cr4) && nrow(cr4)) {
          pb_loaded_floor(4)
          return(cr4)
        }
      }
      pb_loaded_floor(-Inf)
      pb_full_lookup()
    }

    output$binding_ready <- renderText({
      if (isTRUE(binding_ready())) "1" else "0"
    })
    outputOptions(output, "binding_ready", suspendWhenHidden = FALSE)

    output$data_unavailable <- renderUI({
      path_fn <- shared_data$data_path
      p <- if (is.function(path_fn)) {
        path_fn("swissadme_preprocessed.rds")
      } else {
        file.path(shared_data$data_dir %||% "data", "swissadme_preprocessed.rds")
      }
      if (file.exists(p)) return(NULL)
      cpt_data_unavailable_card(
        title = "SwissADME data not available",
        detail = paste(
          "Probe ADME properties need swissadme_preprocessed.rds.",
          "Cancer Dependencies and Functional Cysteines do not require this file."
        ),
        code = "CPT-1009"
      )
    })

    # ---- 2.1 Probe selector ------------------------------------
    # Populate the probe selector from the SwissADME table
    observe({
      adme <- shared_data$swissadme_table()
      req(adme, nrow(adme) > 0, "probe_name" %in% colnames(adme))
      probes <- unique(adme$probe_name)
      probes_sorted <- gtools::mixedsort(probes)
      updateSelectInput(session, "probe_select",
        choices  = probes_sorted,
        selected = probes_sorted[1])
    })

    # Row(s) for the currently selected probe
    selected_adme <- reactive({
      adme <- shared_data$swissadme_table()
      req(adme, input$probe_select)
      filter_swissadme_data(adme, input$probe_select)
    })

    output$probe_loaded_title <- renderUI({
      adme <- selected_adme()
      req(adme, nrow(adme) > 0, input$probe_select)
      tags$div(
        style = "font-size: 2.1rem; font-weight: 700; color: #1f5e8c; letter-spacing: 0.02em;",
        simplify_probe_name(input$probe_select)
      )
    })

  # ---- 2.2 2D structure pane (improved) ----------------------
  # Uses SMILES-based rendering with PubChem → ChemSpider → CDK Depict fallback
  output$structure_img <- renderUI({
    req(selected_adme())

    smiles <- selected_adme()$`Canonical.SMILES`[1]
    req(smiles)

    encoded_smiles <- utils::URLencode(smiles, reserved = TRUE)

    # PubChem (primary, most reliable for simple structures)
    pubchem_url <- paste0(
      "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/",
      encoded_smiles,
      "/PNG?image_size=large"
    )

    # ChemSpider fallback (handles more complex structures)
    chemspider_url <- paste0(
      "https://www.chemspider.com/ImagesHandler.ashx?w=400&h=400&smiles=",
      encoded_smiles
    )

    # CDK Depict fallback (SVG, open-source renderer)
    cdk_url <- paste0(
      "https://www.simolecule.com/cdkdepict/depict/bow/svg?smi=",
      encoded_smiles,
      "&annotate=colmap&zoom=3"
    )

    tags$div(
      tags$img(
        src = pubchem_url,
        alt = "2D Chemical Structure",
        id = ns("chemStructureImg"),
        style = "width:100%; max-width:400px; border: 1px solid #ddd; padding: 5px; background: white;",
        onerror = paste0(
          "this.onerror=function(){",
          "  this.onerror=null;",
          "  this.src='", cdk_url, "';",
          "  document.getElementById('", ns("rendererSource"), "').innerHTML='<strong>Renderer:</strong> CDK Depict (Chemistry Development Kit)';",
          "};",
          "this.src='", chemspider_url, "';",
          "document.getElementById('", ns("rendererSource"), "').innerHTML='<strong>Renderer:</strong> ChemSpider';"
        ),
        onload = paste0(
          "if(this.src.includes('pubchem')) ",
          "document.getElementById('", ns("rendererSource"), "').innerHTML='<strong>Renderer:</strong> PubChem';"
        )
      ),
      tags$div(
        style = "margin-top: 10px; padding: 8px; background-color: #f8f9fa; border-radius: 4px;",
        tags$p(
          style = "margin: 0; font-size: 14px; color: #495057;",
          tags$strong("SMILES: "),
          tags$code(
            style = "background: #e9ecef; padding: 2px 6px; border-radius: 3px; font-size: 13px;",
            smiles
          )
        ),
        tags$p(
          id = ns("rendererSource"),
          style = "margin: 5px 0 0 0; font-size: 13px; color: #6c757d;",
          tags$strong("Renderer: "), "PubChem"
        )
      )
    )
  })
    # ---- 2.3 Molecular property tables (condensed) --------------
    # Helper: extract named rows for a given set of column-name patterns
    # Returns only the first matching occurrence of each pattern to keep
    # the tables compact and focused on key descriptors.
    prop_table <- function(df, col_patterns) {
      if (is.null(df) || nrow(df) == 0) {
        return(data.frame(Property = character(), Value = character()))
      }

      cols <- colnames(df)
      matched <- character(0)

      for (pat in col_patterns) {
        hits <- grep(pat, cols, value = TRUE, ignore.case = TRUE)
        if (length(hits) > 0) {
          # Take only the first hit for each pattern to avoid clutter
          matched <- c(matched, hits[1])
        }
      }

      matched <- unique(matched)

      if (length(matched) == 0) {
        return(data.frame(Property = character(), Value = character()))
      }

      data.frame(
        Property = matched,
        Value    = unlist(df[1, matched, drop = TRUE]),
        stringsAsFactors = FALSE
      )
    }

    # -- Very compact "Formula & Size": just formula and MW -------
    output$tbl_formula <- renderTable({
      req(selected_adme())
      prop_table(
        selected_adme(),
        # Prefer "Formula" and "MW" / "Molecular weight"
        c("^formula$", "^mw$", "molecular.*weight")
      )
    }, rownames = FALSE)

    # -- Lipophilicity: single consensus / primary logP ----------
    output$tbl_logp <- renderTable({
      req(selected_adme())
      prop_table(
        selected_adme(),
        c(
          "consensus.*logp",   # SwissADME consensus, if present
          "^wlogp$",           # WLOGP
          "^xlogp3$",          # XLogP3
          "logp"               # any remaining generic logP
        )
      )
    }, rownames = FALSE)

    # -- Solubility: one ESOL class/value, optionally another ----
    output$tbl_solubility <- renderTable({
      req(selected_adme())
      prop_table(
        selected_adme(),
        c(
          "esol.*class",       # ESOL solubility class
          "esol.*log\\.s",     # ESOL LogS
          "log.*s"             # any generic logS if ESOL not found
        )
      )
    }, rownames = FALSE)

    # -- Permeability & key transport/alert flags ----------------
    output$tbl_permeability <- renderTable({
      req(selected_adme())
      prop_table(
        selected_adme(),
        c(
          "^bbb",              # BBB permeant
          "gi.*absorption",    # GI absorption
          "pgp",               # P-gp substrate
          "caco",              # Caco-2 permeability
          "pains|alert"        # a single PAINS / structural alert summary
        )
      )
    }, rownames = FALSE)

    # ---- 2.4 BOILED-Egg plot (SwissADME geometry) --------------
    # Pure builder: cpt_boiled_egg_gg() in functions.R (unit-tested).
    boiled_gg <- reactive({
      adme <- shared_data$swissadme_table()
      req(adme, nrow(adme) > 0)
      show_mode <- if (is.null(input$boiled_show)) "all" else input$boiled_show
      cpt_boiled_egg_gg(
        adme,
        selected_probe = input$probe_select,
        show_mode = show_mode,
        theme_fn = cpt_theme
      )
    })

    output$boiled_egg <- renderPlotly({
      cpt_plotly(ggplotly(boiled_gg(), tooltip = c("text", "x", "y")),
                 margin = list(t = 30))
    })

    output$dl_boiled <- downloadHandler(
      filename = function() paste0("boiled_egg_", Sys.Date(), ".png"),
      content  = function(file) {
        ggplot2::ggsave(file, plot = boiled_gg(),
          width = 9, height = 6, dpi = 300)
      }
    )

    # ---- 2.5 Radar / drug-likeness plot (SwissADME-style) -----
    normalize <- function(x, min, max) {
      pmin(pmax((x - min) / (max - min), 0), 1)
    }

    radar_values <- reactive({
      req(selected_adme())
      data <- selected_adme()

      # Raw values
      lipo   <- as.numeric(data$WLOGP[1])
      size   <- as.numeric(data$MW[1])
      polar  <- as.numeric(data$TPSA[1])
      insolu <- as.numeric(data$ESOL.Log.S[1])
      insatu <- as.numeric(data$Fraction.Csp3[1])
      flex   <- as.numeric(data$`X.Rotatable.bonds`[1])

      # Normalized 0–1 using SwissADME limits
      vals_norm <- c(
        normalize(lipo,   -0.7, 5.0),   # LIPO
        normalize(size,   150,  500),   # SIZE
        normalize(polar,  20,   130),   # POLAR
        normalize(insolu, -6,   0),     # INSOLU
        normalize(insatu, 0.25, 1.0),   # INSATU
        normalize(flex,   0,    9)      # FLEX
      )

      labels <- c("LIPO", "SIZE", "POLAR", "INSOLU", "INSAT", "FLEX")

      tibble::tibble(
        Axis      = labels,
        NormValue = vals_norm
      )
    })

    # Optimal zone bounds in normalized units (inner/outer)
    radar_optimal_bounds <- reactive({
      tibble::tibble(
        Axis = c("LIPO", "SIZE", "POLAR", "INSOLU", "INSAT", "FLEX"),
        Inner = c(0.1, 0.2, 0.1, 0.3, 0.4, 0.1),
        Outer = c(0.9, 0.9, 0.7, 0.8, 1.0, 0.8)
      )
    })

    radar_plotly_obj <- reactive({
      vals <- radar_values()
      bounds <- radar_optimal_bounds()
      req(vals, bounds)

      labels <- vals$Axis

      probe_r   <- c(vals$NormValue, vals$NormValue[1])
      probe_th  <- c(labels, labels[1])

      outer_r   <- bounds$Outer
      inner_r   <- bounds$Inner

      text_probe <- paste(labels, ":", round(vals$NormValue, 2))
      text_probe <- c(text_probe, text_probe[1])

      plot_ly(
        type = "scatterpolar",
        mode = "lines",
        fill = "toself"
      ) %>%
        add_trace(
          r = c(outer_r, outer_r[1]),
          theta = c(labels, labels[1]),
          fill = "none",
          line = list(
            color = "rgba(125,139,151,0.9)",
            width = 2,
            dash = "dash"
          ),
          name = "Optimal outer",
          hoverinfo = "none"
        ) %>%
        add_trace(
          r = c(inner_r, inner_r[1]),
          theta = c(labels, labels[1]),
          fill = "none",
          line = list(
            color = "rgba(255,182,193,0.9)",
            width = 1.5,
            dash = "dash"
          ),
          name = "Optimal inner",
          hoverinfo = "none"
        ) %>%
        add_trace(
          r = probe_r,
          theta = probe_th,
          fillcolor = "rgba(26,58,92,0.25)",
          line = list(color = CPT_PAL$ink, width = 3),
          marker = list(size = CPT_MARKER_SIZE, color = CPT_PAL$ink),
          mode = "markers+lines",
          name = "Selected probe",
          hoverinfo = "text",
          text = text_probe
        ) %>%
        layout(
          polar = list(
            radialaxis = list(visible = TRUE, range = c(0, 1)),
            angularaxis = list(
              direction = "counterclockwise",
              rotation = 0,
              categoryorder = "array",
              categoryarray = as.list(labels)
            )
          ),
          showlegend = FALSE,
          title = list(text = ""),
          margin = list(l = 80, r = 80, t = 60, b = 40)
        )
    })

    output$radar_plot <- renderPlotly({
      cpt_plotly(radar_plotly_obj())
    })

    output$dl_radar <- downloadHandler(
      filename = function() {
        paste0("radar_", input$probe_select, "_", Sys.Date(), ".png")
      },
      content = function(file) {
        vals   <- radar_values()
        bounds <- radar_optimal_bounds()
        req(vals, bounds)

        axis_levels <- c("LIPO", "SIZE", "POLAR", "INSOLU", "INSAT", "FLEX")
        b <- bounds %>%
          dplyr::mutate(Axis = factor(.data$Axis, levels = axis_levels)) %>%
          dplyr::arrange(.data$Axis)
        v <- vals %>%
          dplyr::mutate(Axis = factor(.data$Axis, levels = axis_levels)) %>%
          dplyr::arrange(.data$Axis)

        df_outer <- radar_ccw_from_east_xy(as.numeric(b$Outer))
        df_inner <- radar_ccw_from_east_xy(as.numeric(b$Inner))
        df_probe <- radar_ccw_from_east_xy(as.numeric(v$NormValue))
        verts <- df_probe[-nrow(df_probe), , drop = FALSE]
        lbl <- radar_axis_label_xy(as.character(axis_levels), r = 1.14)

        grid_rs <- c(0.25, 0.5, 0.75, 1)
        circle_pts <- dplyr::bind_rows(lapply(grid_rs, function(r0) {
          th <- seq(0, 2 * pi, length.out = 144L)
          data.frame(
            x = r0 * cos(th),
            y = r0 * sin(th),
            grp = sprintf("g%.2f", r0)
          )
        }))
        spokes <- dplyr::bind_rows(lapply(seq_along(axis_levels), function(i) {
          th <- (i - 1L) * (2 * pi / length(axis_levels))
          data.frame(
            x = c(0, cos(th)),
            y = c(0, sin(th)),
            grp = i
          )
        }))

        pink_line <- grDevices::adjustcolor("#FF69B4", alpha.f = 0.92)
        # Was hot pink; the reference rings are furniture, not data.
        pink_in <- grDevices::adjustcolor("#aab4bd", alpha.f = 0.92)

        p <- ggplot2::ggplot() +
          ggplot2::geom_path(
            data = circle_pts,
            ggplot2::aes(x = .data$x, y = .data$y, group = .data$grp),
            colour = grDevices::adjustcolor("grey45", alpha.f = 0.35),
            linewidth = 0.25
          ) +
          ggplot2::geom_path(
            data = spokes,
            ggplot2::aes(x = .data$x, y = .data$y, group = .data$grp),
            colour = grDevices::adjustcolor("grey50", alpha.f = 0.4),
            linewidth = 0.35
          ) +
          ggplot2::geom_path(
            data = df_outer,
            ggplot2::aes(x = .data$x, y = .data$y),
            colour = pink_line,
            linetype = "dashed",
            linewidth = 1
          ) +
          ggplot2::geom_path(
            data = df_inner,
            ggplot2::aes(x = .data$x, y = .data$y),
            colour = pink_in,
            linetype = "dashed",
            linewidth = 0.85
          ) +
          ggplot2::geom_polygon(
            data = df_probe,
            ggplot2::aes(x = .data$x, y = .data$y),
            fill = grDevices::adjustcolor(CPT_PAL$ink, alpha.f = 0.3),
            colour = CPT_PAL$ink,
            linewidth = 1.1
          ) +
          ggplot2::geom_point(
            data = verts,
            ggplot2::aes(x = .data$x, y = .data$y),
            colour = CPT_PAL$ink,
            size = 4
          ) +
          ggplot2::geom_text(
            data = lbl,
            ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
            size = 3.2,
            colour = "grey20"
          ) +
          ggplot2::coord_fixed(
            xlim = c(-1.3, 1.3),
            ylim = c(-1.3, 1.3),
            expand = FALSE
          ) +
          ggplot2::labs(x = NULL, y = NULL, title = "Drug-likeness radar") +
          ggplot2::theme_void(base_size = 12) +
          ggplot2::theme(
            plot.title = ggplot2::element_text(
              hjust = 0.5,
              face = "plain",
              margin = ggplot2::margin(b = 6)
            ),
            plot.background = ggplot2::element_rect(fill = "white", colour = NA)
          )

        ggplot2::ggsave(file, plot = p, width = 7, height = 7, dpi = 300, bg = "white")
      }
    )

    # ---- 2.6 Protein binding tables and downloads --------------
    # Toast + compute into pb_loaded_val on this click. Do not wrap the RDS
    # load in eventReactive(input$load_binding): the first click would see NULL.
    observeEvent(input$load_binding, {
      req(input$load_binding > 0)
      cpt_log_usage(
        "swiss_load_binding",
        list(
          probe = input$probe_select,
          cr_cutoff = input$cr_cutoff
        )
      )
      showNotification(
        "Loading protein binding data...",
        type = "message",
        duration = NULL,
        id = "load_swiss_pb"
      )
      cr_floor <- suppressWarnings(as.numeric(input$cr_cutoff))
      if (!length(cr_floor) || !is.finite(cr_floor)) cr_floor <- -Inf
      lookup <- tryCatch(
        pb_lookup(cr_floor),
        error = function(e) {
          msg <- conditionMessage(e)
          if (!inherits(e, "shiny.silent.error") || nzchar(msg)) {
            showNotification(
              paste0("Protein binding error: ", if (nzchar(msg)) msg else class(e)[1]),
              type = "error", duration = 10
            )
          }
          NULL
        }
      )
      removeNotification("load_swiss_pb", session)
      if (is.null(lookup) || !is.data.frame(lookup) || !nrow(lookup)) {
        pb_loaded_val(NULL)
        binding_ready(FALSE)
        showNotification(
          "Protein binding data is empty or missing. See docs/DATA_PROVENANCE.md.",
          type = "warning", duration = 10
        )
      } else {
        pb_loaded_val(lookup)
        binding_ready(TRUE)
        swiss_last_updated(Sys.time())
        showNotification("Protein binding loaded.", type = "message", duration = 4)
      }
    }, ignoreInit = TRUE)

    observeEvent(input$reset_swiss, {
      adme <- shared_data$swissadme_table()
      if (!is.null(adme) && nrow(adme) > 0 && "probe_name" %in% colnames(adme)) {
        probes <- gtools::mixedsort(unique(adme$probe_name))
        if (length(probes) > 0) {
          updateSelectInput(session, "probe_select", selected = probes[1])
        }
      }
      updateNumericInput(session, "cr_cutoff", value = 4)
      updateRadioButtons(session, "boiled_show", selected = "all")
      updateTabsetPanel(session, "swiss_tabs", selected = "Explorer")
      swiss_last_updated(NULL)
      binding_ready(FALSE)
      pb_loaded_val(NULL)
      pb_loaded_floor(Inf)
    })

    # Dropping the cutoff below what the loaded subset covers would silently
    # under-report, so fetch the full table before the lower cutoff is applied.
    observeEvent(input$cr_cutoff, {
      req(isTRUE(binding_ready()))
      cr <- suppressWarnings(as.numeric(input$cr_cutoff))
      if (!length(cr) || !is.finite(cr) || cr >= pb_loaded_floor()) return()
      showNotification("Loading the full binding table for this cutoff...",
                       type = "message", duration = NULL, id = "load_swiss_pb_full")
      full <- tryCatch(pb_full_lookup(), error = function(e) NULL)
      removeNotification("load_swiss_pb_full", session)
      if (is.data.frame(full) && nrow(full)) {
        pb_loaded_val(full)
        pb_loaded_floor(-Inf)
        swiss_last_updated(Sys.time())
      } else {
        showNotification("Could not load the full binding table.",
                         type = "warning", duration = 8)
      }
    }, ignoreInit = TRUE)

    output$swiss_status_line <- renderUI({
      stamp <- swiss_last_updated()
      stamp_txt <- if (is.null(stamp)) "not run yet" else format(stamp, "%Y-%m-%d %H:%M:%S")
      app_status_banner(
        tags$b("Context: "),
        "Probe = ", input$probe_select,
        " | Min CR = ", input$cr_cutoff,
        " | Last updated = ", stamp_txt,
        type = "info"
      )
    })

    related_probe_names <- reactive({
      adme <- shared_data$swissadme_table()
      req(adme, input$probe_select, "probe_name" %in% colnames(adme))
      req("Canonical.SMILES" %in% colnames(adme))

      adme2 <- adme %>%
        dplyr::mutate(
          probe_name_l = tolower(as.character(.data$probe_name)),
          smiles_l = tolower(trimws(as.character(.data$Canonical.SMILES)))
        )
      sel <- adme2 %>% dplyr::filter(.data$probe_name_l == tolower(input$probe_select))
      req(nrow(sel) > 0)
      smiles_sel <- sel$smiles_l[1]

      adme2 %>%
        dplyr::filter(.data$smiles_l == smiles_sel) %>%
        dplyr::pull(.data$probe_name_l) %>%
        unique()
    })

    pb_family_base <- reactive({
      req(isTRUE(binding_ready()))
      lookup <- pb_loaded_val()
      req(lookup, input$probe_select)
      rel <- related_probe_names()
      req(length(rel) > 0)
      indexed <- cpt_pb_subset(lookup, rel, "rows_by_probe")
      if (!is.null(indexed)) return(indexed)
      lookup %>%
        dplyr::filter(tolower(as.character(.data$probe_name)) %in% rel)
    })

    # Rows for selected probe family with CR cutoff applied.
    pb_data <- reactive({
      df <- pb_family_base()
      req(df)
      df %>% dplyr::filter(.data$CR >= input$cr_cutoff)
    })

    # All rows for that probe/SMILES across datasets (same as table)
    pb_all <- reactive({
      df <- pb_data()
      req(df)
      df %>%
        dplyr::select(
          dplyr::any_of(c(
            "probe_name", "proteinid", "gene_name", "CR", "n_targets",
            "cysteineid", "Dataset", "Cell_Line", "SMILES"
          ))
        )
    })

    pb_by_dataset <- reactive({
      pb_data() %>%
        dplyr::group_by(Dataset, gene_name) %>%
        dplyr::summarise(mean_CR = mean(CR, na.rm = TRUE), .groups = "drop") %>%
        dplyr::arrange(Dataset, desc(mean_CR))
    })

    output$pb_all <- renderDT({
      df <- pb_all()
      shiny::validate(need(nrow(df) > 0, "No protein-binding rows passed this CR cutoff for this probe family."))
      df %>%
        datatable(options = list(pageLength = 10, scrollX = TRUE),
                  rownames = FALSE)
    })

    pb_binding_plot_gg <- reactive({
      df <- pb_all()
      req(df, nrow(df) > 0)
      req("CR" %in% names(df), "Dataset" %in% names(df), "gene_name" %in% names(df))
      df <- df %>%
        dplyr::filter(is.finite(.data$CR)) %>%
        dplyr::mutate(
          Dataset = as.character(.data$Dataset),
          gene_name = as.character(.data$gene_name),
          probe_name = as.character(.data$probe_name),
          proteinid = as.character(.data$proteinid),
          Cell_Line = as.character(.data$Cell_Line)
        )
      req(nrow(df) > 0)

      target_rank <- df %>%
        dplyr::group_by(.data$Cell_Line, .data$gene_name) %>%
        dplyr::summarise(target_max_CR = max(.data$CR, na.rm = TRUE), .groups = "drop") %>%
        dplyr::group_by(.data$Cell_Line) %>%
        dplyr::arrange(dplyr::desc(.data$target_max_CR), .data$gene_name, .by_group = TRUE) %>%
        dplyr::mutate(target_rank = dplyr::row_number()) %>%
        dplyr::ungroup()

      df <- df %>%
        dplyr::left_join(target_rank, by = c("Cell_Line", "gene_name"))

      cl_sorted <- sort(unique(df$Cell_Line))
      pal_base <- c(
        "#45669C", "#B12000", "#7E4554", "#304A7E", "#64477C", "#986B76",
        "#65789F", "#8C769D", "#879A7A", "#385723", "#A44A4A", "#FCDC30"
      )
      n_cl <- length(cl_sorted)
      if (n_cl <= length(pal_base)) {
        cols <- pal_base[seq_len(n_cl)]
      } else {
        cols <- c(pal_base, grDevices::colorRampPalette(pal_base)(n_cl - length(pal_base)))
      }
      names(cols) <- cl_sorted
      df$Cell_Line <- factor(df$Cell_Line, levels = cl_sorted)
      df <- df %>%
        dplyr::arrange(.data$Cell_Line, .data$target_rank, dplyr::desc(.data$CR), .data$Dataset) %>%
        dplyr::mutate(
          tip = paste0(
            .data$gene_name,
            "\nTarget rank: ", .data$target_rank,
            "\nDataset: ", .data$Dataset,
            "\nCR: ", round(.data$CR, 3),
            "\nProbe: ", .data$probe_name,
            "\nProtein ID: ", .data$proteinid,
            "\nCell line: ", .data$Cell_Line
          )
        )
      ggplot2::ggplot(
        df,
        ggplot2::aes(
          x = .data$target_rank,
          y = .data$CR,
          colour = .data$Cell_Line,
          text = .data$tip
        )
      ) +
        ggplot2::geom_point(
          position = ggplot2::position_jitter(width = 0.2, height = 0, seed = 1),
          alpha = 0.88,
          size = 2.8
        ) +
        ggplot2::scale_colour_manual(values = cols, drop = FALSE) +
        ggplot2::labs(
          x = "Target Rank",
          y = "Competition Ratio (CR)",
          colour = "Cell Line",
          title = "Probe Protein Binding"
        ) +
        cpt_theme() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(face = "plain", size = 16, hjust = 0.5)
        )
    })

    output$pb_binding_scatter <- plotly::renderPlotly({
      req(binding_ready())
      shiny::validate(need(nrow(pb_all()) > 0, "No points to plot. Lower the CR cutoff and reload protein binding."))
      cpt_plotly(plotly::ggplotly(pb_binding_plot_gg(), tooltip = "text"))
    })

    output$pb_by_dataset <- renderDT({
      req(pb_by_dataset())
      pb_by_dataset() %>%
        datatable(options = list(pageLength = 10, scrollX = TRUE),
                  rownames = FALSE)
    })

    # ---- 2.7 Full SwissADME table download ---------------------
    full_adme_with_targets <- reactive({
      df <- shared_data$swissadme_table()
      if (!is.null(df) && "Molecule" %in% colnames(df)) {
        df <- df[, setdiff(colnames(df), "Molecule"), drop = FALSE]
      }
      order_swissadme_display_columns(df)
    })

    output$full_adme_table <- renderDT({
      input$reset_swiss
      df <- full_adme_with_targets()
      req(df)

      datatable(
        df,
        options = list(
          pageLength = 15,
          scrollX = TRUE,
          columnDefs = list(list(type = "probe_name_nat", targets = 0))
        ),
        rownames = FALSE,
        callback = DT::JS(
          "$.extend($.fn.dataTable.ext.type.order, {",
          "'probe_name_nat-pre': function(d) {",
          "  var s = (d === null || d === undefined) ? '' : d.toString().toUpperCase();",
          "  var m = s.match(/^([A-Z]+)(\\d+)$/);",
          "  if (m) { return m[1] + ('00000000' + m[2]).slice(-8); }",
          "  return s;",
          "},",
          "'probe_name_nat-asc': function(a,b) { return a < b ? -1 : (a > b ? 1 : 0); },",
          "'probe_name_nat-desc': function(a,b) { return a < b ? 1 : (a > b ? -1 : 0); }",
          "});"
        )
      )
    })

    output$dl_protein_binding_pb_tab <- downloadHandler(
      filename = function() paste0("protein_binding_",
        input$probe_select, "_",
        Sys.Date(), ".xlsx"),
      content  = function(file) {
        write_xlsx(list(all = pb_all()), file)
      }
    )

    output$dl_pb_binding_plot <- downloadHandler(
      filename = function() {
        paste0("probe_protein_binding_", input$probe_select, "_", Sys.Date(), ".png")
      },
      content = function(file) {
        req(binding_ready(), nrow(pb_all()) > 0)
        ggplot2::ggsave(
          file,
          plot = pb_binding_plot_gg(),
          width = 10,
          height = 6,
          dpi = 300
        )
      }
    )

    output$dl_swissadme_full_tab <- downloadHandler(
      filename = function() paste0("swissadme_full_", Sys.Date(), ".xlsx"),
      content  = function(file) write_xlsx(full_adme_with_targets(), file)
    )

  }) # end moduleServer
} # end swissadme_server

