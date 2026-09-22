# ============================================================
# Script:   gene_module.R
# Purpose:  Target tab — everything known about ONE target, in one view.
#
# Why this tab exists
#   Everything known about one gene, on one page, without running a subtype
#   analysis first: where it is a dependency, how its effect is distributed
#   across cell lines, which of its cysteines are engaged and at what evidence
#   tier, and the ligands that engage them.
#
#   This reads data/gene_index.rds (built offline by
#   docs/scripts/build_gene_index.R) so every section answers immediately and
#   nothing here touches the 140 MB matrix or the 10.6M-row binding table.
# ============================================================
# Sections:
#   cpt_tier_badge() — Tier badge HTML for a tier integer
#   gene_ui()
#   gene_server()
#   cpt_gene_card()
# ============================================================

# Paper Tiers 1-4, as defined by cpt_evidence_tier() in R/canprotarget_score.R.
# These labels must not drift from that function or from the manuscript: the
# tier is a property of an ENGAGED residue (CR >= 4), not of the gene.
CPT_TIER_LABELS <- c(
  "1" = "Functional, with atlas ligandability",
  "2" = "Functional, without atlas ligandability",
  "3" = "Tested, not functional in assayed contexts",
  "4" = "Absent from atlas (untested)"
)

CPT_TIER_HELP <- paste(
  "Tiers are assigned to each cysteine engaged by a small-molecule covalent",
  "ligand (competition ratio >= 4), then the gene takes its best engaged site.",
  "Tier 1 is an exact functional site that also carries an independent atlas",
  "ligandability annotation; Tier 2 is functional without that annotation;",
  "Tier 3 was tested but did not meet the functional criteria in the assayed",
  "Cys_editing contexts; Tier 4 is absent from the atlas and therefore untested.",
  "Functional evidence elsewhere in the same protein is shown but never",
  "transferred to the engaged residue — that transfer is the error the tiers",
  "exist to prevent.",
  "\n\nNote that two different measurements are both called ligandability.",
  "Competition ratio is engagement measured by a chemoproteomic probe; atlas",
  "ligandability is the Cys_editing atlas's own annotation for the residue, and",
  "it is what separates Tier 1 from Tier 2. A site engaged at CR 5.8 with atlas",
  "ligandability FALSE is a Tier 2 site, not a contradiction."
)

#' Tier badge HTML for a tier integer.
cpt_tier_badge <- function(tier) {
  if (is.null(tier) || is.na(tier)) {
    return(shiny::tags$span(class = "cpt-tier cpt-tier-na", "No cysteine data"))
  }
  shiny::tags$span(
    class = paste0("cpt-tier cpt-tier-", tier),
    paste0("Tier ", tier, " · ", CPT_TIER_LABELS[[as.character(tier)]])
  )
}

gene_ui <- function(id) {
  ns <- NS(id)

  tabItem(
    tabName = "gene_tab",
    fluidRow(
      box(
        title = "Target", width = 3, status = "primary", solidHeader = TRUE,
        cpt_label("Gene", for_id = ns("gene")),
        selectizeInput(
          ns("gene"), label = NULL, choices = character(0),
          options = list(placeholder = "Search a gene symbol…", maxOptions = 1000)
        ),
        cpt_label("Gene effect dataset", for_id = ns("dataset")),
        selectInput(ns("dataset"), label = NULL,
                    choices = c("RNAi" = "RNAi", "CRISPR (23Q4)" = "CRISPR (23Q4)"),
                    selected = "RNAi"),
        hr(),
        cpt_label("Report subtype", for_id = ns("report_subtype")),
        selectInput(ns("report_subtype"), label = NULL, choices = character(0)),
        downloadButton(ns("dl_target_report"), "Download report (HTML)",
                       class = "btn-primary btn-block btn-sm", icon = icon("file-lines")),
        uiOutput(ns("jump_ui"))
      ),
      column(
        9,
        uiOutput(ns("gene_header")),
        fluidRow(
          valueBoxOutput(ns("vb_dependency"), width = 3),
          valueBoxOutput(ns("vb_probes"), width = 3),
          valueBoxOutput(ns("vb_cysteines"), width = 3),
          valueBoxOutput(ns("vb_clinvar"), width = 3)
        ),
        box(
          width = 12, status = "info", solidHeader = FALSE,
          cpt_section("Where this gene is a dependency"),
          shinycssloaders::withSpinner(
            DTOutput(ns("dependency_table")), type = 6, color = "#478EB8"
          ),
          br(),
          downloadButton(ns("dl_dep_table"), "Download table (.xlsx)",
                         class = "btn-sm btn-default")
        ),
        box(
          width = 12, status = "info", solidHeader = FALSE,
          cpt_section("Gene effect across subtypes"),
          selectizeInput(ns("box_subtypes"), "Subtypes", choices = character(0),
                         multiple = TRUE,
                         options = list(placeholder = "Select subtypes to compare…")),
          shinycssloaders::withSpinner(
            shinycssloaders::withSpinner(
              plotlyOutput(ns("effect_box"), height = "420px"),
              type = 6, color = "#478EB8"),
            type = 6, color = "#478EB8"
          ),
          br(),
          downloadButton(ns("dl_effect_box"), "Download PNG",
                         class = "btn-sm btn-default"),
          downloadButton(ns("dl_effect_data"), "Download plot data (.xlsx)",
                         class = "btn-sm btn-default"),
          tags$details(
            class = "cpt-collapse-block",
            tags$summary(class = "cpt-collapse-summary",
                         icon("circle-info"), " Data behind these boxes"),
            tags$div(
              class = "cpt-collapse-body",
              tabsetPanel(
                id = ns("effect_data_tabs"),
                tabPanel(
                  "Individual cell lines",
                  br(),
                  tags$p(class = "text-muted", style = "font-size: 12px;",
                         "Every line behind every box, not only the subtype's: ",
                         "the Group column says which box each line sits in. ",
                         "Read from the gene effect matrix, so the first gene you ",
                         "open in a session takes a few seconds."),
                  DTOutput(ns("effect_points_table")),
                  br(),
                  downloadButton(ns("dl_effect_points"),
                                 "Download cell lines (.xlsx)",
                                 class = "btn-sm btn-default")
                ),
                tabPanel(
                  "Group averages",
                  br(),
                  tags$p(class = "text-muted", style = "font-size: 12px;",
                         "Summary of the same lines listed beside this tab, so ",
                         "the averages cannot drift from the values they cover."),
                  DTOutput(ns("effect_reference_table"))
                )
              )
            )
          ),
        ),
        box(
          width = 12, status = "info", solidHeader = FALSE,
          cpt_section("Compare two subtypes"),
          fluidRow(
            column(5, selectInput(ns("cmp_a"), "Subtype A", choices = character(0))),
            column(5, selectInput(ns("cmp_b"), "Subtype B", choices = character(0))),
            column(2, br(), actionButton(ns("cmp_run"), "Compare",
                                         class = "btn-primary btn-sm btn-block"))
          ),
          uiOutput(ns("cmp_result")),
          downloadButton(ns("dl_cmp"), "Download comparison (.xlsx)",
                         class = "btn-sm btn-default")
        ),
        box(
          width = 12, status = "info", solidHeader = FALSE,
          cpt_section("Functional cysteines"),
          uiOutput(ns("engaged_caption")),
          uiOutput(ns("tier_summary")),
          br(),
          shinycssloaders::withSpinner(
            DTOutput(ns("cys_table")), type = 6, color = "#478EB8"
          ),
          br(),
          downloadButton(ns("dl_cys_table"), "Download table (.xlsx)",
                         class = "btn-sm btn-default")
        ),
        box(
          width = 12, status = "info", solidHeader = FALSE,
          cpt_section("Covalent probes"),
          shinycssloaders::withSpinner(
            DTOutput(ns("probe_table")), type = 6, color = "#478EB8"
          ),
          br(),
          downloadButton(ns("dl_probe_table"), "Download table (.xlsx)",
                         class = "btn-sm btn-default"),
          downloadButton(ns("dl_gene"), "Download everything for this gene (.xlsx)",
                         class = "btn-sm btn-default")
        )
      )
    )
  )
}

#' @param preset Optional reactive carrying a gene key to select (from the
#'   Discover drawer's "Open in Gene tab").
#' @param dep_effect_min Optional reactive: the effect-size cutoff set on the
#'   Discover tab, so both tabs agree on what counts as a dependency.
#' @param dep_context Optional list of reactives from Discover: dataset,
#'   subtype and a one-line context label. The Target tab follows them instead
#'   of keeping its own dataset control.
gene_server <- function(id, shared_data, on_jump = NULL, preset = NULL,
                        dep_effect_min = NULL, dep_context = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # One dataset, two controls. Discover remains the source of truth; this
    # dropdown writes to it and mirrors it back, with identity guards so the
    # two updates cannot chase each other.
    # selectize empties its value while the user types into it, and every output
    # keyed on input$gene went blank mid-search. Hold the last real selection
    # and render from that, so the page only changes once a gene is chosen.
    gene_held <- reactiveVal(NULL)
    observeEvent(input$gene, {
      g <- input$gene
      if (!is.null(g) && nzchar(g)) gene_held(g)
    }, ignoreNULL = FALSE)
    gene_sel <- reactive(gene_held())

    active_dataset <- reactive({
      if (!is.null(dep_context) && is.function(dep_context$dataset)) {
        dep_context$dataset()
      } else {
        if (identical(input$dataset, "CRISPR (23Q4)")) "CRISPR" else "RNAi"
      }
    })

    observeEvent(input$dataset, {
      if (!is.null(dep_context) && is.function(dep_context$set_dataset)) {
        dep_context$set_dataset(input$dataset)
      }
    }, ignoreInit = TRUE)

    observe({
      if (is.null(dep_context) || !is.function(dep_context$dataset)) return()
      want <- if (identical(dep_context$dataset(), "CRISPR")) "CRISPR (23Q4)" else "RNAi"
      if (!identical(want, isolate(input$dataset))) {
        updateSelectInput(session, "dataset", selected = want)
      }
    })

    index <- reactive({
      idx <- shared_data$gene_index()
      validate(need(
        !is.null(idx),
        paste("Gene index not available. Build it with",
              "Rscript docs/scripts/build_gene_index.R")
      ))
      idx
    })

    observeEvent(shared_data$gene_index(), {
      idx <- shared_data$gene_index()
      if (is.null(idx) || !nrow(idx$genes)) return()
      updateSelectizeInput(
        session, "gene",
        choices = stats::setNames(idx$genes$gene_key, idx$genes$symbol),
        server = TRUE
      )
    }, ignoreInit = FALSE)

    if (!is.null(preset)) {
      observeEvent(preset(), {
        key <- preset()
        req(key, nzchar(key))
        idx <- shared_data$gene_index()
        req(idx)
        # Server-side selectize needs the choices resent alongside the selection.
        updateSelectizeInput(
          session, "gene",
          choices = stats::setNames(idx$genes$gene_key, idx$genes$symbol),
          selected = key, server = TRUE
        )
      }, ignoreInit = TRUE)
    }

    gene_row <- reactive({
      idx <- index()
      key <- gene_sel()
      req(key, nzchar(key))
      row <- idx$genes[idx$genes$gene_key == key, , drop = FALSE]
      if (!nrow(row)) return(NULL)
      row[1, , drop = FALSE]
    })

    output$gene_header <- renderUI({
      row <- gene_row()
      if (is.null(row)) {
        return(cpt_empty_state("Pick a gene to see its dependency, cysteine and probe evidence."))
      }
      tags$div(
        class = "cpt-gene-header",
        tags$h2(row$symbol),
        cpt_tier_badge(row$best_engaged_tier)
      )
    })

    output$vb_dependency <- renderValueBox({
      row <- gene_row()
      n <- if (is.null(row) || is.na(row$dep_n_subtypes)) 0 else row$dep_n_subtypes
      valueBox(
        n, "Subtypes where dependent",
        icon = icon("dna"), color = if (n > 0) "navy" else "light-blue"
      )
    })

    # Distinct probes engaging this gene at CR >= 4. n_probes in the index
    # counts every record regardless of CR, which contradicted the probe table
    # below (FLT3 read "1 covalent probe" while its only record is CR 0.93).
    n_engaged_probes <- reactive({
      idx <- index()
      key <- gene_sel()
      req(key, nzchar(key))
      es <- idx$engaged_sites
      if (is.null(es) || !nrow(es)) return(0L)
      keys <- cpt_gene_match_key(es$gene_name)
      length(unique(es$probe_name[keys == key]))
    })

    output$vb_probes <- renderValueBox({
      n <- tryCatch(n_engaged_probes(), error = function(e) 0L)
      row <- tryCatch(gene_row(), error = function(e) NULL)
      total <- if (is.null(row) || is.na(row$n_probes)) 0L else as.integer(row$n_probes)
      valueBox(
        n,
        if (total > n) {
          paste0("Probes at CR \u2265 4 (", format(total, big.mark = ","), " records total)")
        } else {
          "Probes at CR \u2265 4"
        },
        icon = icon("flask"), color = "light-blue"
      )
    })

    output$vb_cysteines <- renderValueBox({
      row <- gene_row()
      nf <- if (is.null(row) || is.na(row$n_functional)) 0 else row$n_functional
      ns_ <- if (is.null(row) || is.na(row$n_sites)) 0 else row$n_sites
      valueBox(
        paste0(nf, " / ", ns_), "Functional cysteine sites",
        icon = icon("bullseye"), color = "light-blue"
      )
    })

    output$vb_clinvar <- renderValueBox({
      row <- gene_row()
      n <- if (is.null(row) || is.na(row$n_clinvar)) 0 else row$n_clinvar
      valueBox(
        n, "ClinVar pathogenic sites",
        icon = icon("triangle-exclamation"), color = "light-blue"
      )
    })

    gene_dependency <- reactive({
      idx <- index()
      key <- gene_sel()
      req(key, nzchar(key))
      out <- idx$dependency %>%
        dplyr::filter(.data$gene_key == key, .data$dataset == active_dataset()) %>%
        dplyr::arrange(.data$effect_size) %>%
        dplyr::select(
          Subtype = "subtype", `Effect size` = "effect_size", `p` = "p_value",
          `Avg (all lines)` = "avg", `Cancer avg` = "cancer_avg",
          `Non-cancer avg` = "noncancer_avg", `p vs non-cancer` = "pval_vs_noncancer"
        )

      # Per-line values stay available in the cell-line table under the plot,
      # which lists every screened line and its gene effect.
      out
    })

    output$dependency_table <- renderDT({
      df <- gene_dependency()
      # emptyTable rather than validate(): a validation error inside renderDT
      # leaves the previous gene's rows on screen.
      datatable(
        df %>% dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 4))),
        rownames = FALSE, selection = "single",
        options = list(
          pageLength = 10, scrollX = TRUE,
          language = list(emptyTable = sprintf(
            "This gene is not a cancer-selective dependency in %s data.",
            active_dataset()))
        )
      )
    })

    # Subtype choices follow the gene and dataset; default to the three
    # strongest dependencies so the plot is populated on arrival.
    observeEvent(list(input$gene, active_dataset()), {
      df <- tryCatch(gene_dependency(), error = function(e) NULL)
      dep <- if (is.null(df) || !nrow(df)) character(0) else df$Subtype

      # effect_distributions stores "Cancer of Interest" rows only for gene x
      # subtype pairs that clear the dependency filter, 86,406 of a possible
      # 3,226,600, so offering only those would limit the plot to subtypes the
      # gene is already a dependency in. Every subtype in the dataset is
      # offered instead; a box the index does not hold is computed from the
      # gene-effect matrix when it is picked.
      st_file <- tryCatch(
        shared_data$cancer_subtype_list_files[[
          if (identical(active_dataset(), "CRISPR")) "CRISPR (23Q4)" else "RNAi"]],
        error = function(e) NULL)
      all_st <- character(0)
      if (!is.null(st_file)) {
        pth <- shared_data$data_path(st_file)
        if (file.exists(pth)) all_st <- read_cancer_subtype_list_file(pth)
      }
      if (!length(all_st)) {
        idx <- tryCatch(index(), error = function(e) NULL)
        if (!is.null(idx) && !is.null(idx$effect_distributions)) {
          d <- idx$effect_distributions
          d <- d[d$gene_key == gene_sel() & d$dataset == active_dataset() &
                   d$group == "Cancer of Interest", , drop = FALSE]
          all_st <- sort(unique(d$subtype[!is.na(d$subtype)]))
        }
      }
      others <- setdiff(all_st, dep)
      choices <- if (length(others)) {
        list(`Dependency in this subtype` = as.list(dep),
             `Other subtypes` = as.list(others))
      } else if (length(dep)) as.list(dep) else character(0)

      # server = FALSE: grouped choices are not supported in server mode, and
      # the list is small enough that client-side filtering is instant.
      updateSelectizeInput(session, "box_subtypes", choices = choices,
                           selected = utils::head(dep, 3), server = FALSE)
    }, ignoreInit = FALSE)

    observeEvent(list(input$gene, active_dataset()), {
      df <- tryCatch(gene_dependency(), error = function(e) NULL)
      ch <- if (is.null(df) || !nrow(df)) character(0) else df$Subtype
      # Any subtype can be compared, not only those with a stored distribution:
      # the per-cell-line values come from the gene effect matrix, so a subtype
      # that never cleared the dependency filter still has lines to summarise.
      st_file <- tryCatch(
        shared_data$cancer_subtype_list_files[[
          if (identical(active_dataset(), "CRISPR")) "CRISPR (23Q4)" else "RNAi"]],
        error = function(e) NULL)
      all_st <- character(0)
      if (!is.null(st_file)) {
        pth <- shared_data$data_path(st_file)
        if (file.exists(pth)) all_st <- read_cancer_subtype_list_file(pth)
      }
      others <- setdiff(all_st, ch)
      grouped <- if (length(others)) {
        list(`Meets the dependency filter` = as.list(ch),
             `Every other subtype in this dataset` = as.list(others),
             `Reference groups` = list("Non-cancer lines", "All lines"))
      } else {
        list(Subtypes = as.list(ch),
             `Reference groups` = list("Non-cancer lines", "All lines"))
      }
      updateSelectInput(session, "cmp_a", choices = grouped,
                        selected = if (length(ch) >= 1) ch[[1]] else "All lines")
      updateSelectInput(session, "cmp_b", choices = grouped,
                        selected = if (length(ch) >= 2) ch[[2]] else "Non-cancer lines")
    }, ignoreInit = FALSE)

    cmp_out <- eventReactive(input$cmp_run, {
      idx <- index()
      key <- gene_sel()
      req(key, nzchar(key), input$cmp_a, input$cmp_b)
      req(!identical(input$cmp_a, input$cmp_b))

      dep <- idx$dependency
      dist <- idx$effect_distributions
      # "Non-cancer lines" and "All lines" are stored with group set and no
      # subtype; a real subtype is stored under "Cancer of Interest". Neither
      # reference group has a limma effect size, because neither is a cohort
      # tested against the rest.
      ref_group <- c("Non-cancer lines" = "Non-Cancer", "All lines" = "All lines")
      grab <- function(st) {
        is_ref <- st %in% names(ref_group)
        d <- if (is_ref) dep[0, , drop = FALSE] else
          dep[dep$gene_key == key & dep$dataset == active_dataset() &
                dep$subtype == st, , drop = FALSE]
        x <- if (is_ref) {
          dist[dist$gene_key == key & dist$dataset == active_dataset() &
                 dist$group == unname(ref_group[[st]]), , drop = FALSE]
        } else {
          dist[dist$gene_key == key & dist$dataset == active_dataset() &
                 dist$group == "Cancer of Interest" & dist$subtype == st, , drop = FALSE]
        }
        n <- if (nrow(x)) x$n[1] else NA_integer_
        mu <- if (nrow(x)) x$mean[1] else NA_real_
        sdv <- if (nrow(x)) x$sd[1] else NA_real_
        # Subtypes outside the dependency results have no stored distribution,
        # so summarise their cell lines from the matrix instead.
        if (is.na(mu) && !is_ref) {
          pts <- all_line_points()
          if (!is.null(pts)) {
            v <- pts$`Gene effect`[!is.na(pts$Subtype) & pts$Subtype == st]
            if (length(v)) { n <- length(v); mu <- mean(v); sdv <- stats::sd(v) }
          }
        }
        list(
          subtype = st, n = n, mean = mu, sd = sdv,
          effect = if (nrow(d)) d$effect_size[1] else NA_real_,
          p = if (nrow(d)) d$p_value[1] else NA_real_
        )
      }
      a <- grab(input$cmp_a); b <- grab(input$cmp_b)
      shiny::validate(need(!is.na(a$mean) || !is.na(b$mean),
                           "Neither subtype has screened cell lines for this gene."))
      list(a = a, b = b)
    })

    output$cmp_result <- renderUI({
      r <- cmp_out()
      fmt <- function(x, d = 3) if (is.null(x) || is.na(x)) "\u2014" else format(round(x, d), nsmall = d)
      fmtp <- function(x) {
        if (is.null(x) || is.na(x)) return("\u2014")
        if (x < 0.001) format(x, scientific = TRUE, digits = 3) else format(round(x, 4), nsmall = 4)
      }
      cell <- function(g) {
        tags$div(
          class = "cpt-cmp-col",
          tags$div(class = "cpt-cmp-name", g$subtype),
          tags$div(class = "cpt-cmp-row", tags$span("cell lines"), tags$b(g$n)),
          tags$div(class = "cpt-cmp-row", tags$span("mean effect"), tags$b(fmt(g$mean))),
          tags$div(class = "cpt-cmp-row", tags$span("sd"), tags$b(fmt(g$sd))),
          tags$div(class = "cpt-cmp-row cpt-cmp-stat", tags$span("limma effect size"), tags$b(fmt(g$effect))),
          tags$div(class = "cpt-cmp-row cpt-cmp-stat", tags$span("limma p (vs all others)"), tags$b(fmtp(g$p)))
        )
      }
      # Difference with a 95% interval: answers how far apart and how sure,
      # without inviting the "significant vs non-significant" misreading that
      # two independent p-values side by side always attract.
      se <- sqrt(r$a$sd^2 / r$a$n + r$b$sd^2 / r$b$n)
      dfree <- (r$a$sd^2 / r$a$n + r$b$sd^2 / r$b$n)^2 /
        ((r$a$sd^2 / r$a$n)^2 / (r$a$n - 1) + (r$b$sd^2 / r$b$n)^2 / (r$b$n - 1))
      diff <- r$a$mean - r$b$mean
      tcrit <- stats::qt(0.975, df = dfree)
      tagList(
        tags$div(class = "cpt-cmp-grid", cell(r$a), cell(r$b)),
        tags$div(
          class = "cpt-cmp-diff",
          tags$span(class = "cpt-cmp-diff-label", "Difference in mean gene effect"),
          tags$span(class = "cpt-cmp-diff-value",
                    sprintf("%s  (95%% CI %s to %s)", fmt(diff),
                            fmt(diff - tcrit * se), fmt(diff + tcrit * se)))
        ),
        tags$p(class = "text-muted", style = "font-size: 11.5px; margin-top: 6px;",
               "Interval from the stored per-subtype means, standard deviations ",
               "and counts. It describes the gap between these two subtypes; the ",
               "p-values above each test their own subtype against all remaining ",
               "models.")
      )
    })

    effect_box_df <- reactive({
      idx <- index()
      key <- gene_sel()
      req(key, nzchar(key))
      sel <- input$box_subtypes
      dist <- idx$effect_distributions
      shiny::validate(need(!is.null(dist), paste(
        "This gene index predates the distribution data. Rebuild it with",
        "Rscript docs/scripts/build_gene_index.R"
      )))
      dist <- dist[dist$gene_key == key & dist$dataset == active_dataset(), , drop = FALSE]
      shiny::validate(need(nrow(dist) > 0, "No effect distribution for this gene."))
      shiny::validate(need(length(sel) > 0, "Select at least one subtype."))

      # Guarded like the reference groups below. Switching to a gene with no
      # rows for the subtype still selected leaves this empty, and assigning a
      # single string to a column of a zero-row frame throws.
      of_interest <- dist[dist$group == "Cancer of Interest" &
                            dist$subtype %in% sel, , drop = FALSE]

      # The index only stores a box for a subtype the gene is a dependency in,
      # so any other pick would silently draw nothing. Compute those from the
      # gene-effect matrix, on the same five-number summary the index uses
      # (build_gene_index.R: min, quartiles, max).
      missing_st <- setdiff(sel, of_interest$subtype)
      if (length(missing_st)) {
        pts <- tryCatch(all_line_points(), error = function(e) NULL)
        if (!is.null(pts) && nrow(pts)) {
          extra <- lapply(missing_st, function(st) {
            v <- pts$`Gene effect`[!is.na(pts$Subtype) & pts$Subtype == st]
            v <- v[!is.na(v)]
            if (!length(v)) return(NULL)
            q <- stats::quantile(v, c(0, 0.25, 0.5, 0.75, 1), names = FALSE)
            data.frame(gene_key = key, dataset = active_dataset(), subtype = st,
                       group = "Cancer of Interest", n = length(v),
                       ymin = q[1], lower = q[2], middle = q[3], upper = q[4],
                       ymax = q[5], mean = mean(v), sd = stats::sd(v),
                       stringsAsFactors = FALSE)
          })
          extra <- do.call(rbind, Filter(Negate(is.null), extra))
          if (!is.null(extra) && nrow(extra)) {
            of_interest <- dplyr::bind_rows(of_interest, extra)
          }
        }
      }
      if (nrow(of_interest)) {
        of_interest$label <- of_interest$subtype
        of_interest$role <- "Subtype"
      }

      noncancer <- dist[dist$group == "Non-Cancer", , drop = FALSE]
      if (nrow(noncancer)) {
        noncancer <- noncancer[1, , drop = FALSE]
        noncancer$label <- "Non-cancer lines"
        noncancer$role <- "Reference"
      }

      # Whole-matrix reference: the population behind the common-essential
      # filter, so a subtype box can be read against every screened line.
      all_lines <- dist[dist$group == "All lines", , drop = FALSE]
      if (nrow(all_lines)) {
        all_lines <- all_lines[1, , drop = FALSE]
        all_lines$label <- "All lines"
        all_lines$role <- "Reference"
      }

      other <- dist[0, , drop = FALSE]
      if (length(sel) == 1) {
        other <- dist[dist$group == "Other Cancers" & dist$subtype == sel[[1]], , drop = FALSE]
        if (nrow(other)) {
          other <- other[1, , drop = FALSE]
          other$label <- "Other cancer lines"
          other$role <- "Reference"
        }
      }

      out <- dplyr::bind_rows(of_interest, other, noncancer, all_lines)
      shiny::validate(need(nrow(out) > 0, "No distribution rows for that selection."))
      # Horizontal layout below, so the name and n sit on one line; reverse the
      # levels because coord_flip() draws the first level at the bottom.
      out$label <- paste0(out$label, "  (n=", out$n, ")")
      out$label <- factor(out$label, levels = rev(out$label))
      out
    })

    # Individual cell lines behind the subtype boxes. Reference groups keep no
    # points: they run to ~1,100 lines, which is not readable as jitter.
    effect_points_df <- reactive({
      idx <- index()
      boxes <- effect_box_df()
      pts <- idx$effect_points
      if (is.null(pts) || !nrow(pts)) return(NULL)
      key <- gene_sel()
      p <- pts[pts$gene_key == key & pts$dataset == active_dataset() &
                 pts$subtype %in% boxes$subtype, , drop = FALSE]
      if (!nrow(p)) return(NULL)
      p$label <- boxes$label[match(p$subtype, boxes$subtype)]
      p <- p[!is.na(p$label), , drop = FALSE]
      if (!nrow(p)) return(NULL)
      p$label <- factor(as.character(p$label), levels = levels(boxes$label))
      p
    })

    effect_box_gg <- reactive({
      df <- effect_box_df()
      pts <- tryCatch(effect_points_df(), error = function(e) NULL)
      ggplot2::ggplot(df, ggplot2::aes(x = .data$label, fill = .data$role)) +
        ggplot2::geom_boxplot(
          ggplot2::aes(ymin = .data$ymin, lower = .data$lower, middle = .data$middle,
                       upper = .data$upper, ymax = .data$ymax),
          stat = "identity", width = 0.6, colour = CPT_PAL$ink
        ) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dotted", colour = "grey55") +
        (if (!is.null(pts)) {
          ggplot2::geom_jitter(
            data = pts,
            mapping = ggplot2::aes(x = .data$label, y = .data$value),
            inherit.aes = FALSE, width = 0.16, height = 0,
            colour = CPT_PAL$ink, alpha = 0.55, size = 1.5
          )
        } else NULL) +
        ggplot2::scale_fill_manual(values = c(Subtype = CPT_PAL$primary, Reference = CPT_PAL$muted)) +
        ggplot2::labs(x = NULL, y = "Gene effect", fill = NULL) +
        # Subtype names are long ("Diffuse Large B-Cell Lymphoma, NOS"); read
        # them down the side rather than colliding along the x-axis.
        ggplot2::coord_flip() +
        cpt_theme() +
        ggplot2::theme(legend.position = "none",
                       panel.grid.major.y = ggplot2::element_blank())
    })

  # The static version could not answer "which line is that point?", which is
  # the first thing anyone asks of a dependency distribution. Boxes are still
  # drawn from the stored quantiles, so the shapes match the deposited numbers
  # exactly rather than being recomputed from the points.
    output$effect_box <- plotly::renderPlotly({
      df <- effect_box_df()
      shiny::validate(need(!is.null(df) && nrow(df) > 0,
                           "Choose at least one subtype."))
      pts <- tryCatch(effect_points_df(), error = function(e) NULL)
      meta <- tryCatch(shared_data$cancer_model_data(), error = function(e) NULL)

      labs <- levels(df$label)
      df$pos <- match(as.character(df$label), labs)

      p <- plotly::plot_ly()
      for (k in seq_len(nrow(df))) {
        r <- df[k, ]
        p <- plotly::add_trace(
          p, type = "box", orientation = "h",
          y = list(r$pos), q1 = list(r$lower), median = list(r$middle),
          q3 = list(r$upper), lowerfence = list(r$ymin), upperfence = list(r$ymax),
          name = as.character(r$label), showlegend = FALSE,
          fillcolor = if (identical(as.character(r$role), "Subtype"))
            CPT_PAL$primary else CPT_PAL$muted,
          line = list(color = CPT_PAL$ink, width = 1.4),
          # A precomputed box emits a hover label for EVERY statistic, so
          # brushing a whisker papered the plot with five rotated labels.
          # hoveron="boxes" collapses them to one label for the whole box.
          hoveron = "boxes", hoverinfo = "text",
          text = paste0(
            as.character(r$label), " (n=", r$n, ")",
            " — median ", signif(r$middle, 3),
            ", quartiles ", signif(r$lower, 3), " to ", signif(r$upper, 3),
            ", whiskers ", signif(r$ymin, 3), " to ", signif(r$ymax, 3)))
      }

      if (!is.null(pts) && nrow(pts)) {
        pts$pos <- match(as.character(pts$label), labs)
        pts <- pts[!is.na(pts$pos), , drop = FALSE]
        set.seed(2)
        jit <- stats::runif(nrow(pts), -0.17, 0.17)
        extra <- rep("", nrow(pts))
        if (!is.null(meta)) {
          m <- match(pts$cell_line, meta$StrippedCellLineName)
          nm <- ifelse(is.na(m), NA, meta$CellLineName[m])
          lin <- ifelse(is.na(m), NA, meta$OncotreeLineage[m])
          extra <- paste0(
            ifelse(is.na(nm) | nm == pts$cell_line, "", paste0("<br>", nm)),
            ifelse(is.na(lin), "", paste0("<br>", lin)))
        }
        p <- plotly::add_trace(
          p, type = "scattergl", mode = "markers",
          x = pts$value, y = pts$pos + jit, showlegend = FALSE,
          marker = list(size = CPT_MARKER_SIZE, color = CPT_PAL$ink,
                        opacity = 0.65, line = list(width = 0)),
          text = paste0("<b>", pts$cell_line, "</b>", extra,
                        "<br>", pts$subtype,
                        "<br>gene effect ", signif(pts$value, 3)),
          hoverinfo = "text")
      }

      cpt_plotly(
        p,
        xaxis = cpt_axis("Gene effect"),
        yaxis = cpt_axis("", tickmode = "array", tickvals = seq_along(labs),
                         ticktext = labs, range = c(0.4, length(labs) + 0.6)),
        shapes = list(list(type = "line", x0 = 0, x1 = 0, y0 = 0, y1 = 1,
                           yref = "paper",
                           line = list(color = "#b0bcc6", dash = "dot", width = 1))),
        margin = list(t = 20, l = 170))
    })

    output$effect_reference_table <- renderDT({
      df <- all_line_summary()
      if (is.null(df)) {
        # Fall back to the stored quantiles when the matrix is not available.
        idx <- index()
        key <- gene_sel()
        req(key, nzchar(key))
        d <- idx$effect_distributions
        ref <- d[d$gene_key == key & d$dataset == active_dataset() &
                   d$group %in% c("Non-Cancer", "All lines", "Other Cancers"), , drop = FALSE]
        shiny::validate(need(nrow(ref) > 0, "No reference distributions for this gene."))
        df <- ref[, c("group", "n", "ymin", "lower", "middle", "upper",
                      "ymax", "mean", "sd"), drop = FALSE]
        names(df) <- c("Group", "n", "Min", "Q1", "Median", "Q3", "Max", "Mean", "SD")
        df[3:9] <- lapply(df[3:9], function(x) round(x, 4))
      }
      datatable(cpt_filter_levels(df), rownames = FALSE, filter = "top",
                options = list(pageLength = 10, scrollX = TRUE))
    })

    # Every line behind every box, not just the subtype's.
    #
    # The index stores per-cell-line points for the subtype group only: the
    # reference groups run to ~1,100 lines per gene, which is not worth holding
    # for all 18k genes. So the values come from the gene-effect matrix, read
    # once per dataset and cached for the session by Discover. First use costs a
    # few seconds; afterwards it is free.
    all_line_points <- reactive({
      key <- gene_sel()
      req(key, nzchar(key))
      getter <- if (!is.null(dep_context)) dep_context$ge_matrix_for else NULL
      if (!is.function(getter)) return(NULL)
      label <- if (identical(active_dataset(), "CRISPR")) "CRISPR (23Q4)" else "RNAi"
      mat <- tryCatch(getter(label), error = function(e) NULL)
      if (is.null(mat)) return(NULL)
      if (!is.matrix(mat)) mat <- as.matrix(mat)
      j <- match(toupper(key), toupper(colnames(mat)))
      if (is.na(j)) return(NULL)
      vals <- mat[, j]
      meta <- tryCatch(shared_data$cancer_model_data(), error = function(e) NULL)
      if (is.null(meta)) return(NULL)
      m <- match(rownames(mat), meta$ModelID)
      subtype_of <- as.character(meta$OncotreeSubtype[m])
      disease_of <- as.character(meta$OncotreePrimaryDisease[m])
      # The same three groups the boxes are drawn from, assigned the same way as
      # in docs/scripts/build_gene_index.R so the table and the plot agree.
      focus <- as.character(input$box_subtypes %||% character(0))
      grp <- ifelse(!is.na(subtype_of) & subtype_of %in% focus, "Cancer of Interest",
                    ifelse(!is.na(disease_of) & disease_of == "Non-Cancerous",
                           "Non-Cancer", "Other Cancers"))
      out <- data.frame(
        Group = grp,
        Subtype = subtype_of,
        `Cell line` = if ("StrippedCellLineName" %in% names(meta))
          as.character(meta$StrippedCellLineName[m]) else rownames(mat),
        `Gene effect` = round(as.numeric(vals), 4),
        `Model ID` = rownames(mat),
        check.names = FALSE, stringsAsFactors = FALSE
      )
      out <- out[!is.na(out$`Gene effect`), , drop = FALSE]
      if (!nrow(out)) return(NULL)
      out[order(factor(out$Group, levels = c("Cancer of Interest", "Other Cancers",
                                             "Non-Cancer")),
                out$`Gene effect`), , drop = FALSE]
    })

    # Group averages over exactly the rows above, so the summary cannot drift
    # from the lines it summarises.
    all_line_summary <- reactive({
      d <- all_line_points()
      if (is.null(d) || !nrow(d)) return(NULL)
      parts <- split(d$`Gene effect`, d$Group)
      parts[["All lines"]] <- d$`Gene effect`
      do.call(rbind, lapply(names(parts), function(g) {
        v <- parts[[g]]
        data.frame(
          Group = g, n = length(v),
          Min = round(min(v), 4), Q1 = round(stats::quantile(v, 0.25, names = FALSE), 4),
          Median = round(stats::median(v), 4),
          Q3 = round(stats::quantile(v, 0.75, names = FALSE), 4),
          Max = round(max(v), 4), Mean = round(mean(v), 4),
          SD = round(stats::sd(v), 4),
          check.names = FALSE, stringsAsFactors = FALSE
        )
      }))
    })

  # The table listed a stripped name and a number, which is not enough to say
  # what a line is. DepMap's model metadata joins cleanly on the stripped name
  # (1,166 of 1,166 lines match), so the line can be described properly.
    effect_points_full <- reactive({
      p <- effect_points_df()
      if (is.null(p) || !nrow(p)) return(NULL)
      meta <- tryCatch(shared_data$cancer_model_data(), error = function(e) NULL)
      df <- p[order(p$subtype, p$value), , drop = FALSE]
      out <- data.frame(
        Subtype = df$subtype,
        `Cell line` = df$cell_line,
        `Gene effect` = round(df$value, 4),
        check.names = FALSE, stringsAsFactors = FALSE
      )
      if (!is.null(meta)) {
        m <- match(df$cell_line, meta$StrippedCellLineName)
        pick <- function(col) if (col %in% names(meta)) meta[[col]][m] else NA
        out$`Model ID` <- pick("ModelID")
        out$`Full name` <- pick("CellLineName")
        out$Lineage <- pick("OncotreeLineage")
        out$`Primary disease` <- pick("OncotreePrimaryDisease")
        out$`Oncotree subtype` <- pick("OncotreeSubtype")
        out$Age <- pick("Age")
        out$`Age category` <- pick("AgeCategory")
        out$Sex <- pick("Sex")
        out$`Primary or metastasis` <- pick("PrimaryOrMetastasis")
      }
      out
    })

    output$effect_points_table <- renderDT({
      df <- all_line_points()
      if (is.null(df)) df <- effect_points_full()
      shiny::validate(need(!is.null(df) && nrow(df) > 0,
                           "No per-cell-line values available for this gene."))
      datatable(cpt_filter_levels(df, max_levels = 40), rownames = FALSE,
                filter = "top", options = list(pageLength = 10, scrollX = TRUE))
    })

    output$dl_effect_points <- downloadHandler(
      filename = function() paste0("cell_lines_", gene_sel() %||% "gene", "_",
                                   Sys.Date(), ".xlsx"),
      content = function(file) {
        df <- effect_points_full()
        rows <- input$effect_points_table_rows_all
        if (!is.null(rows) && length(rows)) df <- df[rows, , drop = FALSE]
        writexl::write_xlsx(list(cell_lines = df), path = file)
      }
    )

    output$dl_effect_box <- downloadHandler(
      filename = function() {
        row <- tryCatch(gene_row(), error = function(e) NULL)
        paste0("gene_effect_", if (is.null(row)) "gene" else row$symbol, "_",
               Sys.Date(), ".png")
      },
      content = function(file) {
        ggplot2::ggsave(file, plot = effect_box_gg(), width = 9, height = 5.5,
                        dpi = 300, bg = "white")
      }
    )

    gene_sites <- reactive({
      atlas <- shared_data$cys_editing_atlas()
      key <- gene_sel()
      req(atlas, key, nzchar(key))
      atlas %>%
        dplyr::mutate(gene_key = cpt_gene_match_key(.data$gene_symbol)) %>%
        dplyr::filter(.data$gene_key == key) %>%
        dplyr::mutate(
          is_func = cpt_is_true(.data$functional),
          is_lig  = cpt_is_true(.data$ligandable),
          # Atlas rows are by definition in the atlas, so these are tiers 1-3;
          # tier 4 exists only for engaged residues the atlas never tested.
          Tier = cpt_evidence_tier(TRUE, .data$is_func, .data$is_lig)
        )
    })

    output$engaged_caption <- renderUI({
      row <- tryCatch(gene_row(), error = function(e) NULL)
      if (is.null(row)) return(NULL)
      n_eng <- row$n_engaged_sites %||% NA
      n_eng <- if (is.na(n_eng)) 0L else as.integer(n_eng)
      n_sites <- row$n_sites %||% NA
      n_sites <- if (is.na(n_sites)) 0L else as.integer(n_sites)
      base <- tags$p(
        class = "text-muted", style = "margin-bottom: 6px;",
        sprintf(
          paste("%d cysteine%s in this protein %s engaged at CR >= 4 — the tiers below",
                "describe those residues. The table lists all %d site%s the atlas tested",
                "in this protein, engaged or not."),
          n_eng, if (n_eng == 1L) "" else "s", if (n_eng == 1L) "is" else "are",
          n_sites, if (n_sites == 1L) "" else "s"
        )
      )
      n_sym <- row$n_symbol_matched %||% NA
      n_sym <- if (is.na(n_sym)) 0L else as.integer(n_sym)
      if (n_sym == 0L) return(base)
      tagList(
        base,
        cpt_note(
          sprintf("%d engaged site%s matched the atlas by gene symbol, not UniProt accession.",
                  n_sym, if (n_sym == 1L) "" else "s"),
          "The atlas row for that residue carries no UniProt accession, so the ",
          "only possible match was gene symbol plus position. Position numbering ",
          "is meaningful within a sequence, not across identifiers, so treat the ",
          "tier for those sites as provisional.",
          tone = "caution"
        )
      )
    })

    output$tier_summary <- renderUI({
      row <- tryCatch(gene_row(), error = function(e) NULL)
      if (is.null(row)) return(NULL)
      counts <- c(
        row$engaged_tier1 %||% NA, row$engaged_tier2 %||% NA,
        row$engaged_tier3 %||% NA, row$engaged_tier4 %||% NA
      )
      counts[is.na(counts)] <- 0L
      if (sum(counts) == 0) {
        return(tags$p(
          class = "text-muted",
          "No cysteine in this gene is engaged by a ligand at CR >= 4, so no ",
          "engaged-site tier applies. Any atlas sites for the gene are listed below."
        ))
      }
      tags$div(
        class = "cpt-tier-row",
        lapply(1:4, function(t) {
          tags$div(
            class = paste0("cpt-tier-chip cpt-tier-", t),
            tags$span(class = "cpt-tier-count", as.integer(counts[[t]])),
            tags$span(class = "cpt-tier-label", CPT_TIER_LABELS[[as.character(t)]])
          )
        })
      )
    })

    output$cys_table <- renderDT({
      sites <- gene_sites()
      req(nrow(sites) > 0)
      df <- sites %>%
        dplyr::arrange(.data$Tier, dplyr::desc(.data$ligandability_score)) %>%
        dplyr::transmute(
          Tier = .data$Tier,
          Site = .data$cysteine_position,
          Context = .data$study_context,
          Editor = .data$editor_support,
          `ABE LFC` = round(.data$abe_mean_lfc, 3),
          `CBE LFC` = round(.data$cbe_mean_lfc, 3),
          `Functional (atlas)` = .data$is_func,
          `Atlas ligandability` = .data$is_lig,
          `Atlas ligandability score` = .data$ligandability_score,
          `Accessibility (RSA)` = round(.data$relative_solvent_accessibility, 3),
          pLDDT = .data$alphafold_plddt,
          ClinVar = cpt_is_true(.data$clinvar_pathogenic)
        )
      datatable(df, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE))
    })

    output$probe_table <- renderDT({
      idx <- index()
      key <- gene_sel()
      req(key, nzchar(key))
      # Every covalent probe recorded for the gene. A CR floor here hid the
      # whole table for genes with no engagement at CR >= 4, which is a fact
      # about the gene worth seeing rather than an empty panel.
      df <- idx$probes %>% dplyr::filter(.data$gene_key == key)

      # DT's own empty text rather than validate(): a validation error inside
      # renderDT leaves the previous gene's widget on screen.
      empty_msg <- "No covalent probes recorded for this gene."

      datatable(
        df %>%
          dplyr::arrange(dplyr::desc(.data$CR)) %>%
          dplyr::select(dplyr::any_of(c("probe_name", "CR", "n_targets",
                                        "cysteineid", "ligandable", "SMILES"))),
        rownames = FALSE,
        options = list(pageLength = 10, scrollX = TRUE,
                       language = list(emptyTable = empty_msg))
      )
    })

    # Hand the selected subtype to the Discover tab instead of making the user
    # re-find it there.
    output$jump_ui <- renderUI({
      df <- tryCatch(gene_dependency(), error = function(e) NULL)
      sel <- input$dependency_table_rows_selected
      if (is.null(df) || !nrow(df) || is.null(sel) || !length(sel)) {
        # A muted line rather than a bare help icon, which would sit at the
        # foot of the box with no control attached to it.
        return(tagList(
          hr(),
          tags$p(class = "text-muted", style = "font-size: 12px; margin: 0;",
                 "Select a row in the dependency table to send that subtype to ",
                 "Discover with the analysis already set up.")
        ))
      }
      subtype <- df$Subtype[sel[1]]
      tagList(
        hr(),
        actionButton(ns("jump_discover"), paste0("Open ", subtype, " in Discover"),
                     icon = icon("arrow-right"), class = "btn-primary btn-block btn-sm")
      )
    })

    observeEvent(input$jump_discover, {
      df <- tryCatch(gene_dependency(), error = function(e) NULL)
      sel <- input$dependency_table_rows_selected
      if (is.null(df) || !nrow(df) || is.null(sel) || !length(sel)) return()
      if (is.function(on_jump)) {
        on_jump(list(subtype = df$Subtype[sel[1]], dataset = active_dataset(),
                     gene = gene_sel()))
      }
    }, ignoreInit = TRUE)

    # Report subject: the comparison's subtype A if set, else the gene's
    # strongest dependency in the selected dataset.
    # Choices: every subtype where this gene is a dependency, defaulting to the
    # one Discover has selected when that is among them.
    observeEvent(list(input$gene, active_dataset(),
                      if (is.null(dep_context)) NULL else dep_context$subtype()), {
      df <- tryCatch(gene_dependency(), error = function(e) NULL)
      ch <- if (is.null(df) || !nrow(df)) character(0) else df$Subtype
      # A report should be producible for any cancer, not only those that clear
      # the dependency filter. The index stores distributions only for the
      # gene x subtype pairs that passed, which for most genes is one or none,
      # so the list comes from the dataset's full subtype roster instead. The
      # report template reads the gene effect matrix directly, so it can render
      # any of them.
      all_st <- character(0)
      st_file <- tryCatch(
        shared_data$cancer_subtype_list_files[[
          if (identical(active_dataset(), "CRISPR")) "CRISPR (23Q4)" else "RNAi"]],
        error = function(e) NULL)
      if (!is.null(st_file)) {
        pth <- shared_data$data_path(st_file)
        # read_cancer_subtype_list_file() strips the file's comment header; a
        # plain readLines() put those lines in the dropdown as subtypes.
        if (file.exists(pth)) all_st <- read_cancer_subtype_list_file(pth)
      }
      others <- setdiff(all_st, ch)
      choices <- if (length(others)) {
        list(`Meets the dependency filter` = as.list(ch),
             `Every other subtype in this dataset` = as.list(others))
      } else if (length(ch)) as.list(ch) else character(0)
      from_discover <- if (!is.null(dep_context) && is.function(dep_context$subtype)) {
        dep_context$subtype()
      } else NULL
      sel <- if (!is.null(from_discover) && from_discover %in% c(ch, others)) {
        from_discover
      } else if (length(ch)) ch[[1]] else if (length(others)) others[[1]] else NULL
      updateSelectInput(session, "report_subtype", choices = choices, selected = sel)
    }, ignoreInit = FALSE)

    report_subtype <- reactive({
      if (!is.null(input$report_subtype) && nzchar(input$report_subtype)) {
        return(input$report_subtype)
      }
      df <- tryCatch(gene_dependency(), error = function(e) NULL)
      if (is.null(df) || !nrow(df)) return(NULL)
      as.character(df$Subtype[[1]])
    })

    output$dl_target_report <- downloadHandler(
      filename = function() {
        row <- tryCatch(gene_row(), error = function(e) NULL)
        sym <- if (is.null(row)) "gene" else row$symbol
        paste0("CanProTarget_", sym, "_", Sys.Date(), ".html")
      },
      content = function(file) {
        row <- tryCatch(gene_row(), error = function(e) NULL)
        st <- report_subtype()
        if (is.null(row) || is.null(st)) {
          writeLines(paste0(
            "<!DOCTYPE html><html><body><p>No dependency subtype for this gene in ",
            active_dataset(), " data, so there is nothing to report on.</p></body></html>"), file)
          return(invisible(NULL))
        }
        withProgress(message = "Generating report...", value = 0.3, {
          out <- tryCatch(
            cpt_report_gene_dependency(
              gene = row$symbol, subtype = st,
              dataset = if (identical(active_dataset(), "CRISPR")) "CRISPR" else "RNAi",
              project_root = getwd()
            ),
            error = function(e) {
              showNotification(paste("Report error:", conditionMessage(e)), type = "error")
              NULL
            }
          )
          incProgress(0.6)
          if (is.null(out) || !file.exists(out)) {
            writeLines("<!DOCTYPE html><html><body><p>Report could not be rendered.</p></body></html>", file)
          } else {
            file.copy(out, file, overwrite = TRUE)
          }
        })
      }
    )

    # Everything the plot is drawn from: the box summaries, the individual cell
    # lines, and the two-subtype comparison if one has been run.
    output$dl_effect_data <- downloadHandler(
      filename = function() {
        row <- tryCatch(gene_row(), error = function(e) NULL)
        paste0("gene_effect_data_", if (is.null(row)) "gene" else row$symbol, "_",
               Sys.Date(), ".xlsx")
      },
      content = function(file) {
        sheets <- list(
          box_summaries = tryCatch(
            effect_box_df()[, setdiff(names(effect_box_df()), "role")],
            error = function(e) NULL),
          cell_lines = tryCatch(
            effect_points_df()[, c("subtype", "cell_line", "value")],
            error = function(e) NULL),
          comparison = tryCatch({
            r <- cmp_out()
            data.frame(
              subtype = c(r$a$subtype, r$b$subtype),
              n = c(r$a$n, r$b$n), mean = c(r$a$mean, r$b$mean),
              sd = c(r$a$sd, r$b$sd),
              limma_effect_size = c(r$a$effect, r$b$effect),
              limma_p_vs_all_others = c(r$a$p, r$b$p),
              stringsAsFactors = FALSE
            )
          }, error = function(e) NULL)
        )
        sheets <- sheets[!vapply(sheets, is.null, logical(1))]
        if (!length(sheets)) sheets <- list(empty = data.frame())
        writexl::write_xlsx(sheets, file)
      }
    )

    output$dep_context <- renderText({
      if (!is.null(dep_context) && is.function(dep_context$context_label)) {
        dep_context$context_label()
      } else {
        active_dataset()
      }
    })


    dl_table <- function(stub, data_fn) {
      downloadHandler(
        filename = function() {
          row <- tryCatch(gene_row(), error = function(e) NULL)
          paste0(stub, "_", if (is.null(row)) "gene" else row$symbol, "_",
                 Sys.Date(), ".xlsx")
        },
        content = function(file) {
          df <- tryCatch(data_fn(), error = function(e) NULL)
          writexl::write_xlsx(if (is.null(df)) data.frame() else df, file)
        }
      )
    }

    output$dl_dep_table <- dl_table("dependency", function() gene_dependency())
    output$dl_cys_table <- dl_table("cysteine_sites", function() {
      gene_sites() %>% dplyr::select(-dplyr::any_of(c("is_func", "is_lig")))
    })
    output$dl_probe_table <- dl_table("probes", function() {
      idx <- index()
      df <- idx$probes %>% dplyr::filter(.data$gene_key == gene_sel())
      df
    })
    output$dl_cmp <- dl_table("subtype_comparison", function() {
      r <- cmp_out()
      data.frame(
        subtype = c(r$a$subtype, r$b$subtype), n = c(r$a$n, r$b$n),
        mean = c(r$a$mean, r$b$mean), sd = c(r$a$sd, r$b$sd),
        limma_effect_size = c(r$a$effect, r$b$effect),
        limma_p_vs_all_others = c(r$a$p, r$b$p), stringsAsFactors = FALSE
      )
    })

    output$dl_gene <- downloadHandler(
      filename = function() {
        row <- tryCatch(gene_row(), error = function(e) NULL)
        sym <- if (is.null(row)) "gene" else row$symbol
        paste0("canprotarget_", sym, "_", Sys.Date(), ".xlsx")
      },
      content = function(file) {
        sheets <- list(
          dependency = tryCatch(gene_dependency(), error = function(e) NULL),
          cysteine_sites = tryCatch(
            gene_sites() %>% dplyr::select(-dplyr::any_of(c("is_func", "is_lig"))),
            error = function(e) NULL
          ),
          probes = tryCatch({
            idx <- index()
            idx$probes %>% dplyr::filter(.data$gene_key == gene_sel())
          }, error = function(e) NULL)
        )
        sheets <- sheets[!vapply(sheets, is.null, logical(1))]
        if (!length(sheets)) sheets <- list(empty = data.frame())
        writexl::write_xlsx(sheets, file)
      }
    )

    # Expose the current gene so other tabs can follow along.
    list(current_gene = gene_sel)
  })
}

# ============================================================
# Shared gene card
# ============================================================
#' Compact, static summary of one gene, built from the precomputed index.
#'
#' Used by the Discover drawer so gene information has one source rather than
#' a second implementation living next to the table. The Gene tab remains the
#' full interactive view; this is the at-a-glance version.
#'
#' @param idx data/gene_index.rds contents
#' @param gene_key Upper-case bare symbol
#' @param dataset "CRISPR" or "RNAi"
#' @param n_show Rows per section
#' @return shiny.tag
cpt_gene_card <- function(idx, gene_key, dataset = "RNAi", n_show = 5) {
  if (is.null(idx) || is.null(gene_key) || !nzchar(gene_key)) return(NULL)
  row <- idx$genes[idx$genes$gene_key == gene_key, , drop = FALSE]
  if (!nrow(row)) return(shiny::tags$p("No index entry for this gene."))
  row <- row[1, , drop = FALSE]

  num <- function(x, digits = 3) {
    if (is.null(x) || length(x) == 0 || is.na(x)) "—" else format(round(x, digits))
  }
  count <- function(x) if (is.null(x) || length(x) == 0 || is.na(x)) 0L else as.integer(x)

  dep <- idx$dependency
  dep <- dep[dep$gene_key == gene_key & dep$dataset == dataset, , drop = FALSE]
  dep <- dep[order(dep$effect_size), , drop = FALSE]
  dep <- utils::head(dep, n_show)

  pr <- idx$probes
  pr <- pr[pr$gene_key == gene_key & !is.na(pr$CR) & pr$CR >= 4, , drop = FALSE]
  pr <- pr[order(-pr$CR), , drop = FALSE]
  pr <- utils::head(pr, n_show)

  stat <- function(label, value) {
    shiny::tags$div(
      class = "cpt-card-stat",
      shiny::tags$span(class = "cpt-card-stat-n", value),
      shiny::tags$span(class = "cpt-card-stat-l", label)
    )
  }

  mini_rows <- function(df, left, right, empty) {
    if (is.null(df) || !nrow(df)) {
      return(shiny::tags$p(class = "text-muted", style = "font-size: 12px;", empty))
    }
    shiny::tags$div(
      class = "cpt-card-rows",
      lapply(seq_len(nrow(df)), function(i) {
        shiny::tags$div(
          class = "cpt-card-row",
          shiny::tags$span(class = "cpt-card-row-l", left(df[i, ])),
          shiny::tags$span(class = "cpt-card-row-r", right(df[i, ]))
        )
      })
    )
  }

  shiny::tagList(
    shiny::tags$div(
      class = "cpt-gene-header",
      shiny::tags$h2(row$symbol),
      cpt_tier_badge(row$best_engaged_tier)
    ),
    shiny::tags$div(
      class = "cpt-card-stats",
      stat("subtypes", count(row$dep_n_subtypes)),
      stat("engaged sites", count(row$n_engaged_sites)),
      stat("probes", count(row$n_probes)),
      stat("ClinVar", count(row$n_clinvar))
    ),
    shiny::tags$h5(class = "cpt-card-h", paste0("Strongest dependencies (", dataset, ")")),
    mini_rows(dep, function(r) r$subtype, function(r) num(r$effect_size),
              "Not a selective dependency in this dataset."),
    shiny::tags$h5(class = "cpt-card-h", "Top engaged probes (CR ≥ 4)"),
    mini_rows(pr, function(r) r$probe_name, function(r) num(r$CR, 2),
              "No cysteine engaged at CR ≥ 4.")
  )
}
