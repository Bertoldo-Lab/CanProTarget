# ============================================================
# Script:   cys_editing_module.R
# Purpose:  Explore the Cravatt Lab base-editing atlas and
#           intersect it with the current CanProTarget analysis.
# ============================================================
# Sections:
#   cys_editing_table_ui()
#   cys_editing_chart_ui()
#   cys_editing_server()
# ============================================================

#' @param wrap TRUE returns a standalone tabItem; FALSE returns the body so it
#'   can be embedded (it now lives inside Discover).
# The atlas is the residue evidence layer's contribution to Discover's shared
# table and chart spaces, so it is exposed as pieces the parent can place rather
# than as one page of its own.

#' The cysteine-level evidence table, its download and its provenance note.
cys_editing_table_ui <- function(id) {
  ns <- NS(id)
  tagList(
    div(class = "cpt-scroll-x",
        shinycssloaders::withSpinner(DTOutput(ns("atlas_table")), type = 6,
                                     color = "#478EB8", hide.ui = FALSE)),
    br(),
    downloadButton(ns("download_atlas"), "Download filtered CSV",
                   class = "btn-sm btn-default"),
    tags$p(
      class = "text-muted", style = "margin-top: 12px; font-size: 12px;",
      "A site is called functional on the publication thresholds: mean ",
      "LFC \u2264 -0.6, empirical p < 0.05, and within-gene FDR < 0.1 after ",
      "cell-context filtering. Ligandability and structural annotations are ",
      "joined per residue, never per gene \u2014 functional evidence at one ",
      "cysteine is not evidence for another in the same protein. Source: ",
      tags$a("cravattlab/Cys_editing",
             href = "https://github.com/cravattlab/Cys_editing",
             target = "_blank", rel = "noopener noreferrer"),
      " \u00b7 ",
      tags$a("Li et al., Nature Chemical Biology (2023)",
             href = "https://doi.org/10.1038/s41589-023-01428-w",
             target = "_blank", rel = "noopener noreferrer")
    )
  )
}

#' One atlas chart, named. The parent's chart switcher chooses between them.
cys_editing_chart_ui <- function(id, which = c("abe_cbe", "overview")) {
  ns <- NS(id)
  which <- match.arg(which)
  out <- if (identical(which, "abe_cbe")) "lfc_plot" else "gene_plot"
  shinycssloaders::withSpinner(
    plotlyOutput(ns(out), height = "440px"),
    type = 6, color = "#478EB8", hide.ui = FALSE)
}


cys_editing_server <- function(id, shared_data, dependency_genes = NULL,
                               filters = NULL) {
  moduleServer(id, function(input, output, session) {
    # This view is filtered by the Cysteine function box and the action bar, not
    # by controls of its own, so one control drives one thing. `filters` is the
    # list of reactives supplying them; each has a fallback for a standalone
    # render where no filters are passed.
    fval <- function(name, default = NULL) {
      f <- if (is.null(filters)) NULL else filters[[name]]
      if (is.null(f)) return(default)
      v <- tryCatch(f(), error = function(e) default)
      if (is.null(v)) default else v
    }
    atlas <- reactive({
      data <- shared_data$cys_editing_atlas()
      if (is.null(data)) {
        cpt_notify_error("CPT-4001", "Cysteine atlas not available",
                         detail = "cys_editing_atlas.rds missing or failed to load")
        return(NULL)
      }
      data
    })

    observeEvent(atlas(), {
      data <- atlas()
      if (is.null(data) || nrow(data) == 0) return()
      genes <- sort(unique(as.character(data$gene_symbol)))
      updateSelectizeInput(session, "genes", choices = genes, server = TRUE)
    }, ignoreInit = FALSE)

    current_dependency_genes <- reactive({
      if (is.null(dependency_genes)) {
        return(character(0))
      }
      # Dependency gene names arrive with DepMap's Entrez suffix ("KRAS (3845)")
      # while the atlas uses bare symbols, so toupper() alone never matches.
      # cpt_gene_match_key() strips the suffix and upper-cases. See functions.R.
      tryCatch(
        unique(cpt_gene_match_key(dependency_genes())),
        error = function(e) character(0)
      )
    })

    filtered_atlas <- reactive({
      data <- atlas()
      if (is.null(data) || nrow(data) == 0) {
        return(data.frame(stringsAsFactors = FALSE))
      }

      genes <- fval("genes")
      if (length(genes)) {
        data <- data[cpt_gene_match_key(data$gene_symbol) %in% genes, , drop = FALSE]
      }
      # Empty means any context here, matching the Cysteine function box's own
      # wording. The old box defaulted to all three selected and read an empty
      # selection as "show nothing", which is the opposite.
      ctx <- fval("contexts")
      if (length(ctx)) {
        data <- data[data$study_context %in% ctx, , drop = FALSE]
      }

      # Functional vs not is carried by the evidence tier below: tiers 1 and 2
      # are functional, tier 3 is tested and not functional.

      editor <- fval("editor", "Any editor")
      if (identical(editor, "ABE")) {
        data <- data[data$abe_functional, , drop = FALSE]
      } else if (identical(editor, "CBE")) {
        data <- data[data$cbe_functional, , drop = FALSE]
      } else if (identical(editor, "ABE and CBE")) {
        data <- data[data$abe_functional & data$cbe_functional, , drop = FALSE]
      }

      if (isTRUE(fval("atlas_lig", FALSE))) {
        data <- data[data$ligandable, , drop = FALSE]
      }

      if (isTRUE(fval("clinvar", FALSE))) {
        if ("clinvar_pathogenic" %in% colnames(data)) {
          data <- data[data$clinvar_pathogenic, , drop = FALSE]
        }
      }

      if (isTRUE(fval("dep_only", FALSE))) {
        dep <- current_dependency_genes()
        data <- data[toupper(data$gene_symbol) %in% dep, , drop = FALSE]
      }

      # Evidence tier, from the same helper Discover and the Gene tab use.
      data$evidence_tier <- cpt_evidence_tier(
        TRUE, cpt_is_true(data$functional), cpt_is_true(data$ligandable)
      )
      keep_tiers <- fval("tiers", c("1", "2", "3"))
      data <- data[as.character(data$evidence_tier) %in% keep_tiers, , drop = FALSE]

      # Engagement is a property of the probe table, not the atlas, so it comes
      # from the precomputed engaged-site list (CR >= 4) by atlas site id.
      if (isTRUE(fval("engaged", FALSE))) {
        idx <- tryCatch(shared_data$gene_index(), error = function(e) NULL)
        engaged_ids <- if (!is.null(idx) && !is.null(idx$engaged_sites)) {
          unique(idx$engaged_sites$cys_site_id)
        } else {
          character(0)
        }
        data <- data[data$site_id %in% engaged_ids, , drop = FALSE]
      }

      data
    })


    output$lfc_plot <- renderPlotly({
      data <- filtered_atlas()
      shiny::validate(need(
        is.data.frame(data) && nrow(data) > 0,
        "No cysteines match the current filters."
      ))

      data$evidence_class <- ifelse(
        data$functional_ligandable,
        "Functional + ligandable",
        ifelse(data$functional, "Functional", "Not significant")
      )

      # Build hover text safely (all character, no NAs in display)
      conservation_text <- if ("conservation_score" %in% colnames(data)) {
        ifelse(is.na(data$conservation_score), "",
               paste0("<br>Conservation: ", round(data$conservation_score * 100, 0), "%"))
      } else {
        rep("", nrow(data))
      }

      clinvar_text <- if ("clinvar_pathogenic" %in% colnames(data)) {
        annotation <- as.character(data$clinvar_annotation)
        annotation[is.na(annotation)] <- ""
        ifelse(data$clinvar_pathogenic,
               paste0("<br>ClinVar: ", annotation),
               "")
      } else {
        rep("", nrow(data))
      }

      data$hover_text <- paste0(
        "<b>", as.character(data$site_id), "</b>",
        "<br>UniProt: ", ifelse(is.na(data$uniprot_accession), "not mapped", as.character(data$uniprot_accession)),
        "<br>Context: ", as.character(data$study_context),
        "<br>Editor support: ", as.character(data$editor_support),
        "<br>ABE mean LFC: ", round(data$abe_mean_lfc, 3),
        "<br>CBE mean LFC: ", round(data$cbe_mean_lfc, 3),
        "<br>Ligandability: ", ifelse(data$ligandable, "yes", "no"),
        conservation_text,
        clinvar_text
      )

      plot <- plotly::plot_ly(
        data = data,
        x = ~abe_mean_lfc,
        y = ~cbe_mean_lfc,
        color = ~evidence_class,
        colors = c(
          "Functional + ligandable" = CPT_PAL$ink,
          "Functional" = CPT_PAL$primary,
          "Not significant" = "#c3ccd3"
        ),
        text = ~hover_text,
        hoverinfo = "text",
        type = "scatter",
        mode = "markers",
        marker = list(size = CPT_MARKER_SIZE, opacity = 0.75,
                      line = list(width = 0))
      )

      cpt_plotly(
        plot,
        xaxis = cpt_axis("ABE mean log2 fold-change", zeroline = TRUE,
                         zerolinecolor = "#d5dde2"),
        yaxis = cpt_axis("CBE mean log2 fold-change", zeroline = TRUE,
                         zerolinecolor = "#d5dde2"),
        legend = list(orientation = "h", x = 0, y = -0.18),
        margin = list(b = 95)
      )
    })

    # Two further readings of the same filtered set: which proteins carry the
    # most evidence, and how hard the dropout is in each tier. Both use the
    # shared tier palette so a colour means the same thing in every chart.
    # What the current filters actually left: how many sites, at which tiers,
    # in which study contexts. A per-protein ranking answered a narrower
    # question than the one a reader has when they open this view.
    output$gene_plot <- renderPlotly({
      data <- filtered_atlas()
      shiny::validate(need(is.data.frame(data) && nrow(data) > 0,
                           "No cysteines match the current filters."))
      data$tier <- as.character(data$evidence_tier)
      ctx <- as.character(data$study_context)
      ctx[is.na(ctx) | !nzchar(ctx)] <- "unspecified"
      counts <- as.data.frame(table(context = ctx, tier = data$tier),
                              stringsAsFactors = FALSE)
      counts <- counts[counts$Freq > 0, , drop = FALSE]
      order_ctx <- names(sort(tapply(counts$Freq, counts$context, sum)))
      p <- plot_ly()
      for (tr in sort(unique(counts$tier))) {
        sub <- counts[counts$tier == tr, , drop = FALSE]
        p <- add_trace(
          p, data = sub, x = ~Freq, y = ~context, type = "bar", orientation = "h",
          name = paste("Tier", tr),
          marker = list(color = unname(CPT_TIER_PAL[tr]),
                        line = list(color = "#ffffff", width = 0.5)),
          hovertemplate = paste0("%{y}<br>Tier ", tr, ": %{x} sites<extra></extra>")
        )
      }
      cpt_plotly(
        layout(p, barmode = "stack",
               xaxis = cpt_axis("Cysteine sites"),
               yaxis = cpt_axis("", automargin = TRUE,
                                categoryorder = "array",
                                categoryarray = order_ctx),
               legend = list(orientation = "h", y = -0.15))
      )
    })

    # Column order, and the heading each one is shown under. The rest of the
    # app labels its tables in plain English; the atlas was the last place
    # still showing the source frame's snake_case names.
    display_columns <- c(
      evidence_tier                  = "Tier",
      site_id                        = "Site",
      gene_symbol                    = "Gene",
      protein_name                   = "Protein",
      uniprot_accession              = "UniProt",
      residue_mapping_status         = "Residue mapping",
      study_context                  = "Study context",
      editor_support                 = "Editor support",
      abe_mean_lfc                   = "ABE mean LFC",
      cbe_mean_lfc                   = "CBE mean LFC",
      abe_guide_count                = "ABE guides",
      cbe_guide_count                = "CBE guides",
      abe_neg_log10_fdr              = "ABE -log10 FDR",
      cbe_neg_log10_fdr              = "CBE -log10 FDR",
      functional                     = "Functional",
      ligandability_score            = "Ligandability score",
      ligandable                     = "Ligandable",
      functional_ligandable          = "Functional + ligandable",
      proteomic_accessibility        = "Proteomic accessibility",
      relative_solvent_accessibility = "Solvent accessibility",
      alphafold_plddt                = "AlphaFold pLDDT",
      conservation_score             = "Conservation",
      ortholog_cys_count             = "Orthologs with cysteine",
      ortholog_total                 = "Orthologs compared",
      clinvar_pathogenic             = "ClinVar pathogenic",
      clinvar_annotation             = "ClinVar annotation",
      clinvar_phenotype              = "ClinVar phenotype"
    )

    output$atlas_table <- renderDT({
      data <- filtered_atlas()
      available_cols <- intersect(names(display_columns), colnames(data))
      data <- data[, available_cols, drop = FALSE]
      numeric_columns <- vapply(data, is.numeric, logical(1))
      data[numeric_columns] <- lapply(data[numeric_columns], function(x) round(x, 3))
      # No Scroller extension here. With deferRender the body is virtualised
      # and the per-column filter boxes render but never filter, which is why
      # search worked in every other table and not this one. Plain paging, the
      # same as the rest of the app.
      datatable(
        data,
        rownames = FALSE,
        colnames = unname(display_columns[available_cols]),
        filter = "top",
        options = list(
          pageLength = 15,
          scrollX = TRUE,
          language = list(emptyTable =
            "No cysteine site matches the current filters.")
        )
      )
    })

    output$download_atlas <- downloadHandler(
      filename = function() paste0("CanProTarget_Cys_editing_", Sys.Date(), ".csv"),
      content = function(file) {
        data <- filtered_atlas()
        if (is.null(data) || nrow(data) == 0) {
          utils::write.csv(data.frame(message = "No data matches current filters"),
                          file, row.names = FALSE)
        } else {
          utils::write.csv(data, file, row.names = FALSE, na = "")
        }
      }
    )

  })
}
