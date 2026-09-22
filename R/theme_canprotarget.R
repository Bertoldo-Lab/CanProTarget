# ============================================================
# Script:   theme_canprotarget.R
# Purpose:  Publication-quality ggplot2 theme and color palette
#           for all CanProTarget figures (in-app, reports, exports).
# ============================================================
# Sections:
#   Color Palette (derived from logo)
#   Publication Theme
#     theme_canprotarget()
#   Scale helpers
#     scale_fill_dependency()
# ============================================================

# ---- Color Palette (derived from logo) -------------------------

cpt_colors <- list(
  dark_blue   = "#1a3a5c",
  medium_blue = "#4a7fb5",
  light_blue  = "#7bafd4",
  pale_blue   = "#b8d4e8",
  accent      = "#478EB8",
  danger      = "#cc3340",
  warning     = "#e0a800",
  grey_dark   = "#434343",
  grey_medium = "#666666",
  grey_light  = "#999999",
  white       = "#FFFFFF"
)

# Named vector for quick programmatic access
cpt_palette <- unlist(cpt_colors)

# Discrete color scales (ordered for common plot use cases)
cpt_discrete_colors <- c(
  "#478EB8", "#cc3340", "#4a7fb5", "#e0a800",
  "#1a3a5c", "#7bafd4", "#434343", "#b8d4e8"
)

# Sequential blue palette for heatmaps (low to high)
cpt_seq_blue <- c("#f7fbff", "#deebf7", "#b8d4e8", "#7bafd4",
                  "#4a7fb5", "#2171b5", "#1a3a5c", "#08306b")

# Diverging palette: blue (not dependent) -> white -> red (strongly dependent)
cpt_diverging <- c("#1a3a5c", "#4a7fb5", "#7bafd4", "#b8d4e8",
                   "#FFFFFF",
                   "#fcbba1", "#fc9272", "#cc3340", "#8b1f27")

# ---- Publication Theme -----------------------------------------

#' CanProTarget ggplot2 theme for publication-quality figures.
#'
#' Clean, minimal theme with branded typography. Suitable for
#' journal figures, reports, and in-app display.
#'
#' @param base_size Base font size in points (default 11)
#' @param base_family Base font family (default "")
#' @return A ggplot2 theme object
theme_canprotarget <- function(base_size = 11, base_family = "") {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      # Title and subtitle
      plot.title = ggplot2::element_text(
        face = "bold", color = cpt_colors$dark_blue,
        size = ggplot2::rel(1.2), margin = ggplot2::margin(b = 6)
      ),
      plot.subtitle = ggplot2::element_text(
        color = cpt_colors$medium_blue,
        size = ggplot2::rel(0.95), margin = ggplot2::margin(b = 10)
      ),
      plot.caption = ggplot2::element_text(
        color = cpt_colors$grey_medium,
        size = ggplot2::rel(0.8), hjust = 1
      ),

      # Axes
      axis.title = ggplot2::element_text(
        color = cpt_colors$grey_dark, size = ggplot2::rel(0.95)
      ),
      axis.text = ggplot2::element_text(
        color = cpt_colors$grey_dark, size = ggplot2::rel(0.85)
      ),
      axis.ticks = ggplot2::element_line(color = "#cccccc", linewidth = 0.3),

      # Grid
      panel.grid.major = ggplot2::element_line(color = "#eeeeee", linewidth = 0.4),
      panel.grid.minor = ggplot2::element_blank(),

      # Facets
      strip.text = ggplot2::element_text(
        face = "bold", color = cpt_colors$dark_blue,
        size = ggplot2::rel(0.9)
      ),
      strip.background = ggplot2::element_rect(fill = "#f0f4f8", color = NA),

      # Legend
      legend.position = "bottom",
      legend.title = ggplot2::element_text(
        face = "bold", size = ggplot2::rel(0.85)
      ),
      legend.text = ggplot2::element_text(size = ggplot2::rel(0.8)),
      legend.key.size = ggplot2::unit(0.8, "lines"),

      # Plot margins
      plot.margin = ggplot2::margin(10, 10, 10, 10)
    )
}

# ---- Scale helpers ---------------------------------------------

#' Diverging fill for dependency heatmaps (blue = not dep, red = strongly dep)
#' @param ... Additional arguments passed to scale_fill_gradient2
scale_fill_dependency <- function(midpoint = 0, ...) {
  ggplot2::scale_fill_gradient2(
    low = cpt_colors$medium_blue,
    mid = cpt_colors$white,
    high = cpt_colors$danger,
    midpoint = midpoint,
    ...
  )
}
