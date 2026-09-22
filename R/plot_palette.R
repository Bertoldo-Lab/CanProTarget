# ============================================================
# Script:   plot_palette.R
# Purpose:  One visual language for every chart in the app.
#
#           One palette, marker size and grid for every chart, built from
#           the site's own blue so a figure reads as part of the page.
#           A given meaning ("this row was excluded") has one colour
#           throughout.
# ============================================================
# Sections:
#   cpt_plotly() — Shared plotly layout + legend-hover highlighting
#   cpt_axis()
# ============================================================

# Site blues, dark to pale. Sequence matters: use them in order for ordered
# categories (evidence tiers), so darker always means stronger evidence.
CPT_PAL <- list(
  ink       = "#1a3a5c",  # darkest — the series that matters most
  deep      = "#2C6A94",
  primary   = "#478EB8",  # the header blue
  light     = "#7bafd4",
  pale      = "#b8d4e8",
  muted     = "#c8d4dc",  # excluded / de-emphasised
  grid      = "#e6ecf0",
  alert     = "#cc3340",  # cutoffs, out-of-range values
  caution   = "#e0a800",  # needs attention but not an error
  text      = "#1a3a5c",
  text_soft = "#7b8a94"
)

# Evidence tiers 1-4. Hue carries the meaning rather than lightness alone:
# blue means the residue has functional evidence (tiers 1-2), grey means it
# does not (tier 3 tested-and-negative, tier 4 untested). Navy and mid blue
# are far enough apart to separate tiers 1 and 2 at a glance.
CPT_TIER_PAL <- c(
  "1" = CPT_PAL$ink,      # navy
  "2" = CPT_PAL$primary,  # mid blue - clearly separable from navy
  "3" = "#7d8b97",        # mid grey, solid
  "4" = "#aab6c0"         # light grey, still readable on white
)

# Tier 4 is about 82% of all points, so it is drawn lighter than the rest and
# behind them. It was #dde3e8 at half opacity, which on a white panel left the
# points barely visible at all -- in a scatter where tier 4 is most of the
# data. #aab6c0 is the lightest grey that still keeps a perceptual distance of
# about 11 from every other colour in the palette, the nearest being the muted
# grey used for excluded rows, and it stays clearly lighter than tier 3's
# solid #7d8b97. The opacity comes up to match.
CPT_TIER_ALPHA <- c("1" = 0.9, "2" = 0.9, "3" = 0.85, "4" = 0.7)

# One marker size everywhere, so a dot in one chart means the same as a dot
# in another, and large enough to hit with a cursor.
CPT_MARKER_SIZE <- 7
CPT_MARKER_SIZE_LG <- 11   # emphasis markers (out-of-range flags)

#' Shared plotly layout + legend-hover highlighting.
#'
#' Hovering a legend entry dims every other trace, which is the only way to
#' read a dense overplotted scatter without clicking entries off one at a time.
#'
#' @param p A plotly object
#' @param ... Passed to plotly::layout()
cpt_plotly <- function(p, ...) {
  p <- plotly::layout(
    p,
    font = list(family = "'Source Sans Pro', 'Helvetica Neue', Arial, sans-serif",
                size = 13, color = CPT_PAL$text),
    paper_bgcolor = "rgba(0,0,0,0)",
    plot_bgcolor = "rgba(0,0,0,0)",
    hoverlabel = list(bgcolor = "#ffffff", bordercolor = CPT_PAL$primary,
                      font = list(size = 12.5, color = CPT_PAL$text)),
    hovermode = "closest",
    ...
  )
  p <- plotly::config(p, displaylogo = FALSE, responsive = TRUE,
                      modeBarButtonsToRemove = c("select2d", "lasso2d",
                                                 "hoverCompareCartesian"))
  htmlwidgets::onRender(p, "
    function(el) {
      var full = el.data.map(function(t) {
        return (t.marker && t.marker.opacity !== undefined) ? t.marker.opacity : 1;
      });
      el.on('plotly_legendhover', function(d) {
        var dim = el.data.map(function(_, i) {
          return i === d.curveNumber ? full[i] : 0.08;
        });
        Plotly.restyle(el, {'marker.opacity': dim});
        return false;
      });
      el.on('plotly_legenddoubleclick', function() { return false; });
      el.on('plotly_unhover', function() {
        Plotly.restyle(el, {'marker.opacity': full});
      });
      el.addEventListener('mouseleave', function() {
        Plotly.restyle(el, {'marker.opacity': full});
      });
      // A chart drawn inside a hidden tab or conditionalPanel measures a
      // container that is not laid out yet, so it keeps that stale width when
      // the panel is shown and overhangs its box. Re-measure on any size
      // change of the parent, and once when it first becomes visible.
      var fit = function() {
        if (el.offsetParent !== null) { Plotly.Plots.resize(el); }
      };
      if (window.ResizeObserver && el.parentNode) {
        new ResizeObserver(fit).observe(el.parentNode);
      }
      if (window.IntersectionObserver) {
        new IntersectionObserver(function(entries) {
          entries.forEach(function(e) { if (e.isIntersecting) { fit(); } });
        }).observe(el);
      }
      window.addEventListener('resize', fit);
    }
  ")
}

#' Standard axis spec, so every chart's gridlines and zero line match.
cpt_axis <- function(title, ...) {
  list(title = title, gridcolor = CPT_PAL$grid, zeroline = FALSE,
       linecolor = CPT_PAL$grid, tickcolor = CPT_PAL$grid, ...)
}
