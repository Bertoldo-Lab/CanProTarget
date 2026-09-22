#!/usr/bin/env Rscript --vanilla
# ============================================================
# Script:   capture_screenshots.R
# Purpose:  Re-shoot the screenshots in docs/screenshots from a running
#           CanProTarget, so the images in the README and the docs match
#           what the app currently looks like.
#
# Usage:    Rscript --vanilla docs/scripts/capture_screenshots.R [url] [outdir]
#           Defaults to the live app and docs/screenshots.
#
# Notes:    Each shot drives the real interface rather than a fixture, so the
#           waits are generous: a dependency run against the live app takes
#           the better part of a minute. Shots are full-page, which means the
#           viewport is resized to the document height before capture.
# ============================================================

suppressPackageStartupMessages(library(chromote))

args   <- commandArgs(trailingOnly = TRUE)
url    <- if (length(args) >= 1) args[[1]] else "https://bertoldolab.shinyapps.io/CanProTarget/"
outdir <- if (length(args) >= 2) args[[2]] else "docs/screenshots"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

WIDTH <- 1600L
# Shots are taken at twice the CSS width so the text stays readable when the
# images are viewed at full size. chromote ignores deviceScaleFactor here; the
# scale argument to $screenshot() is the one that works.
SCALE <- 2
b <- ChromoteSession$new()
on.exit(try(b$close(), silent = TRUE), add = TRUE)

set_size <- function(h) {
  b$Emulation$setDeviceMetricsOverride(width = WIDTH, height = as.integer(h),
                                       deviceScaleFactor = 1, mobile = FALSE)
}
js <- function(code) {
  out <- b$Runtime$evaluate(code, awaitPromise = FALSE, returnByValue = TRUE)
  invisible(out$result$value)
}
jsv <- function(code) {
  b$Runtime$evaluate(code, awaitPromise = FALSE, returnByValue = TRUE)$result$value
}
# Shiny's own busy/idle events, not the DOM.
#
# Counting .recalculating or .load-container does not work in this app: most
# outputs live inside conditionalPanels that are closed, so they never render
# and never lose the class. On a quiet page 65 elements still claim to be
# recalculating. The shiny:busy / shiny:idle events are the supported signal
# and track the actual request queue.
say <- function(...) { cat(format(Sys.time(), "%H:%M:%S"), "", ..., "\n"); flush.console() }

install_idle_hook <- function() {
  js("(function(){
        if (window.__cptHook) return 'already';
        window.__cptHook = true; window.__cptIdle = true;
        $(document).on('shiny:busy', function(){ window.__cptIdle = false; });
        $(document).on('shiny:idle', function(){ window.__cptIdle = true; });
        return 'ok';
      })()")
}
wait_idle <- function(timeout = 180, quiet_for = 3) {
  t0 <- Sys.time(); quiet <- 0
  repeat {
    idle <- tryCatch(isTRUE(jsv("window.__cptIdle === true")), error = function(e) FALSE)
    quiet <- if (idle) quiet + 1 else 0
    if (quiet >= quiet_for) return(TRUE)
    if (as.numeric(difftime(Sys.time(), t0, units = "secs")) > timeout) return(FALSE)
    Sys.sleep(1)
  }
}
shoot <- function(name) {
  say("  waiting for idle before", name)
  wait_idle()
  h <- tryCatch(jsv("Math.max(document.body.scrollHeight, document.documentElement.scrollHeight)"),
                error = function(e) 1200)
  h <- max(800L, min(as.integer(h) + 40L, 6000L))
  set_size(h)
  Sys.sleep(2)
  path <- file.path(outdir, name)
  b$screenshot(filename = path, scale = SCALE)
  cat(sprintf("  %-34s %s  (%d x %d)\n", name,
              format(file.size(path) / 1024, digits = 4), WIDTH * SCALE, h * SCALE))
  set_size(1100)
}
nav <- function(tab) {
  js(sprintf("(function(){var a=document.querySelector('a[href=\"#shiny-tab-%s\"]'); if(a) a.click(); return 'ok';})()", tab))
  Sys.sleep(4)
}
# Layer checkboxes are the only unnamed-id inputs on Discover, so they are
# addressed by value rather than by id.
set_layers <- function(vals) {
  say("  layers ->", paste(vals, collapse = "+"))
  js(sprintf(
    "(function(){var want=%s;
       [...document.querySelectorAll('#shiny-tab-gene_effect input[type=checkbox]')]
         .filter(c=>['dependency','ligandability','functionality'].includes(c.value))
         .forEach(function(c){ if(c.checked !== want.includes(c.value)) c.click(); });
       return 'ok';})()",
    paste0("[", paste(sprintf("'%s'", vals), collapse = ","), "]")))
  Sys.sleep(3)
}
run_analysis <- function() {
  say("  run analysis")
  js("(function(){var b=document.getElementById('deps-run_analysis'); if(b && !b.disabled) b.click(); return 'ok';})()")
  Sys.sleep(6)
  ok <- wait_idle()
  say("  analysis settled:", ok)
  ok
}

# Pick a named view out of the table or chart switcher. The switcher is a row
# of plain buttons driving one Shiny input, so a click is the whole interaction.
pick <- function(kind, label) {
  sw <- if (identical(kind, "table")) "deps-tp_table_switch" else "deps-tp_chart_switch"
  jsv(sprintf("(function(){var b=[...document.querySelectorAll('#%s .cpt-switch-btn')]
       .find(x=>x.textContent.trim()===%s); if(b){b.click(); return 'ok';} return 'missing';})()",
       sw, shQuote(label, type = "cmd")))
  Sys.sleep(3); wait_idle()
}
sub_tab <- function(pane, label) {
  jsv(sprintf("(function(){var a=[...document.querySelectorAll('#shiny-tab-%s .nav-tabs a')]
       .find(x=>x.textContent.trim()===%s); if(a){a.click(); return 'ok';} return 'missing';})()",
       pane, shQuote(label, type = "cmd")))
  Sys.sleep(4); wait_idle()
}

# The dependency layer needs a cohort before it means anything. Neuroblastoma
# on RNAi is the worked example the README and the report screenshots follow.
# selectInput is a plain <select>; Cancer Subtype is a shinyWidgets pickerInput,
# which is a styled dropdown over a hidden <select>, so both are set through
# the underlying element and the change event is raised by hand.
set_cohort <- function(dataset = "RNAi", subtype = "Neuroblastoma") {
  say("  cohort ->", dataset, "/", subtype)
  jsv(sprintf("(function(){
      var d = document.getElementById('deps-dataset');
      if (d) { d.value = %s; $(d).trigger('change'); }
      var p = document.getElementById('deps-cancer_subtypes');
      if (!p) return 'no subtype widget';
      $(p).selectpicker('val', %s);
      $(p).trigger('change');
      return $(p).val();
    })()", shQuote(dataset, type = "cmd"), shQuote(subtype, type = "cmd")))
}

cat("Capturing from", url, "\n")
set_size(1100)
invisible(b$Page$navigate(url))
invisible(b$Page$loadEventFired(wait_ = TRUE))
Sys.sleep(12)
install_idle_hook()
wait_idle()

shoot("01_home.png")

nav("gene_effect")
set_cohort()
# The opening state, with no layer switched on and every filter disclosure
# opened, so one image carries the whole filter vocabulary.
js("(function(){[...document.querySelectorAll('#shiny-tab-gene_effect input[type=checkbox]')]
     .filter(c=>['dependency','ligandability','functionality'].includes(c.value))
     .forEach(function(c){ if(c.checked) c.click(); }); return 'ok';})()")
Sys.sleep(2); wait_idle()
say("  disclosures opened:",
    jsv("(function(){var d=[...document.querySelectorAll('#shiny-tab-gene_effect details')];
         d.forEach(function(x){ x.open = true; }); return d.length;})()"))
Sys.sleep(2)
shoot("02_discover_guide.png")
js("(function(){[...document.querySelectorAll('#shiny-tab-gene_effect details')]
     .forEach(function(x){ x.open = false; }); return 'ok';})()")
Sys.sleep(1)

# --- Dependency on its own -------------------------------------------------
set_layers("dependency"); run_analysis()
shoot("03_dependency_summary.png")
pick("table", "Cancer-selective genes");        shoot("04_dependency_ranked.png")
# Shots are full-page, so a table shot and the chart shot that follows it come
# out identical unless the chart underneath is different. Each pair below moves
# one of the two switchers, never neither.
pick("chart", "Subtype vs other lines")
pick("table", "Every gene in this subtype");    shoot("05_dependency_all_genes.png")
pick("chart", "Effect vs significance");        shoot("06_dependency_volcano.png")
pick("table", "Cancer-selective genes")
pick("chart", "Subtype vs other lines");        shoot("07_dependency_vs_others.png")

# --- Ligandability on its own ----------------------------------------------
set_layers("ligandability"); run_analysis()
shoot("08_ligandability_summary.png")
pick("table", "Best probe per protein");        shoot("09_ligandability_best_probe.png")
pick("chart", "Competition ratio")
pick("table", "Every engaged cysteine");        shoot("10_ligandability_engaged.png")
pick("chart", "Engagement vs selectivity");     shoot("11_ligandability_scatter.png")
pick("table", "Best probe per protein")
pick("chart", "Competition ratio");             shoot("12_ligandability_cr.png")

# --- Cysteine function on its own ------------------------------------------
set_layers("functionality"); run_analysis()
shoot("13_cysteine_function_engaged.png")
pick("chart", "ABE versus CBE dropout")
pick("table", "Cysteine atlas");                shoot("14_cysteine_atlas.png")
pick("chart", "Tier composition");              shoot("15_cysteine_tiers.png")

# --- All three, where the joins appear -------------------------------------
set_layers(c("dependency", "ligandability", "functionality")); run_analysis()
shoot("16_all_layers_summary.png")
pick("table", "Dependency + ligandability");    shoot("17_combined_dep_lig.png")
pick("table", "Dependency + cysteine function");shoot("18_combined_dep_cys.png")
pick("table", "All combined");                  shoot("19_all_combined.png")
pick("chart", "Dependency by tier");            shoot("20_chart_dependency_by_tier.png")
pick("chart", "Ligandable volcano");            shoot("21_chart_ligandable_volcano.png")
pick("chart", "Tiers of engaged sites");        shoot("22_chart_tiers_engaged.png")

# --- The other tabs --------------------------------------------------------
nav("gene_tab")
set_gene <- function(sym) {
  jsv(sprintf("(function(){
      var s = document.getElementById('gene-gene');
      if (!s || !s.selectize) return 'no widget';
      var sel = s.selectize;
      sel.addOption({value: '%s', label: '%s'});
      // Not addItem(v, true): the second argument is `silent`, which suppresses
      // the change event, so the widget shows the new gene while the server
      // still holds the old one.
      sel.setValue('%s');
      return sel.getValue();
    })()", sym, sym, sym))
}
say("  gene ->", set_gene("IGF2BP3"))
Sys.sleep(15); shoot("23_target_igf2bp3.png")

nav("chemistry_tab")
jsv("(function(){var b=document.getElementById('swiss-load_binding'); if(b) b.click(); return 'ok';})()")
Sys.sleep(12); wait_idle()
shoot("24_chemistry_explorer.png")
sub_tab("chemistry_tab", "Protein Binding")
shoot("25_chemistry_protein_binding.png")

nav("about_tab"); shoot("26_about.png")

# --- The generated report, rendered outside the app ------------------------
# cpt_report_gene_dependency() writes standalone HTML, so it is shot from the
# file rather than from a tab. Sourcing the app's R/ is enough to call it.
report <- tryCatch({
  op <- options(warn = -1)
  on.exit(options(op), add = TRUE)
  for (f in list.files("R", pattern = "[.][Rr]$", full.names = TRUE)) {
    try(suppressPackageStartupMessages(source(f)), silent = TRUE)
  }
  cpt_report_gene_dependency("IGF2BP3", "Neuroblastoma", "RNAi",
                             output_file = "igf2bp3_report.html",
                             output_dir = file.path(tempdir(), "cptreport"),
                             project_root = ".")
}, error = function(e) { say("  report skipped:", conditionMessage(e)); NULL })
if (!is.null(report)) {
  invisible(b$Page$navigate(paste0("file://", normalizePath(report))))
  invisible(b$Page$loadEventFired(wait_ = TRUE))
  Sys.sleep(15)
  h <- tryCatch(jsv("Math.max(document.body.scrollHeight, document.documentElement.scrollHeight)"),
                error = function(e) 1200)
  set_size(max(800L, min(as.integer(h) + 40L, 9000L)))
  Sys.sleep(2)
  path <- file.path(outdir, "27_report_igf2bp3_html.png")
  b$screenshot(filename = path, scale = SCALE)
  cat(sprintf("  %-34s %s\n", "27_report_igf2bp3_html.png",
              format(file.size(path) / 1024, digits = 4)))
}

cat("Done ->", normalizePath(outdir), "\n")
