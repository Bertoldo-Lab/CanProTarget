# ============================================================
# Script:   deploy_shinyapps.R
# Purpose:  Republish the live app (shinyapps.io application 16957714).
#
#   Usage:  Rscript docs/scripts/deploy_shinyapps.R
#
#   Whatever is checked out is what gets deployed, so a rollback is:
#       git checkout <ref>
#       Rscript docs/scripts/deploy_shinyapps.R
#
# What this script exists to get right
#
#   The file list is explicit. The repo holds ~1.5 GB under data/; the bundle
#   needs ~1.04 GB of it and must exclude the 400 MB source CSV, data/raw/,
#   and any truncated-cache backup.
#
# There is no renv lockfile to work around any more. The project does not use
# renv -- the app runs from the system library -- and rsconnect 1.8 used to
# copy renv.lock into the bundle unconditionally, then refuse to deploy with
# "parseRenvDependencies(): Library and lockfile are out of sync". The script
# moved the lockfile aside for the duration and restored it afterwards, which
# was a trap of its own: a `git add -A` inside that window once committed the
# hold file and untracked the lockfile. Removing renv removed both problems.
# ============================================================

app_dir <- normalizePath(".", mustWork = TRUE)
if (!file.exists(file.path(app_dir, "app.R"))) {
  stop("Run this from the repository root (app.R not found here).")
}

# gene_index.rds is included only when present, so checking out a ref from
# before the Gene tab and redeploying still produces a valid bundle.
app_files <- c(
  "app.R", "R", "www", "inst/report_templates",
  "data/CRISPRGeneEffect_23Q4_clean.rds",
  "data/CRISPRGeneEffect_23Q4_clean_modelids.rds",
  "data/d2_gene_effect_headers_refined.rds",
  "data/cancer_model_data.rds",
  "data/cancer_subtypes_CRISPR.txt",
  "data/cancer_subtypes_RNAi.txt",
  "data/cys_editing_atlas.rds",
  # Both forms of the binding table ship. The factored copy is the one the app
  # reads (555 MB rather than 936 MB in RAM, 1.8 s rather than 22.8 s), but
  # several paths still key on the original: the bindCache signatures here and
  # the MCP worker's SMCL index. Dropping it to save 24 MB on a 555 MB bundle
  # took the ligandability layer offline in production, which is not a trade
  # worth 4%.
  "data/protein_binding_lookup_preprocessed.rds",
  if (file.exists(file.path(app_dir, "data/protein_binding_lookup_factored.rds")))
    "data/protein_binding_lookup_factored.rds",
  if (file.exists(file.path(app_dir, "data/protein_binding_cr4.rds")))
    "data/protein_binding_cr4.rds",
  "data/swissadme_preprocessed.rds",
  "data/data_versions.yaml",
  if (file.exists(file.path(app_dir, "data/gene_index.rds"))) "data/gene_index.rds",
  "data/precomputed_effectsizes"
)
app_files <- app_files[!vapply(app_files, is.null, logical(1))]

missing <- app_files[!file.exists(file.path(app_dir, app_files))]
if (length(missing)) {
  stop("Missing bundle inputs:\n  ", paste(missing, collapse = "\n  "),
       "\nSee data/README.md — the DepMap matrices are not in git.")
}

ref <- tryCatch(system("git rev-parse --short HEAD", intern = TRUE), error = function(e) "unknown")
message("Deploying ", ref, " to application 16957714 ...")

# Guard against a lockfile reappearing: rsconnect bundles it whatever appFiles
# says, and the deploy then fails on a library/lockfile mismatch.
if (file.exists(file.path(app_dir, "renv.lock"))) {
  stop("renv.lock is present. This project does not use renv and rsconnect ",
       "will bundle the lockfile and refuse the deploy. Remove it first.",
       call. = FALSE)
}

rsconnect::deployApp(
  appDir      = app_dir,
  appFiles    = app_files,
  appName     = "CanProTarget",
  account     = "johnpaulong",     # public slug is bertoldolab
  server      = "shinyapps.io",
  appId       = 16957714,
  forceUpdate = TRUE,
  logLevel    = "normal",          # verbose trips an rsconnect/httr2 bug
  lint        = FALSE,
  launch.browser = FALSE
)

message("Deployed ", ref, " -> https://bertoldolab.shinyapps.io/CanProTarget/")
