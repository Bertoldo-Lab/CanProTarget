# R environment

The project does not use renv: a lockfile in the project root breaks the
shinyapps deploy, so it was removed rather than left disabled. This file is
the substitute. It records the exact versions the released build was
developed, benchmarked and deployed against, so a reader can reconstruct the
environment without a lockfile.

Captured 19 September 2026 on the machine that produced the
v1.0.0 release, the 19 September 2026 benchmark run and the screenshots.

## Platform

| | |
|---|---|
| R | R version 4.5.1 (2025-06-13) |
| Platform | aarch64-apple-darwin20 |
| CRAN mirror | https://cran.rstudio.com/ (set in `.Rprofile`) |

## Packages the application loads

Every package reached by `library()` or `::` from `app.R` and `R/`.

| Package | Version |
|---|---|
| `shiny` | 1.11.1 |
| `shinydashboard` | 0.7.3 |
| `DT` | 0.33 |
| `dplyr` | 1.2.1 |
| `magrittr` | 2.0.5 |
| `tidyr` | 1.3.2 |
| `ggplot2` | 4.0.3 |
| `plotly` | 4.11.0 |
| `tibble` | 3.3.1 |
| `writexl` | 1.5.4 |
| `shinyjs` | 2.1.0 |
| `shinycssloaders` | 1.1.0 |
| `shinyWidgets` | 0.9.0 |
| `colourpicker` | 1.3.0 |
| `ggrepel` | 0.9.6 |
| `gtools` | 3.9.5 |
| `jsonlite` | 2.0.0 |
| `readr` | 2.2.0 |
| `rmarkdown` | 2.31 |
| `yaml` | 2.3.10 |
| `htmlwidgets` | 1.6.4 |
| `remotes` | 2.5.0 |
| `cdsrmodels` | 0.1.0 |

`cdsrmodels` is not on CRAN. `R/functions.R` calls its `lin_associations()`
for the per-gene linear-association p-values. It imports limma, WGCNA, ashr,
corpcor, data.table, gausscov, plyr, ranger, tidyverse and useful, so it pulls
a large stack behind it — `Hmisc` arrives through WGCNA. limma is a
Bioconductor package, which `install.packages()` cannot reach:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("limma")
remotes::install_github("broadinstitute/cdsrmodels")
```

`limma` resolves to 3.64.3 in this environment.

## Packages the rebuild scripts and the release tooling need

Not needed to run the app. `readxl` is used by the `docs/scripts/`
preprocessing that rebuilds the RDS files from the source workbooks;
`chromote` drives `docs/scripts/capture_screenshots.R`; `rsconnect` deploys.

| Package | Version |
|---|---|
| `readxl` | 1.4.5 |
| `chromote` | 0.5.1 |
| `rsconnect` | 1.8.0 |

## Installing

```r
install.packages(c(
  "shiny", "shinydashboard", "DT", "dplyr", "magrittr", "tidyr", "ggplot2",
  "plotly", "tibble", "writexl", "shinyjs", "shinycssloaders", "shinyWidgets",
  "colourpicker", "ggrepel", "gtools", "jsonlite", "readr", "rmarkdown",
  "yaml", "remotes",
  # rebuild scripts and release tooling, not needed to run the app
  "readxl", "chromote", "rsconnect"
))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("limma")
remotes::install_github("broadinstitute/cdsrmodels")
```

`htmlwidgets` arrives with `plotly` and `DT`, and `Hmisc` with WGCNA; neither
needs a separate install.

## Full session

`analysis/software_benchmark/outputs/session_info.txt` carries the complete
`sessionInfo()` of the benchmark run, including transitive dependencies.
