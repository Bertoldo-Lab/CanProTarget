# Maintenance notes

Working notes for people developing CanProTarget: the history of data problems
and how they were resolved, operational traps, branch conventions and the
issue log. Written for maintainers, not for users — see the
[README](../README.md) for what the project is and how to run it.

Some sections describe states the project has since moved past. They are kept
because the reasoning behind a fix is easy to lose and expensive to rediscover.

---

**CanProTarget** (Cancer chemoproteomics-based target prioritization pipeline) prioritizes **cancer gene dependencies** and links them to **covalent handles** (functional / ligandable cysteines), **chemoproteomic probes**, and **SwissADME** drug-likeness.

- **Shiny app** for interactive analysis and branded HTML reports  
- **MCP server** so AI agents can query the same logic programmatically  
- **CPT Score** (0-100) combining dependency, selectivity, and cysteine evidence  

Original Shiny foundation: John-Paul Ong, `cell_cpt`. Platform extensions (MCP, CPT score, reports, agent docs) are developed here.

Live app: https://bertoldolab.shinyapps.io/CanProTarget/  
Release repository: [Bertoldo-Lab/CanProTarget](https://github.com/Bertoldo-Lab/CanProTarget), branch `main`. License: GNU AGPL v3.0.  
Development history stays in the private repository [Bertoldo-Lab/CanProTarget_dev](https://github.com/Bertoldo-Lab/CanProTarget_dev).

---

## Status of previously open issues

The August 2026 README listed eight open items (local vs live CRISPR mismatch,
`n_targets`, ADME weight, small-n warning, report-test branding, orphan
`precomputed_ligandable/`, uncommitted cache, published truncated cache). Those
items are answered here in the same order. Work that was *not* on that list is
under [Changes in this update](#changes-in-this-update).

| # | Original issue | Status |
|---|----------------|--------|
| 1 | Local CRISPR matrix disagreed with the live app (Acral Melanoma: 984 vs 870 genes; CCND1 present only on the live app; inferred n = 4 vs 3) | **Resolved.** The disagreement was a truncated CRISPR RDS (Excel 16,384-column limit: 1,100 × 16,383, last gene `TNFRSF10C`, no TP53) plus retracted Portal `Model.csv`. Official Public 23Q4 was redownloaded, cleaned in R, and written to RDS (**1,100 × 18,443**, last gene `ZZZ3`). Metadata is **`Model_v2.csv`**. Acral Melanoma CRISPR is **n = 3**: `ACH-002509` / WM4235 is in metadata but not in 23Q4 CRISPR, so n = 4 was never a valid screen count. The live app was republished on 20 August 2026 from this matrix. Do not ingest Portal `Model.csv` or 26Q1 CRISPR. |
| 2 | `n_targets` did not reproduce from the shipped probe table (450 / 998 probes) | **Addressed in the UI; RDS not rebuilt.** Discover’s “max probe targets” now uses the stored `n_targets` column (distinct proteins with CR ≥ 4), not a count of cysteine rows. That column was computed upstream, before gene-mapping losses, so recomputing it from `protein_binding_lookup_preprocessed.rds` still will not match every probe. The chemoproteomics RDS was not regenerated. |
| 3 | ADME CPT dimension provisional, weight 0 | **Open.** Scored for ~4,470 genes and shown in the UI; it does not enter `cpt_score` until the weight is raised. Formula excludes reactivity filters (Brenk alert on 971/1000 probes *is* the warhead). Spread is narrow (IQR 92.6–98.7). See [docs/CPT_SCORE.md](docs/CPT_SCORE.md) §4.6. |
| 4 | No warning on tiny subtypes (Acral Melanoma CRISPR n = 3) | **Open.** The app still allows n = 3 with no UI or report flag for n < 5. |
| 5 | Two report tests fail on branding (`github.com/danielxb` not in the HTML templates) | **Addressed in tests.** Templates still do not print the repo URL in HTML; tests now look for `github.com/Bertoldo-Lab` in the MCP `source=` field. |
| 6 | `precomputed_ligandable/` referenced but never generated | **Resolved.** The unused TSV fast path was removed. The ligandable volcano is always probe results ∩ the dependency table via `cpt_gene_match_key()`. |
| 7 | No precompute cache in git (~20–30 s per subtype) | **Unchanged by design for git.** Caches are gitignored (derived, large, and a truncated copy fails silently). The **live app** now includes the full per-subtype TSVs in the shinyapps bundle, so Discover there is seconds, not limma-on-server. Local clones still build the cache with `docs/scripts/precompute_effectsizes_all_subtypes.R` if they want that speed. |

**Known limitation (not a bug).** KRAS has no ligandable cysteine and no probe rows
here: wild-type KRAS has glycine at codon 12 (so G12C cannot appear in a wild-type
atlas), and KRAS was not detected in the six probe screens (HRAS and NRAS were).
Missing evidence is not a negative result. The same limit is stated in
[docs/DATA_PROVENANCE.md](DATA_PROVENANCE.md).

Subtype names must match OncotreeSubtype exactly, with ≥ 3 cell lines
(`Colon Adenocarcinoma` is valid; `Colorectal Adenocarcinoma` is not).

---

## Changes in this update

Numbered items 1–5 are the data and copy changes relative to the previous
branch. The rest is application behaviour that went out with the same live
republish (20 August 2026).

1. **Full CRISPR 23Q4 matrix.** The previous RDS was the Excel-truncated 16,383-gene
   file. Official Public 23Q4 was redownloaded, cleaned in R (never opened in Excel),
   and written to `data/CRISPRGeneEffect_23Q4_clean.rds`. Shape **1,100 × 18,443**,
   last gene `ZZZ3`, TP53 present. [`docs/scripts/rebuild_depmap_23q4.R`](docs/scripts/rebuild_depmap_23q4.R)
   refuses a last column of `TNFRSF10C`.
2. **`Model.csv` replaced with `Model_v2.csv`.** DepMap retracted Portal `Model.csv`
   because it listed models with no screen data. Acral Melanoma was the example:
   four models in metadata, but WM4235 / `ACH-002509` is absent from 23Q4 CRISPR.
   Subtype *n* is now metadata ∩ matrix (Acral CRISPR = 3).
3. **HTML Gene Dependency Report “Is selective”** uses effect size **< −0.1** (was
   −0.2) in the Metric table, matching the app default. The Yes/No flag is still a
   one-sided t-test on mean difference, not limma / `lin_associations`.
4. **HTML dependency reports use the Group Comparison gene.** Opening a report
   previously always used the top CPT-scored gene. It now uses the gene selected on
   the open Group Comparison subtab.
5. **“Cancer-specific” → “Cancer-selective”.** The Dependencies hit list is labelled
   **Cancer-Selective Genes** (essential in this subtype relative to other cancers,
   not “tumour-only”).

Also in this update:

- **Live app republished** (20 August 2026) with the full CRISPR/RNAi RDS files and
  per-subtype TSVs. Git push does not update shinyapps.io. Verified on
  https://bertoldolab.shinyapps.io/CanProTarget/: Acral Melanoma All Genes includes
  **ZZZ3** (17,787 rows; small-*n* NA drop, not the Excel cutoff); BRAF Compare
  Subtypes Melanoma −1.015 (n = 52) vs Colon Adenocarcinoma −0.282 (n = 45); Load
  Probes, Find Probes, and Load Protein Binding complete on the first click.
- **Graphical abstract** on the Home tab (`www/graphical_abstract.png`); RNAi labelled
  **shRNA**, not miRNA. The figure now shows the new dependency–ligand pairs discussed
  in the manuscript.
- **Run Analysis / Load Probes / Load Protein Binding** complete on the first click.
  `observeEvent` wrapping `eventReactive` on the same button made click 1 return
  `NULL`. Results now sit in `reactiveVal`s and are computed in the observer
  (`R/dependencies_module.R`, `R/swissadme_module.R`).
- **Probe lookups no longer scan ~10.6M rows per click.** `cpt_index_protein_binding()`
  / `cpt_pb_subset()` (`R/app_helpers.R`, built once in `app.R`) index by gene and
  probe. Load Probes filters to cancer-gene keys before isoform/CR scans.
- **All-genes volcano** interactive plotly is capped at 4,000 points (hits kept).
  PNG download is still the full set.
- **Compare Subtypes plot** is `plotly::plot_ly` (wrapped labels, y-axis padded
  around the bars and −0.5) instead of `ggplotly`.
- **HTML reports:** pandoc error 99 came from `includes: in_header: null` and a
  missing logo path. Templates now skip that include and prefer `www/cpt_logo.png`.
  Reports are generated on demand; no sample HTML is committed.

### Data products (for maintainers)

Gene-effect matrices and per-subtype TSVs are **gitignored**. The app reads them at
runtime; rebuild scripts write them.

| Product | Produced by | Notes |
|---------|-------------|--------|
| `data/CRISPRGeneEffect_23Q4_clean.rds` | `docs/scripts/rebuild_depmap_23q4.R` from `CRISPRGeneEffect_23Q4.csv` (`utils::read.csv`; never Excel). Strips DepMap `GENE (Entrez)` suffixes. | Also in the shinyapps bundle. |
| `data/cancer_model_data.rds` and subtype `.txt` lists | Same script, from `Model_v2.csv` | Pickers are Model_v2 ∩ matrix. |
| `data/precomputed_effectsizes/<CRISPR\|RNAi>_<subtype>.tsv` | `docs/scripts/precompute_effectsizes_all_subtypes.R` | If present, Discover skips limma. The app indexes **`.tsv` filenames**, not CSV logs. In the shinyapps bundle as of August 2026. |
| `data/protein_binding_lookup_preprocessed.rds` | Existing preprocess | Independent of the CRISPR rebuild; do not remake it for this update. |

```bash
# Inputs (gitignored) under data/: CRISPRGeneEffect_23Q4.csv, Model_v2.csv
Rscript docs/scripts/rebuild_depmap_23q4.R
Rscript docs/scripts/precompute_effectsizes_all_subtypes.R   # optional, ~1 h
```

Earlier in August 2026 (already on `main`): five DepMap ↔ probe/atlas joins now
use `cpt_gene_match_key()` (`tests/test_gene_key_matching.R`); chemoproteomics
and SwissADME RDS files restored (998 probes × 8,425 genes × 37,774 sites from
**six** studies); ADME CPT dimension implemented at weight 0; AGENTS.md no longer
uses `Colorectal Adenocarcinoma` (1 cell line).

---

## Current status

| Area | Status | Notes |
|------|--------|-------|
| **Discover** | Working | CRISPR 23Q4 + RNAi (DEMETER2); CPT columns, user weight sliders, radar, weight sensitivity, Compare Subtypes, Share This View |
| **Discover: residue evidence** | Working | Cys_editing atlas in repo (Li et al. 2023) |
| **Reports** | Working | Gene dependency + cysteine HTML (download in app; MCP `generate_report`) |
| **MCP / agent API** | Working | 14 tools; see [AGENTS.md](AGENTS.md) and [mcp/README.md](mcp/README.md) |
| **About** | Working | Data versions from `data/data_versions.yaml` |
| **Discover: ligandability** | Working | `data/protein_binding_lookup_preprocessed.rds` bundled |
| **SwissADME** | Working | `data/swissadme_preprocessed.rds` bundled |

Discover, Target and Chemistry run from bundled chemoproteomics / atlas / SwissADME data.
Local clones still download DepMap gene-effect matrices separately (see
[data/README.md](../data/README.md)). The live shinyapps bundle includes those
matrices plus the subtype TSV cache.

### Precomputed subtype cache

`data/precomputed_effectsizes/` is **gitignored**.
A local clone without those TSVs recomputes each subtype (limma + linear
association + three t-tests over ~18,400 genes) in roughly **20–30 seconds**.
That is expected, not a hang.

The **live shinyapps.io app includes the TSV cache in the published bundle**,
so Discover there should return in seconds. That cache is not on
GitHub on purpose: it is derived, large (~850 MB uncompressed), and a truncated
copy fails silently (missing genes, no error).

To build the cache locally:

```bash
Rscript docs/scripts/precompute_effectsizes_all_subtypes.R
```

The app picks the TSVs up from `data/precomputed_effectsizes/` (or
`data/precomputed_effectsizes/`) and skips recomputation. Expect the full run to take
around an hour. `data-assets-v1` no longer ships a stats tarball (deleted 20 Aug 2026).

---

## Features

### Shiny app

| Tab / capability | What it does |
|------------------|--------------|
| **Discover** | Subtype analysis (limma + t-tests), volcano, group comparison, probe overlay when data exist, **CPT Score** table columns, dimension radar, weight sensitivity, **Compare Subtypes** |
| **Discover: residue evidence** | Functional / ligandable cysteine sites; intersect with dependency results |
| **Discover: ligandability** | Probes for entered genes (requires binding RDS) |
| **SwissADME** | ADME properties, BOILED-Egg (WLOGP x, TPSA y), radar (requires SwissADME RDS) |
| **About** | Data versions, workflow, citation |
| **Share This View** | URL bookmark restores tab + subtype |
| **Reports** | Branded HTML with methods notes and data provenance |

### Agent / MCP (14 tools)

Prefer `assess_target` and `rank_targets` for multi-step questions. Full table and examples: **[AGENTS.md](AGENTS.md)**.

```json
{
  "mcpServers": {
    "canprotarget": {
      "command": "mcp/.venv/bin/python",
      "args": ["mcp/canprotarget_server.py"],
      "cwd": "/path/to/CanProTarget"
    }
  }
}
```

Setup: Python 3.10+ with `mcp`, R 4.x+ with `jsonlite`, data under `data/`. Details in [mcp/README.md](mcp/README.md).

### CPT Score

Composite 0–100 prioritization score (rank-based dimensions): dependency
strength, cancer selectivity, cysteine ligandability, conservation, ClinVar.
The ADME dimension is **computed but weighted 0** by default (see Known
limitations). Single-gene scores rank dependency/selectivity against **all genes
in the subtype**. Research aid only, not a clinical recommendation.

**Full methods (weights, percentiles, priority guards, rank_targets pool):** **[docs/CPT_SCORE.md](docs/CPT_SCORE.md)**.

---

## Roadmap (summary)

The development plan is kept outside the public repository.

| Phase | What | Status |
|-------|------|--------|
| 0 | Chemoproteomics / SwissADME files + atlas | Done (bundled RDS; six probe studies) |
| 1 | Cleanup, error UX, data versions, bookmarking | Done on `dev` |
| 2 | CPT Score, heatmap, reports, compare UI (core) | Done on `dev` (ADME dim + list-diff/Venn later) |
| 3 | MCP server + AGENTS.md | Done on `dev` |
| 4 | Performance / precomputation | Live app ships per-subtype TSVs (August 2026); local cache still optional |
| 5 | Paper prep, case studies, landing page | Recovery fixtures in `tests/fixtures/` |
| 6 | Docker, CI/CD, Zenodo | Planned |
| 7 | Post-publication (e.g. SL explorer) | Future |

---

## Quick start (local)

### 1. Install R

Developed with R 4.5+; R 4.6.x works with `Rscript --vanilla`. The project
does not use renv; packages come from the system library.

### 2. Install packages

```r
install.packages(c(
  "shiny", "shinydashboard", "dplyr", "tidyr", "ggplot2",
  "plotly", "DT", "tibble", "writexl", "magrittr",
  "shinyjs", "shinycssloaders", "shinyWidgets", "colourpicker",
  "ggrepel", "gtools", "jsonlite", "readr", "rmarkdown", "yaml", "remotes",
  # rebuild scripts and release tooling, not needed to run the app
  "readxl", "chromote", "rsconnect"
))

# limma is on Bioconductor, so install.packages() cannot reach it, and
# cdsrmodels imports it.
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("limma")

remotes::install_github("broadinstitute/cdsrmodels")
```

`cdsrmodels` is what `R/functions.R` calls for `lin_associations()`. It imports
limma, WGCNA, ashr, gausscov, ranger and tidyverse, so it pulls a large stack;
`Hmisc` arrives through WGCNA and needs no separate line.

The exact versions the v1.0.0 release was built, benchmarked and deployed
against are in [docs/R_ENVIRONMENT.md](R_ENVIRONMENT.md), with the full
`sessionInfo()` in `analysis/software_benchmark/outputs/session_info.txt`.

The GitHub package is `cdsrmodels`, not `cdsr_models`.

### 3. Data

DepMap-style matrices and model metadata live under `data/` (many large RDS files are gitignored). Preprocess from raw downloads:

```bash
Rscript --vanilla docs/scripts/preprocess_data.R
```

Or recover from an old ChemProTarget data folder:

```bash
Rscript --vanilla docs/scripts/recover_original_data.R /path/to/ChemProTarget/data
```

**Missing chemoproteomics / SwissADME:** see **[docs/DATA_PROVENANCE.md](DATA_PROVENANCE.md)**.  
Protein-binding offline build (when you have `table-s2.xlsx` + mapping): `docs/scripts/preprocess_protein_binding.R` (pointer: `R/preprocess_protein_binding.r`).

Data inventory: `data/data_versions.yaml`. Provenance notes: [docs/DATA_PROVENANCE.md](docs/DATA_PROVENANCE.md).

### 4. Run the app

```bash
Rscript --vanilla -e 'shiny::runApp(".", port = 3838, host = "127.0.0.1", launch.browser = TRUE)'
```

### 5. Tests (from project root)

```bash
Rscript --vanilla tests/test_cpt_agent_workflows.R
Rscript --vanilla tests/test_phase1_cleanup.R
Rscript --vanilla tests/test_finalize_today.R
Rscript --vanilla tests/test_reports.R
# MCP (needs mcp/.venv)
mcp/.venv/bin/pytest tests/test_mcp_server.py -v
```

---

## Project structure

```
CanProTarget/
├── app.R                         # Shiny entry + shared data
├── AGENTS.md                     # MCP / agent tool guide
├── R/
│   ├── api_functions.R           # Pure MCP/API queries
│   ├── canprotarget_score.R      # CPT Score + annotate / radar / sensitivity
│   ├── mcp_worker.R              # Warm R worker for MCP
│   ├── report_generator.R        # HTML report renderer
│   ├── pancancer_profile.R
│   ├── theme_canprotarget.R
│   ├── error_codes.R             # CPT-XXXX notifications
│   ├── functions.R               # Analysis helpers (ge_analysis, BOILED-Egg, ...)
│   ├── dependencies_module.R     # Discover tab
│   ├── cys_editing_module.R
│   ├── protein_lookup_module.R
│   ├── swissadme_module.R
│   └── preprocess_protein_binding.r  # pointer → docs/scripts/
├── mcp/                          # Python FastMCP server + R bridge
├── inst/report_templates/        # Rmd + CSS for reports
├── data/                         # Runtime RDS + rebuild sources
│   ├── cys_editing_atlas.rds     # Bundled (conservation / ortholog / ClinVar)
│   ├── protein_binding_lookup_preprocessed.rds   # ligandability + Chemistry
│   ├── swissadme_preprocessed.rds                # SwissADME
│   ├── raw/                      # Backus / SwissADME rebuild inputs (xlsx via LFS)
│   ├── data_versions.yaml        # Provenance manifest (About tab + reports)
│   └── README.md
├── docs/
│   ├── DATA_PROVENANCE.md
│   ├── ARCHITECTURE.md
│   └── scripts/                  # Offline preprocess / recover
├── tests/
│   └── fixtures/                 # Known-target recovery fixtures (KRAS, BRAF, EGFR_797)
└── www/                          # Logo, graphical abstract
```

---

## Data sources

| Dataset | Source | Used by |
|---------|--------|---------|
| CRISPR Gene Effect (23Q4) | [DepMap](https://depmap.org/portal/) | Dependencies, CPT, MCP |
| RNAi Gene Effect (DEMETER2) | [DepMap](https://depmap.org/portal/) | Dependencies, CPT, MCP |
| Cell line metadata | DepMap 23Q4 `Model_v2.csv` (Portal `Model.csv` is retracted) | Subtype grouping |
| Functional cysteine atlas | [Li et al., Nat Chem Biol 2023](https://doi.org/10.1038/s41589-023-01428-w) | Cys tab, CPT, MCP |
| Chemoproteomics (CR) | Six cysteine-reactive probe studies, lab reformatted (see below) | Discover: ligandability, Chemistry |
| SwissADME | [swissadme.ch](http://www.swissadme.ch/) | Chemistry tab |

### Chemoproteomics coverage

`data/protein_binding_lookup_preprocessed.rds` is an aggregate, not a single
screen: 10,588,541 probe-cysteine measurements over **998 probes, 8,425 genes
and 37,774 cysteine sites**, pooled from six studies across seven cell lines.

| Source dataset | Probes | Genes | Cell lines |
|----------------|--------|-------|------------|
| `kuljanin_gygi_ligandable` | 858 | 6,956 | HCT116, HEK293T, PaTu-8988T |
| `backus_cravatt_ligandable` | 86 | 6,075 | MDA-MB-231, Ramos |
| `vinogradova_cravatt_ligandable` | 21 | 6,461 | T Cell |
| `yang_wang_ligandable` | 19 | 3,573 | Ramos |
| `cao_backus_ligandable` | 10 | 6,170 | HEK293T |
| `yan_backus_ligandable` | 1 | 4,402 | Jurkat |

996 of the 998 probes have matching SwissADME descriptors.

**Coverage is per-cysteine and wild-type.** A gene absent from this table was
not detected in these screens; that is not evidence against ligandability.
Mutation-created cysteines (KRAS G12C being the obvious case) cannot appear
here at all, because the screens profile the wild-type cysteinome. See
[docs/DATA_PROVENANCE.md](DATA_PROVENANCE.md).

---

## Biological / methods notes

- **Gene effect**: more negative = more essential. Common working thresholds: mean &lt; -0.5 (dependency), effect size &lt; -0.1 with p &lt; 0.05 (selective). Conventions for prioritization, not clinical cutoffs.
- **CPT Score**: research prioritization aid; dimensions and caveats are listed in JSON / report output.
- **Cysteines**: base-editing dropout (ABE/CBE) + ligandability flags from the atlas.
- **Specialist scope**: prioritization + covalent handles + links out (AlphaFold, UniProt, PDB, DepMap). Not a general science workbench.

---

## Citation

> Ong JP, Martins D, Bertoldo JB. CanProTarget. Zenodo. DOI pending release.

Also credit the original app: **John-Paul Ong**, `cell_cpt`.

---

## Publishing (shinyapps.io)

Live URL: https://bertoldolab.shinyapps.io/CanProTarget/

The app is published manually (`rsconnect`); pushing git does not update it.
The 20 August 2026 bundle includes the full CRISPR/RNAi RDS files and
`data/precomputed_effectsizes/*.tsv`.

### Republishing

One command, from the repository root of a checkout that has the runtime data:

```sh
Rscript docs/scripts/deploy_shinyapps.R
```

Whatever is checked out is what gets deployed, so **rolling back is a checkout
plus the same command**:

```sh
git checkout 13ce167          # the bundle live before the IA/performance work
Rscript docs/scripts/deploy_shinyapps.R
```

The script stages `app.R`, `R/`, `www/`, `inst/report_templates/`, the runtime
RDS files, `data_versions.yaml`, the subtype lists, `gene_index.rds` (when the
checked-out ref has it) and the TSV cache; it excludes `data/raw/`, the 400 MB
source CSV and any truncated-cache backup. Target app id **16957714**, account
`johnpaulong` (public URL slug `bertoldolab`). It uses `logLevel = "normal"`;
`verbose` hits an rsconnect/httr2 bug.

**Why the script rather than a bare `deployApp()` call:** the bundle contents
have to be named explicitly, or rsconnect uploads the whole of `data/`,
including the 400 MB source CSV and the rebuild inputs under `data/raw/`.

The project also carries no `renv.lock`. rsconnect 1.8 copies one into the
bundle *unconditionally* — neither `appFiles` nor `.rscignore` prevents it —
and then aborts with `parseRenvDependencies(): Library and lockfile are out of
sync`. The script used to move the lockfile aside and restore it afterwards,
which created a window in which a `git add -A` could commit the hold file and
untrack the lockfile; that happened. renv is gone, and the script now simply
stops if a lockfile reappears.

Deploying needs the ~1 GB of runtime data that is not in git (see
[data/README.md](../data/README.md)); the script fails early and names anything
absent rather than publishing a broken bundle.

---

## Branch strategy

This repository is the release snapshot. It has one branch, `main`.

Ongoing development, including the full commit history, stays in the private
repository [Bertoldo-Lab/CanProTarget_dev](https://github.com/Bertoldo-Lab/CanProTarget_dev).
Do not push release fixes to personal forks as the source of truth.
`johnpaul-ong/cell_cpt` is the historical Shiny foundation and is not the live project.

### Clone

```bash
git clone https://github.com/Bertoldo-Lab/CanProTarget.git
cd CanProTarget
git lfs pull
```

Git push does **not** update https://bertoldolab.shinyapps.io/CanProTarget/ — that is a manual republish.

### Paper drafts

`paper/` is gitignored. Keep Word/LaTeX manuscripts off GitHub.

### Before changing the release

1. App starts: `Rscript --vanilla -e 'shiny::runApp(".", port=3838, host="127.0.0.1", launch.browser=FALSE)'`
2. Tabs load
3. Tests: `Rscript --vanilla tests/test_site_evidence.R` and `mcp/.venv/bin/pytest tests/test_mcp_server.py -v`
4. README and data/README.md still accurate

---

## Contributing

Science / paper: Bertoldo Lab. Original Shiny app: John-Paul Ong. MCP, CPT Score, reports: this repo.
Development history and pull requests live in `Bertoldo-Lab/CanProTarget_dev`.

See also: [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md), [docs/LOGGING.md](docs/LOGGING.md).
