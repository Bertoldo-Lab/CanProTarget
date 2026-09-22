# Data provenance (RDS inputs to CanProTarget)

The **Shiny app only reads compiled `.rds` files** (and plain-text subtype lists, TSV precomputes, etc.). It does **not** load DepMap CSV exports, `swissadme.csv`, `Model_v2.csv`, retracted `Model.csv`, or Excel chemoproteomics sources at runtime.

To **rebuild** those RDS files from the original downloads, use the **reference scripts under `docs/scripts/`** from the **project root** (`cell_cpt/`).

---

## Runtime files under `data/` (what the app expects)

| File | Role |
|------|------|
| `CRISPRGeneEffect_23Q4_clean.rds` | CRISPR gene-effect matrix (models × genes). |
| `d2_gene_effect_headers_refined.rds` | RNAi (D2) gene-effect matrix. |
| `cancer_model_data.rds` | DepMap model metadata (lineage, `ModelID`, `OncotreeSubtype`, …). |
| `cancer_subtypes_CRISPR.txt`, `cancer_subtypes_RNAi.txt` | Allowed subtype labels for the Discover subtype picker. |
| `swissadme_preprocessed.rds` | Chemistry table with `probe_name` and descriptors (e.g. canonical SMILES). |
| `protein_binding_lookup_preprocessed.rds` | Long chemoproteomics binding table (CR == 0 omitted; includes `gene_name_key`). |
| `precomputed_effectsizes/*.rds` | Per-subtype precomputed dependency statistics (fast path). **Not tracked in git** — build locally, see below. |
| `*_modelids.rds` | Optional sidecars (model IDs = matrix row names); used only by the offline subtype-list step. |

---

## Column / structure reference

### Gene effect matrices (`CRISPRGeneEffect_23Q4_clean.rds`, `d2_gene_effect_headers_refined.rds`)

| Part | Meaning |
|------|--------|
| **Object type** | `matrix` (numeric) |
| **Row names** | DepMap **`ModelID`** strings, one row per cell line / model in the screen. |
| **Column names** | Gene symbols (Hugo-style), one column per gene. |
| **Cell values** | Dependency scores (negative values indicate stronger dependency / “essentiality” in that model). Same convention as DepMap **Gene Effect** or **DEMETER2** downloads for the file you used. |

The app aligns rows with `cancer_model_data.rds` by **`ModelID`** (see `R/functions.R`).

---

### Cancer model metadata (`cancer_model_data.rds`)

| Part | Meaning |
|------|--------|
| **Object type** | `data.frame` (from DepMap 23Q4 **`Model_v2.csv`**. Current Portal **`Model.csv`** is retracted.) |
| **Key columns used in code** | See below. Additional columns from DepMap are preserved but need not be present for every analysis path. |

| Column | Meaning |
|--------|--------|
| **`ModelID`** | Unique model identifier; must match **row names** of the gene-effect matrix. |
| **`OncotreeSubtype`** | Oncotree cancer subtype label (used to define “cancer” vs “other” cell lines in analyses). |
| **`OncotreePrimaryDisease`** | Broad disease bucket; values **Non-Cancerous** identify non-tumor controls (used for `NonCancer_Avg` and `pval_vs_NonCancer`). |

Other DepMap columns often present (e.g. `PatientID`, `CellLineName`, `OncotreeLineage`, `SangerModelID`, …) carry lineage and sample identity; see the DepMap **Model** file documentation for the release you imported.

---

### Subtype lists (`cancer_subtypes_CRISPR.txt`, `cancer_subtypes_RNAi.txt`)

| Part | Meaning |
|------|--------|
| **Format** | Plain text; **one OncotreeSubtype string per line** (exact match to `cancer_model_data$OncotreeSubtype`). |
| **Lines starting with `#`** | Treated as comments and ignored. |
| **Purpose** | Populate the Discover subtype picker; lists are filtered by feasibility and precomputed availability in `cancer_subtype_choices_build()` (`R/functions.R`). |

---

### Chemistry table (`swissadme_preprocessed.rds`)

| Part | Meaning |
|------|--------|
| **Object type** | `data.frame` |
| **Origin** | Built offline from **`swissadme.csv`** (SwissADME export) with optional **Backus** link join (`docs/scripts/reference_swissadme_from_csv.R`). |

| Column | Meaning |
|--------|--------|
| **`probe_name`** | Canonical probe label used across the app (after `simplify_probe_name()` at load time). May come from **`Molecule`** or from the Backus link’s **`updated_compound_names`**. |
| **`Molecule`** | May appear if the Backus join kept a SwissADME molecule column; the Chemistry tab export can drop it for display. |
| **`Canonical.SMILES`** (or similar) | SwissADME canonical SMILES (column name may vary slightly; code matches any column whose name contains `canonical` and `smiles`, case-insensitive). |
| **Other columns** | SwissADME physicochemical / ADME descriptors as in the SwissADME web export (e.g. **LogP**, **TPSA**, `ESOL` classes, **LogS**, etc.). Names often use **dots** instead of spaces (`check.names = TRUE` on CSV import). |
| **`compound_names`** / **`compound.names`** | If present, extra compound name fields (see `order_swissadme_display_columns()` in `R/functions.R`). |

Exact descriptor set depends on your **swissadme.csv** export version.

---

### Protein binding lookup (`protein_binding_lookup_preprocessed.rds`)

| Part | Meaning |
|------|--------|
| **Object type** | `data.frame` (long format: one row per probe–protein–cell-line–dataset observation). |
| **Build** | `docs/scripts/preprocess_protein_binding.R`: rows with **`CR == 0`** removed; **`gene_name_key`** added for case-insensitive matching. |
| **Size as bundled** | 10,588,541 rows × 11 columns: 998 probes, 8,425 genes, 37,774 cysteine sites. |
| **Composition** | Aggregate of **six** cysteine-reactive probe studies across **seven** cell lines, not a single screen (see table below). |
| **ID mapping** | All 8,512 distinct `proteinid` values resolve through `data/raw/id_mapping.tsv` (8,707 rows), so unmapped protein IDs are not a source of gene loss. |

#### Source studies pooled into the `Dataset` column

| `Dataset` value | Probes | Genes | Cell lines |
|-----------------|--------|-------|------------|
| `kuljanin_gygi_ligandable` | 858 | 6,956 | HCT116, HEK293T, PaTu-8988T |
| `backus_cravatt_ligandable` | 86 | 6,075 | MDA-MB-231, Ramos |
| `vinogradova_cravatt_ligandable` | 21 | 6,461 | T Cell |
| `yang_wang_ligandable` | 19 | 3,573 | Ramos |
| `cao_backus_ligandable` | 10 | 6,170 | HEK293T |
| `yan_backus_ligandable` | 1 | 4,402 | Jurkat |
| (unlabelled) | 3 | 1,799 | — |

`ligandable` is `"yes"` on 3,421,577 rows and `NA` on 7,166,964; `NA` means the
source sheet carried no call for that site, which CPT scoring treats as missing
evidence rather than as a negative.

**Coverage is wild-type and per-cysteine.** Genes absent from this table were not
detected in these screens. Mutation-created cysteines cannot appear at all:
the screens profile the wild-type cysteinome, so a cysteine that exists only
because of a somatic mutation, such as KRAS G12C, is outside their scope.

| Column | Meaning |
|--------|--------|
| **`probe_name`** | Chemical probe identifier (Chemoproteomic naming; **ACRYL_** / **CL_** / **OTHER_** rules applied in preprocessing). |
| **`proteinid`** | Protein identifier from the ligandable sheet (e.g. UniProt-style). |
| **`gene_name`** | Mapped gene symbol via **`id_mapping.tsv`**. |
| **`gene_name_key`** | `tolower(gene_name)` for fast lookups. |
| **`CR`** | **Competition ratio**: chemoproteomic engagement score for that probe-target pair (higher = stronger binding signal in the assay). |
| **`n_targets`** | Count of distinct **`proteinid`** with **CR ≥ 4** for that probe (selectivity context). |
| **`cysteineid`** | Residue / site identifier from the ligandable sheet. |
| **`ligandable`** | Whether the site is classified ligandable in the source sheet (e.g. `"yes"` / `"no"`). |
| **`Dataset`** | Experimental dataset label from **Compound Keys** (e.g. cell line panel). |
| **`Cell_Line`** | Cell line label from **Compound Keys**. |
| **`SMILES`** | Structure string; after offline preprocessing, may be **overwritten** with SwissADME **canonical** SMILES per `probe_name` when `swissadme_preprocessed.rds` exists (`overlay_canonical_smiles_from_swissadme()` in `R/functions.R`). |

---

### Optional sidecar: `*_modelids.rds` (e.g. `CRISPRGeneEffect_23Q4_clean_modelids.rds`)

| Part | Meaning |
|------|--------|
| **Object type** | `character` vector |
| **Content** | Same strings as **`rownames()`** of the corresponding gene-effect matrix (DepMap **ModelID**). |
| **Use** | Offline only: avoids reloading the full matrix when regenerating **`cancer_subtypes_*.txt`**. |

---

### Precomputed per-subtype statistics (`precomputed_effectsizes/*.rds`)

| Part | Meaning |
|------|--------|
| **Object type** | Data frame stored as RDS. The app prefers `.rds` and falls back to a `.tsv` of the same name, so an older build still works: see `cpt_effectsize_path()` in `R/functions.R`. RDS replaced TSV because it is 850 MB against 279 MB and three to sixteen times faster to read. |
| **Content** | Same rows as **`ge_analysis()`** output **`all_gene_ge_df`** in `R/functions.R` for one **OncotreeSubtype**, produced by **`docs/scripts/precompute_effectsizes_all_subtypes.R`** or by the in-app pipeline. |
| **Naming** | `CRISPR_<Sanitized_Subtype>.rds` or `RNAi_<Sanitized_Subtype>.rds` (sanitization matches `sanitize_subtype()` in `R/dependencies_module.R`). |

| Column | Meaning |
|--------|--------|
| **`gene_name`** | Gene symbol. |
| **`EffectSize`** | Limma-style contrast for **cancer of interest vs rest** (from `cdsrmodels::run_lm_stats_limma` via `run_lm_ge()`); more negative = stronger selective dependency in the subtype. |
| **`Avg`** | Mean effect across **all** cell lines in the matrix (used for “common essentiality” filters). |
| **`Cancer_Avg`** | Mean **gene effect** in cell lines whose **`OncotreeSubtype`** is the selected cancer subtype. |
| **`NonCancer_Avg`** | Mean in lines with **`OncotreePrimaryDisease == "Non-Cancerous"`**. |
| **`Other_Avg`** | Mean in all lines **except** those in the selected subtype (cancer + non-cancer “other” pool). |
| **`p_value`** | Linear-association p-value from **`lin_associations`** (`lin_ass_pval()`), merged with limma results. |
| **`adj.P.Val`** | FDR-adjusted p-value from limma (when present in limma output). |
| **`P.Value`** | Raw limma p-value (if present in upstream output). |
| **`neg_log10_p_value`** | `-log10(p_value)` for volcano plots (may be NA if `p_value` invalid). |
| **`pval_vs_Other_Avg`** | One-sided **t-test** p-value: selected subtype vs **all other** lines (non-subtype pool). |
| **`pval_vs_NonCancer`** | One-sided **t-test** p-value: subtype vs **non-cancerous** lines only. |
| **`pval_vs_OtherCancers`** | One-sided **t-test** p-value: subtype vs **other tumor** lines (excluding non-cancer). |

Not every column is guaranteed to be non-NA for every subtype (depends on sample counts ≥3 per group in `ge_analysis()`).

**Optional extra columns:** The merged limma table may include additional fields from **`cdsrmodels::run_lm_stats_limma`** (e.g. moderated **t**, **logFC**, **B**, confidence intervals) depending on package version; the UI prioritizes **`EffectSize`**, **`p_value`**, **`adj.P.Val`**, and the **`Avg`** / **`Cancer_Avg`** / contrast columns above.

---

## How each RDS was produced (old / source files)

### Functional-cysteine atlas (`cys_editing_atlas.rds`)

| Part | Meaning |
|------|--------|
| **Object type** | Compact `data.frame`, one row per tested cysteine site. |
| **Origin** | Public [`cravattlab/Cys_editing`](https://github.com/cravattlab/Cys_editing) repository, commit `89bc6a268ea4385a687f2e193046b4486edd0185`. |
| **Source objects** | `Cys_dropout.Rdat`, `ceres_used.Rdat`, and `AF_RSA.Rdat` under `Part5_global_analysis/Rdat/`. |
| **Build** | `Rscript docs/scripts/import_cys_editing_atlas.R /path/to/Cys_editing` |
| **Study** | Li et al., *Nature Chemical Biology* (2023), DOI `10.1038/s41589-023-01428-w`. |
| **License/provenance** | Source code is MIT licensed; the license is preserved under `docs/licenses/`. Scientific use must cite the study and source repository. |

The importer reproduces the publication code's within-gene Benjamini–Hochberg
correction and its cell-context consistency filters. A site is marked
`functional` when either editor has mean LFC ≤ -0.6, empirical p < 0.05
(`-log10(p) > 1.3`), and within-gene FDR < 0.1 (`-log10(FDR) > 1`). A site is
marked `ligandable` when `KB_engage > 50`. The source plotting code explicitly
annotates EGFR C797 with `KB_engage = 100`; the imported table retains a
`ligandability_manual_annotation` flag for that record.

UniProt accession and residue checks are joined from the source `AF_RSA` table
using `gene_symbol + cysteine_position`. `residue_mapping_status` distinguishes
source-sequence matches, mismatches, and unmapped sites; it should not be
interpreted as a validation against future UniProt releases.

In **Discover, with the ligandability layer on**, chemoproteomic probe rows are first
annotated at the exact probe site using `UniProt accession + cysteine position`,
with a labelled `gene + position` fallback when accession mapping is unavailable.
The **Functional Cysteine Targets** table additionally links functional sites in
the same dependency gene. Its `cys_probe_site_match` and `cys_relationship`
columns distinguish exact-site evidence from same-gene prioritisation; these two
evidence levels must not be interpreted as equivalent.

### Recovery from the historical ChemProTarget folder

When the original data directory is available, rebuild every core deployable
data product without staging the large raw files inside the repository:

```bash
Rscript docs/scripts/recover_original_data.R /path/to/ChemProTarget/data data
```

This produces both gene-effect RDS matrices and model-ID sidecars,
`cancer_model_data.rds`, both subtype lists,
`protein_binding_lookup_preprocessed.rds`, and
`swissadme_preprocessed.rds`.

### Gene effect matrices (`*_clean.rds`, `d2_*.rds`)

- **CRISPR 23Q4 source:** official `CRISPRGeneEffect.csv` (Public 23Q4). Expected
  shape after Entrez stripping: **1,100 × 18,443**, last gene **`ZZZ3`**, TP53
  present. Fail if the last gene is `TNFRSF10C` (Excel 16,384-column truncation).
  Never open this CSV in Excel.
- **Command:** `Rscript docs/scripts/rebuild_depmap_23q4.R` (or
  `docs/scripts/preprocess_data.R` step 1, which sources that script).
- **RNAi:** D2 refined CSV via `preprocess_data.R` (`d2_gene_effect_headers_refined.csv`).

### Cancer model metadata (`cancer_model_data.rds`)

- **Source:** DepMap 23Q4 **`Model_v2.csv`** (1,921 models × 36 columns). Do **not**
  use current Portal **`Model.csv`**; DepMap retracted it because it includes
  models that do not exist (later freeze, 2,154 models, 6 Acral lines).
- **Command:** `rebuild_depmap_23q4.R`, or `preprocess_data.R` step 3.
- **Acral Melanoma:** metadata lists 4 lines including WM4235 (`ACH-002509`).
  23Q4 CRISPR contains 3 of them. Count n as matrix ∩ subtype, not metadata rows.

### Subtype picker lists (`cancer_subtypes_*.txt`)

- **Source:** `cancer_model_data.rds` + gene-effect RDS + precomputed TSV index (via `precomputed_index_dataframe()`).
- **Command:** `docs/scripts/preprocess_data.R` step 3b, calls `cancer_subtype_choices_build()` in `R/functions.R`.

### SwissADME (`swissadme_preprocessed.rds`)

- **Source:** **`swissadme.csv`** plus optional **`backus_swiss_link.xlsx`** or **`backus_swiss_link.rds`** (probe / Molecule alignment).
- **Helpers:** `docs/scripts/reference_swissadme_from_csv.R` defines `read_swissadme_from_raw_files()` used only offline.
- **Optional:** `Rscript docs/scripts/backus_swiss_link_xlsx_to_rds.R` to refresh `backus_swiss_link.rds` from the xlsx before building SwissADME.
- **Command:** `docs/scripts/preprocess_data.R` step 2.

### Protein binding (`protein_binding_lookup_preprocessed.rds`)

- **Source:** **`table-s2.xlsx`** (Ligandable Dataset + Compound Keys sheets) and **`id_mapping.tsv`** (proteinid → gene); uses **`swissadme_preprocessed.rds`** for canonical SMILES alignment inside the script.
- **Command:** Invoked from `docs/scripts/preprocess_data.R` step 4, which sources **`docs/scripts/preprocess_protein_binding.R`**. That script also `source()`s `R/functions.R` for `overlay_canonical_smiles_from_swissadme()`.
- **Probe naming:** the source columns are renumbered to the scheme SwissADME
  uses. Verified probe by probe against the `Compound Keys` sheet by
  heavy-atom formula: all 998 agree. One boundary previously sent `ACRYL_6`
  onto `ACRYL_4`, merging publication compounds 15 and 14. That is corrected
  here by the CanProTarget authors, not upstream — see `data/raw/README.md`.


### Optional: bulk precomputed effect sizes (`precomputed_effectsizes/*.rds`)

- **Needs:** Gene-effect RDS + `cancer_model_data.rds`.
- **Command:** `Rscript docs/scripts/precompute_effectsizes_all_subtypes.R`  
  The app discovers available subtypes from **the filenames** in `data/precomputed_effectsizes/`.
- **Not distributed.** Both directories are gitignored. Without them the app
  recomputes each subtype in ~20-30 s, which is correct behaviour, not a hang.
- **Filenames must come from `sanitize_subtype()`** (`gsub("[^A-Za-z0-9]+", "_")`,
  then collapse repeats). A cache built with another convention is never found and
  the app silently falls back to recomputation. Generate with the script above
  rather than by hand.
- **Verify any cache you did not generate yourself.** A truncated TSV reads
  without error and simply omits genes. Row count should be close to `ncol()` of
  the matrix, and the last `gene_name` should be a Z gene:

  ```r
  tsv <- readr::read_tsv("data/precomputed_effectsizes/CRISPR_Melanoma.tsv")
  nrow(tsv); tail(tsv$gene_name, 1)
  ```

  A batch recovered in July 2026 stopped at `TNFRSF10C`, dropping 2,508 CRISPR
  genes (TP53, WEE1, TOP2A, VCP and WRN among them) while every value it did
  contain matched the live pipeline to ~1e-9. The cutoff was already in the
  CRISPR matrix (Excel column limit), not in the TSV export.
- **`gene_name` format differs between paths.** The live pipeline emits matrix
  column names carrying the Entrez suffix (`KRAS (3845)`); some external caches
  carry bare symbols. CPT scoring strips the suffix before matching so both work,
  but do not rely on exact string equality across the two paths.

---

## One-shot full rebuild (typical)

From the repo root, with raw files present in `data/` as expected by the scripts:

```bash
# CRISPR 23Q4 + Model_v2 (never Excel; gitignored CSVs)
Rscript docs/scripts/rebuild_depmap_23q4.R

Rscript docs/scripts/preprocess_data.R
```

Order matters inside `preprocess_data.R`: SwissADME RDS before protein binding; model RDS before subtype lists. Do not remake `protein_binding_lookup_preprocessed.rds` unless the chemoproteomics sources changed.

---

## What was removed from the app

- No runtime fallback to `swissadme.csv`, `Model.csv`, `Model_v2.csv`, or gene-effect CSV.
- No reading of CSV log files for precomputed indexing. The app builds the subtype index from the **filenames** in `data/precomputed_effectsizes/`.

All CSV/Excel **ingest** for building RDS lives in **`docs/scripts/`** only.
