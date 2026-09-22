# CanProTarget Score (CPT Score) — Methods

> **Status:** Implemented in `R/canprotarget_score.R`. A research prioritization aid — **not** a clinical recommendation.

This document describes exactly how the composite CPT Score is computed, so that it can be audited, cited, and revised transparently.

---

## 1. What it is

The **CPT Score** is a **0–100 composite prioritization metric** for a gene in a cancer subtype. It blends:

1. How essential the gene is in that cancer (DepMap gene effect)
2. How selective that essentiality is vs other cancers
3. Whether the gene has covalent cysteine handles (editing atlas)
4. Optional conservation / ClinVar evidence on those sites
5. ADME/drug-likeness of engaging probes (**implemented, provisional, weight 0 by default**)

It is used in:

| Surface | Function |
|---------|----------|
| Shiny Dependencies table | `cpt_annotate_gene_table()` — CPT_Score / CPT_Rank columns |
| Single-gene radar / sensitivity | `cpt_score_single()`, `cpt_dimension_radar_gg()`, `cpt_weight_sensitivity()` |
| MCP | `canprotarget_score`, `rank_targets`, `assess_target` |
| HTML reports | Gene dependency report embeds score + dimension breakdown |

**Code of record:** `R/canprotarget_score.R`  
**API wrappers:** `R/api_functions.R` (`api_canprotarget_score`, `api_rank_targets`, …)

---

## 2. Design principles

1. **Rank-based, not raw-scale fusion**  
   Dependency and selectivity live on gene-effect scales; cysteine ligandability is already ~0–100. We convert DepMap axes to **percentile ranks within the subtype**, then weight.

2. **Percentiles need a full gene universe**  
   Scoring one gene alone would always place it at the median (50). Single-gene CPT scores rank dependency/selectivity against **all genes in that subtype** (`cpt_subtype_background()`).

3. **Missing evidence is omitted, not zero**  
   NA dimensions are dropped; remaining weights are **renormalized** over active dimensions. Do not compare scores with different `n_dimensions_used` without checking that field.

4. **High composite ≠ automatic high priority**  
   Priority labels are guarded so non-dependencies cannot be labeled “high”.

5. **Configurable weights**  
   Defaults are below. The Shiny **Discover** sidebar exposes sliders for every dimension (**Defaults** / **Equal** presets). Scores, radar, table ranks, and the sensitivity table recompute from the user’s weights. MCP / reports still use platform defaults unless a `weights` argument is passed in code.

---

## 3. Inputs and raw quantities

### 3.1 Dependency matrix

- **CRISPR:** DepMap 23Q4 gene effect (models × genes)  
- **RNAi:** DEMETER2 v6 (same role when selected)  
- **Metadata:** `OncotreeSubtype` per `ModelID`

For subtype \(S\) and gene \(g\):

\[
\mu_{\text{cancer}}(g,S) = \text{mean gene effect of } g \text{ over cell lines in } S
\]

\[
\mu_{\text{other}}(g,S) = \text{mean gene effect of } g \text{ over cell lines not in } S
\]

\[
\text{effect size}(g,S) = \mu_{\text{cancer}}(g,S) - \mu_{\text{other}}(g,S)
\]

- More **negative** \(\mu_{\text{cancer}}\) → more essential  
- More **negative** effect size → more selectively essential in \(S\)

**Minimum subtype size:** scoring requires **≥ 3** cell lines in the subtype (`cpt_subtype_background`).

### 3.2 Working thresholds (raw, not the composite)

These are used for flags, pool selection, and priority guards (same spirit as MCP `query_dependency`):

| Flag | Rule |
|------|------|
| **Dependency** | \(\mu_{\text{cancer}} < -0.5\) |
| **Selective by effect size** (CPT priority guards) | effect size \(< -0.2\) |
| **Selective dependency** (`query_dependency` / HTML report) | effect size \(< -0.1\) **and** one-sided t-test \(p < 0.05\) |

CPT composite uses the first two for **priority demotion**; it does **not** require the t-test for the score itself.

### 3.3 Cysteine editing atlas

From `data/cys_editing_atlas.rds` (Li et al., Nat Chem Biol 2023 ± lab enrichments). Per-gene dimensions use **max** (or any-true) across sites on that gene.

---

## 4. Dimension scores (each 0–100 or NA)

### 4.1 `dependency_strength`

- **Raw:** \(\mu_{\text{cancer}}\) for all genes in subtype  
- **Score:** percentile rank of that gene among all genes, with **lower (more negative) = better**  
- Formula (implementation):  
  \(\text{percentile} = (\text{rank}(-\mu) - 1) / (n_{\text{valid}} - 1) \times 100\)  
  (`ties.method = "average"`; single valid value → 50)

### 4.2 `cancer_selectivity`

- **Raw:** effect size for all genes in subtype  
- **Score:** same percentile construction, **more negative effect size = better**

### 4.3 `cysteine_ligandability`

For gene \(g\), take all atlas rows with matching `gene_symbol`:

| Condition | Score |
|-----------|--------|
| Not in atlas | **NA** (dimension dropped) |
| `ligandability_score` present | **max** of that column (already 0–100) |
| Else `ligandable` TRUE any site | **50** |
| Else `ligandable` all FALSE | **0** |
| Else `functional_ligandable` TRUE | **75** |
| Else `functional_ligandable` all FALSE | **0** |
| In atlas but no usable columns | **NA** |

### 4.4 `conservation`

- Requires atlas column `conservation_score`  
- Max across sites for the gene  
- If values are on 0–1, scaled to 0–100: `max(cons) * 100`  
- Missing column / no sites / all NA → **NA**

### 4.5 `clinical_evidence`

- Requires `clinvar_pathogenic`  
- Any TRUE → **100**; annotated but none pathogenic → **0**  
- Not in atlas or all NA → **NA**

### 4.6 `adme_druggability`

**Implemented, provisional, weighted 0 by default.**

- **Scored for ~4,470 genes.** `cpt_build_adme_gene_scores()` collapses the 10.6M-row
  binding table plus SwissADME into one row per gene; `score_adme_druggability()` looks
  the gene up. Genes with no covalent probe engaging them at **CR >= 4** score `NA`
  (no credible engagement, therefore no ADME claim) rather than 0.
- **Default weight is 0**, so composite scores are byte-identical to the
  pre-implementation baseline. Verified in `tests/test_adme_dimension.R`. The value is
  still reported in `dim_adme`, the radar and the MCP JSON, so it can be inspected
  before it is allowed to count.
- **Raising the weight is a methods decision, not a display setting.** At the proposed
  0.5 the PDAC `rank_targets` order already changes
  (`RPL34 > DNA2 > DCAF13 > SNRPD3` becomes `RPL34 > SNRPD3 > DCAF1 > USO1`), and that
  default governs reports and MCP output for everyone who never touches a slider.
  `CPT_ADME_PROPOSED_WEIGHT` holds the 0.5 candidate.

#### Formula

For each probe, a 0-100 developability score is the mean of four
fragment-appropriate windows (`cpt_window_score`, full marks inside the window,
linear taper to zero):

| Descriptor | Full marks | Zero at | Rationale |
|------------|-----------|---------|-----------|
| `MW` | <= 300 | 500 | Fragment rule of three |
| `Consensus.Log.P` | 0 to 3 | -3 / 6 | Too polar cannot cross membranes; too greasy brings solubility and promiscuity problems |
| `TPSA` | <= 90 | 140 | Veber oral-absorption guidance |
| `Synthetic.Accessibility` | <= 3 | 7 | Ease of synthesis (1 easy, 10 hard) |

The gene score is the **mean** across qualifying probes.

#### Deliberate exclusions

**Reactivity filters are excluded.** `Brenk..alerts` is 1 for **971 of 1000** probes and
that alert *is* the covalent warhead (acrylamide / chloroacetamide). Rewarding fewer
alerts would penalise compounds for being covalent probes, inverting the platform's
purpose. `PAINS..alerts` is near-constant.

**Pass/fail drug-likeness rules are excluded** because they are saturated over a
fragment library (median MW 235) and would contribute a near-constant term:

| Descriptor | Spread across 998 probes |
|------------|--------------------------|
| `Bioavailability.Score` | 993/1000 identical at 0.55 (2 distinct values) |
| `Lipinski..violations` | 991/1000 are zero |
| `GI.absorption` | 989/1000 "High" |
| `Pgp.substrate` | 965/1000 "No" |

Descriptors that do vary usefully: `MW` (285 distinct), `Consensus.Log.P` (210),
`TPSA` (107), `Synthetic.Accessibility` (163).

#### Why mean and not best probe

Taking the **max** across probes is the intuitive choice and is wrong here. It saturates
at 100 for 75.8% of genes and correlates **+0.52** (Spearman) with the number of probes
hitting the gene, so it measures how well-screened a target is rather than how
developable its chemistry is. The mean decouples from coverage (**+0.015** measured) and
yields 957 distinct values instead of 37. `tests/test_adme_dimension.R` asserts
`|rho| < 0.2` so this cannot silently regress.

#### Known limitation

The spread is narrow: IQR ~92.6-98.7, 19.4% at 100. Fragments genuinely are
physicochemically benign, so this axis **shifts scores more than it reorders them**. It
should not be expected to discriminate strongly. An alternative weighting the mean toward
selective probes via `n_targets` spreads slightly wider but mixes pharmacology into a
pharmacokinetics axis; left unimplemented pending review.

#### Open before the weight is raised

Standard drug-likeness filters carry almost no signal over this probe library, which is
fragment-scale (median MW 235):

| Descriptor | Spread across 998 probes |
|------------|--------------------------|
| `Bioavailability.Score` | 993/1000 identical at 0.55 (2 distinct values) |
| `Lipinski..violations` | 991/1000 are zero |
| `GI.absorption` | 989/1000 "High" |
| `Pgp.substrate` | 965/1000 "No" |
| `PAINS..alerts` | mostly zero (2 distinct values) |

A score built on those would be near-constant and would dilute the biological signal.

Worse, `Brenk..alerts` is **1 for 971 of 1000 probes**, and that alert *is* the covalent
warhead (acrylamide / chloroacetamide). Treating fewer structural alerts as better would
systematically penalise compounds for being covalent probes, which inverts the platform's
purpose.

Descriptors that do vary usefully: `MW` (285 distinct), `Consensus.Log.P` (210), `TPSA`
(107), `Synthetic.Accessibility` (163). A fragment-appropriate window (rule of three,
MW <= 300 and LogP <= 3) is a better yardstick than Lipinski here. Probe selectivity via
`n_targets` from our own binding table is arguably the strongest available signal and is
not a SwissADME field at all.

#### Open before activation

1. Aggregation across probes: best probe, mean, or covalent-only? EGFR is hit by 570
   probes and VCP by 966, and the choice encodes different scientific claims.
2. Default weight. The dimension is computed but weighted 0, so it does not
   move any score today. That is deliberate: at the proposed 0.5, about 9% of
   total weight, a stand-in ADME dimension was already enough to move VCP above
   KRAS for top rank, so switching it on would reshape reports and MCP output
   for users who never touch a slider.
3. Whether reactivity filters are excluded explicitly, and documented as such.

---

## 5. Composite formula

Default weights (`CPT_DEFAULT_WEIGHTS`):

| Dimension | Weight | Rationale (heuristic) |
|-----------|--------|------------------------|
| `dependency_strength` | **3.0** | Core: is it essential? |
| `cancer_selectivity` | **3.0** | Core: is it selective? |
| `cysteine_ligandability` | **2.0** | Covalent handle evidence |
| `conservation` | **1.5** | Functional importance of Cys |
| `clinical_evidence` | **1.0** | Pathogenic site annotations |
| `adme_druggability` | **0.0** | Provisional developability. Computed but weighted 0 so composites match the pre-implementation baseline; `CPT_ADME_PROPOSED_WEIGHT` = 0.5 |

For gene \(g\), let \(D_g\) be the set of dimensions with non-NA scores, and \(w_d\) the default weight.

**Requirement:** \(|D_g| \ge 2\) (`min_dimensions = 2`), else composite = NA.

\[
w'_d = \frac{w_d}{\sum_{d \in D_g} w_d}
\qquad
\text{CPT}(g) = \sum_{d \in D_g} w'_d \cdot s_{g,d}
\]

Rounded to **2 decimals** for display. Rank = dense-style order by descending CPT among the table being ranked (`ties.method = "min"` where applicable).

### Example (schematic)

Active dims only: dep=90, sel=80, lig=70 (weights 3, 3, 2).

\[
w' = (3/8,\ 3/8,\ 2/8) \Rightarrow \text{CPT} = 0.375\cdot90 + 0.375\cdot80 + 0.25\cdot70 = 81.25
\]

ADME is NA for genes with no probe at CR >= 4, and for those genes its weight is renormalised away entirely. Where it is present, the default weight of 0 makes its contribution exactly zero.

---

## 6. Priority labels (heuristic bands + guards)

### 6.1 Score bands (`cpt_priority_label`)

| Composite CPT | Band |
|---------------|------|
| ≥ 70 | `high` |
| ≥ 40 and &lt; 70 | `moderate` |
| &lt; 40 | `low` |
| NA | `unscored` |

### 6.2 Guards (`cpt_apply_priority_guards`)

Applied after the band, using **raw** thresholds:

1. If **not** a dependency (\(\mu_{\text{cancer}} \ge -0.5\)) and band is `high` → demote to **`moderate`**  
2. If **not** a dependency **and** not selective by effect size (ES ≥ −0.2) and band is `moderate` → demote to **`low`**

**Consequence:** non-dependencies are **never** labeled high priority, even if percentile ranks look strong. This is intentional fairness (see e.g. BRCA1-style cases where mid-pack non-essentials can still get middling–high percentiles).

---

## 7. How different call paths use the score

### 7.1 Single gene — `cpt_score_single` / MCP `canprotarget_score`

1. Build full-subtype background (all genes’ means + effect sizes)  
2. Percentile dep + sel against that background  
3. Cysteine dims for that gene only  
4. Weighted composite (min 2 dims)  
5. Priority + guards + caveats + interpretation string  

### 7.2 Shiny results table — `cpt_annotate_gene_table`

1. Percentiles for dep/sel computed on the **full analysis background** gene table  
2. Cysteine dims only for genes in the **filtered** results table (speed)  
3. Composite per row; **CPT_Rank** is rank **within the filtered table**, not genome-wide  

### 7.3 Ranking targets — `cpt_rank_targets` / MCP `rank_targets`

1. Full-background percentiles for all genes  
2. **Candidate pool** (size `pool`, default 100):  
   - Prefer genes with \(\mu_{\text{cancer}} < -0.5\), ordered by most negative effect size  
   - If fewer than 10 such dependencies, fall back to all genes ordered by effect size  
3. Score pool with full CPT dimensions  
4. Optional: `require_ligandable` drops NA ligandability  
5. If ≥ 10 scored dependencies remain, **restrict final ranking to dependencies only**  
6. Sort by CPT descending (then effect size), return top `n`  

**Caveat:** rank_targets is **not** an exhaustive genome-wide CPT sort of all ~18k genes; it re-ranks a selective-dependency pool. For a specific gene’s score, use `canprotarget_score` / `assess_target`.

### 7.4 Site CPT (`cpt_rank_site_targets` / MCP `rank_site_targets`)

Gene-level CPT is unchanged. **Site CPT** applies the same weighted composite to one engaged cysteine:

- Dependency and selectivity percentiles remain **gene-level** (the residue does not have its own DepMap score)
- Ligandability, conservation and ClinVar are taken from **that residue** in the Cys_editing atlas, not the gene-wise max
- ADME, when present, is still the gene-level probe developability score and stays at weight 0 by default
- Missing site dimensions are dropped and remaining weights renormalised

**Evidence tiers** (exact-site accounting):

| Tier | Meaning |
|------|---------|
| 1 | In atlas, functional, and atlas-ligandable |
| 2 | In atlas and functional, without that extra ligandability flag |
| 3 | In atlas, tested, not functional in the Cys_editing contexts |
| 4 | Absent from the atlas (untested). Missing coverage, not a negative result |

MCP `rank_site_targets` orders by **tier first**, then Site CPT within the tier, so a high-scoring untested site cannot outrank a functional Tier 1 or 2 site. `rank_targets` remains the gene-level ranking and may attach `best_engaged_site` as annotation only.

**Gene pool.** Default is the CPT dependency ranking (subtype mean < −0.5). To match the manuscript RNAi discovery filters, pass `effect_size_max = -0.1`, `p_max = 0.05`, `exclude_common_essentials = TRUE` (whole-matrix mean > −0.5), `min_cr = 4`, `max_targets = 20`, and `dataset = "RNAi"`. `max_targets` is applied to every CR ≥ 4 probe **before** the display cap, so a site is not dropped just because its hottest probes are promiscuous.

Engagement is a **CR ≥ 4** row in the bundled chemoproteomic table (six CysDB-indexed studies). That is not a live CysDB query.

---

## 8. User-chosen weights (Shiny)

| Control | Location | Effect |
|---------|----------|--------|
| Six sliders (0–10) | Dependencies left panel | Relative weights for composite |
| **Defaults** | Same panel | Restore `CPT_DEFAULT_WEIGHTS` (3 / 3 / 2 / 1.5 / 1 / **0**) |
| **Equal** | Same panel | All weights = 1 |
| **Reset Controls** | Same panel | Also restores default weights |

Implementation: `cpt_user_weights()` in `R/dependencies_module.R` → `cpt_coerce_weights()` → `cpt_annotate_gene_table(..., weights=)`, `cpt_score_single(..., weights=)`, `cpt_weight_sensitivity(..., weights=)`.

**Note:** Changing weights re-ranks the cancer-selective gene table without re-running limma. Dimension *scores* (percentiles / atlas) stay the same; only the composite blend changes.

## 9. Weight sensitivity (exploratory)

`cpt_weight_sensitivity()` multiplies each active weight by factors (default 0.5, 1, 2) and recomputes the composite. In the UI it starts from **your current slider weights**, not only the product defaults.

---

## 10. What is *not* in the score (yet)

| Item | Status |
|------|--------|
| ADME / SwissADME of probes | **Implemented** for ~4,470 genes, provisional formula, weight 0 by default — see 4.6 |
| limma / FDR from Shiny volcano | Used in UI analysis path; **not** inside CPT composite |
| Precomputed subtype TSV stats | Optional speed path for analysis; CPT background recomputes means from matrix |
| Clinical trial counts | Roadmap only |
| Absolute “druggability probability” | Not claimed — relative prioritization only |

---

## 11. Caveats for interpretation and publication

1. **Scores are subtype-relative.** A 90 in SCLC is not comparable to a 90 in melanoma without context.  
2. **`n_dimensions_used` matters.** Genes with only dep+sel vs genes with full cysteine evidence are not strictly commensurate.  
3. **Atlas coverage is uneven.** Genes absent from the cysteine atlas drop ligandability/conservation/ClinVar.  
4. **Common essentials** can score high on dependency strength but should lose on selectivity; still inspect both dims.  
5. **Priority labels are UX heuristics**, not validated cutoffs.  
6. Always cite methods + data versions (`data/data_versions.yaml` / report provenance).

---

## 12. How to recompute / audit

```r
# From project root
source("R/app_helpers.R")
source("R/api_functions.R")
source("R/canprotarget_score.R")

data_env <- list(
  data_dir = normalizePath("data"),
  crispr_matrix = readRDS("data/CRISPRGeneEffect_23Q4_clean.rds"),
  rnai_matrix = readRDS("data/d2_gene_effect_headers_refined.rds"),
  cancer_model_data = readRDS("data/cancer_model_data.rds"),
  cys_atlas = readRDS("data/cys_editing_atlas.rds")
)

# Single gene
cpt_score_single("TXN", "Small Cell Lung Cancer", "CRISPR", data_env)

# Ranked list
cpt_rank_targets("Small Cell Lung Cancer", "CRISPR", n = 15L, data_env = data_env)
```

Unit / workflow tests: `tests/test_cpt_agent_workflows.R`, `tests/test_finalize_today.R`.

---

## 13. Open methodological questions

The score is a transparent heuristic, not a fitted model. These are the choices
most open to revision, and the reasons they are worth revisiting:

- **Weights.** The defaults (3/3/2/1.5/1, with ADME at 0) encode a judgement
  that dependency and selectivity matter most. They are not data-driven, and equal weighting or
  a fitted alternative would be a defensible substitute.
- **Dependency gate.** \(\mu < -0.5\) demotes priority. The same cut is applied
  to CRISPR and RNAi, whose effect distributions differ.
- **Selectivity significance.** CPT uses effect size alone for the priority
  guards; it does not require \(p < 0.05\).
- **Ligandability fallbacks.** Where a continuous competition ratio is missing,
  binary evidence maps to fixed values (50 / 75 / 0).
- **Ranking pool.** `rank_targets` scores a candidate pool rather than the full
  genome, which bounds compute but also bounds recall.
- **ADME.** Probe quality is not yet a dimension; if ADME activates, whether to
  use the best probe, the mean, or only covalent probes is unresolved.

---

## 14. Related docs

- [DATA_PROVENANCE.md](DATA_PROVENANCE.md) — data sources  
- [AGENTS.md](../AGENTS.md) — agent-facing interpretation notes  
- Implementation: `R/canprotarget_score.R`
