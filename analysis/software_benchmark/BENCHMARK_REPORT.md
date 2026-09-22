# CanProTarget CPT Score software benchmark

Run date: 19 September 2026. The 3 August 2026 run it replaces is kept in
`outputs_2026-08-03/` and `figures_2026-08-03/`.

## Why this was re-run

The ADME dimension was implemented on 4 August 2026, one day after the first
benchmark ran. That run therefore pinned `adme_druggability` to `NA`, every ADME
weight it recorded was inert, and it could report neither an ADME coverage
figure nor what activating the dimension does. The manuscript reports both, so
its numbers could not have come from the archived run, and the two disagreed on
equal weighting because the archived "equal" scenario equalised five dimensions
while the manuscript equalised six.

Re-running against the release implementation reproduces the manuscript to the
precision it quotes. The archived run was the stale artefact, not the paper.

Building the dimension for real also exposed a defect in the application. The
chemoproteomics tables spell probes with SwissADME's short codes (`AC5`,
`CL174`) and the SwissADME file on disk uses the publication names (`ACRYL_5`,
`CL_64`). The Shiny loader normalised on read; `cpt_build_adme_gene_scores()`
did not, so any caller reading the RDS directly — the MCP worker — joined on the
7 of 998 names that coincide by accident and scored 363 genes instead of 4,468.
The scorer now canonicalises its own inputs. Composite scores are unaffected at
the default weight of 0; the MCP `dim_adme` field was wrong and is now right.

## Outcome

The benchmark supports three claims. First, the independent benchmark
implementation reproduced every one of the 14,877 CPT Scores displayed by the
application across 20 cancer cohorts after two-decimal rounding. Second,
candidate ranking was stable under equal weights and under 0.5-fold or twofold
perturbations of individual evidence weights. Third, exact-residue matching
removed many gene-level associations that would otherwise imply unsupported
convergence between ligandability and cysteine function.

The benchmark does not establish external validation of the composite CPT Score.
Cysteine-derived evidence was available for only 7.5-7.9% of the retained
dependency associations, and computed ADME for 22.8%. In a CovPDB positive-set
analysis, the ligandability dimension enriched structurally observed
cysteine-covalent targets within its covered subset, whereas the full CPT Score
was approximately neutral. This result identifies evidence coverage, rather than
score arithmetic or weight instability, as the current limiting factor.

## Benchmark design

The analysis used the 14,877 gene-cohort RNAi dependency associations retained
by the prespecified manuscript criteria: 5,570 adult and 9,307 paediatric
associations across ten cohorts in each panel.

Four tests were performed:

1. **Implementation concordance.** Scores were independently reconstructed from
   the benchmark input matrices and compared with direct calls to
   `cpt_annotate_gene_table()` for every cohort.
2. **Weight sensitivity.** Nineteen configurations were evaluated against the
   shipped default, which weights ADME 0: ADME activated at the proposed 0.5,
   equal weights over six dimensions, equal weights over the five active ones,
   and leave-one-out or 0.5-fold/twofold perturbations of each active dimension.
3. **Positive-set enrichment.** Human cysteine-covalent targets were obtained
   from the CovPDB cysteine index and mapped through the UniProt REST API.
   Recovery was measured at the top 5%, 10%, and 20% within each cohort. CovPDB
   membership was treated as a positive set, not as a complete binary ground
   truth.
4. **Evidence-resolution benchmark.** Gene-level co-occurrence of ligandability
   and cysteine function was compared with exact-residue matching using the same
   workflow inputs and molecule selectivity rules.

Missing dimensions were excluded and the contributing weights renormalised, as
in the application. At least two finite, positively weighted dimensions were
required. Rank correlations were calculated only among candidates scoreable in
both configurations; top-set overlap was calculated over the full cohort
candidate universe so that loss of score coverage was not hidden.

## Main results

| Test | Result | Interpretation |
|---|---:|---|
| Implementation concordance | 14,877/14,877 displayed scores matched | The benchmark reproduces the application implementation at its displayed precision. |
| ADME activated at 0.5 | Median Spearman rho 0.999; minimum 0.998 | Activating the sixth dimension at the proposed weight moves the ranking very little. |
| ADME activated at 0.5, top 25 | Median overlap 0.96; minimum 0.92 | Shortlist movement is confined to the selection boundary. |
| Equal weights over six dimensions | Median Spearman rho 0.979; minimum 0.963 | Ranking is not dependent on the exact default weight ratios. |
| Equal weights over six dimensions, top 25 | Median overlap 0.88; minimum 0.68 | Equalising ADME with the evidence dimensions does move the shortlist. |
| Equal weights over five dimensions | Median Spearman rho 0.994; minimum 0.974; top-25 median 1.00, minimum 0.92 | Leading candidates were largely preserved when ADME stays inactive. |
| 0.5-fold/twofold single-weight changes | Median Spearman rho 0.988-1.000; minimum 0.966 | Reasonable one-dimension perturbations caused limited rank movement. |
| 0.5-fold/twofold single-weight changes, top 25 | Median overlap 0.92-1.00; minimum 0.80 | Some cohort-specific movement remains near the selection boundary. |
| Dependency and selectivity coverage | 100% each | These dimensions currently anchor the composite score. |
| Computed ADME coverage | 22.79% (3,391/14,877) | Available wherever a covalent probe engages the gene at CR >= 4. |
| Ligandability, conservation, clinical coverage | 7.85%, 7.49%, 7.85% | Most candidates have only the two dependency-derived dimensions. |
| Leave out dependency or selectivity | Median score coverage 8.05%; minimum 3.74% | The apparent high conditional correlation after omission should not be interpreted without the accompanying coverage collapse. |
| CovPDB positive set | 91 human targets; 86 in the RNAi background; 98 gene-cohort positives | The set provides independent structural positives but not cancer-specific negatives. |
| CovPDB top-10% enrichment, CPT default | 1.02-fold, cohort-bootstrap 95% CI 0.28-1.86 | No evidence of enrichment for the composite score. |
| CovPDB top-10% enrichment, ligandability only | 2.44-fold, 95% CI 1.58-3.39 | The ligandability dimension recovers independent structural positives within its covered subset. |
| Exact-site retention, <=20-target molecules | Adult 17/63 (27.0%); paediatric 23/61 (37.7%) | Gene-level evidence often joined different cysteines in the same gene. The discarded records are unsupported evidence upgrades, not experimentally proven biological false positives. |

## What changed against the 3 August 2026 run

Only the ADME-dependent results moved. Concordance, the five-dimension weight
scenarios, dimension coverage for the five evidence dimensions, and the
exact-site benchmark are identical to the digit. CovPDB fold enrichments are
identical; their bootstrap intervals shift in the third digit because the extra
ADME ranker changes how many draws precede each resample, not because the data
changed.

| Quantity | 3 Aug 2026 run | 19 Sep 2026 run | Manuscript |
|---|---:|---:|---:|
| Scores reproduced | 14,877/14,877 | 14,877/14,877 | 14,877 |
| ADME computed for | not computed | 22.79% | 22.8% |
| ADME at 0.5, median rho | not evaluated | 0.999 | 0.999 |
| ADME at 0.5, minimum rho | not evaluated | 0.998 | 0.998 |
| ADME at 0.5, median top-25 | not evaluated | 0.96 | 96% |
| ADME at 0.5, minimum top-25 | not evaluated | 0.92 | 92% |
| Equal weights, median rho | 0.994 (five dims) | 0.979 (six dims) | 0.979 |
| Equal weights, minimum rho | 0.974 (five dims) | 0.963 (six dims) | 0.963 |
| Equal weights, median top-25 | 1.00 (five dims) | 0.88 (six dims) | 88% |
| Equal weights, minimum top-25 | 0.92 (five dims) | 0.68 (six dims) | 68% |
| Ligandability / conservation / clinical coverage | 7.85 / 7.49 / 7.85% | 7.85 / 7.49 / 7.85% | 7.5-7.9% |
| CovPDB top-10%, CPT default | 1.02-fold | 1.02-fold | not quoted |
| CovPDB top-10%, ligandability only | 2.44-fold | 2.44-fold | not quoted |
| Exact-site retention, <=20 targets | 27.0% / 37.7% | 27.0% / 37.7% | 27.0% / 37.7% |

## Resource comparison

CanProTarget was compared qualitatively with Open Targets, CysDB, DrugMap,
CovPDB, and CovalentInDB 2.0. These resources answer related but non-identical
questions and do not emit directly comparable cohort-specific CanProTarget
rankings. A forced quantitative head-to-head would therefore be misleading. The
capability matrix in `outputs/software_resource_capability_comparison.csv`
records the comparison: CanProTarget uniquely combines cohort-specific cancer
dependency, exact-residue engagement, residue-function evidence, and
configurable integrated ranking, whereas the other resources provide
complementary disease, chemoproteomic, or structural evidence.

## Reproducibility

From the repository root, with the 2026-07-13 dependency snapshot available:

```bash
Rscript analysis/software_benchmark/work/parse_covpdb_reference.R
Rscript analysis/software_benchmark/work/run_software_benchmark.R \
  analysis/software_benchmark \
  /path/to/manuscript_revision_2026-07-13/data/benchmark
```

The second argument, or `CPT_BENCHMARK_INPUTS`, points at the directory holding
`rna_all_gene_statistics.rds`, `gene_level_vs_exact_summary.csv` and
`cohort_workflow_counts.csv`. These are outside the repository because they are
the frozen manuscript dependency universe rather than application data. The run
records input paths, MD5 checksums, the R version, external endpoints, and the
random seed in `outputs/benchmark_provenance.csv`. The benchmark used seed
`20260803` and 2,000 cohort-level bootstrap resamples. Publication should use a
frozen repository release and archive the external CovPDB and UniProt snapshots
with the output tables.

## Key files

- `outputs/cpt_implementation_validation.csv`: cohort-level score concordance.
- `outputs/cpt_weight_robustness_summary.csv`: scenario-level coverage, correlation, and top-set overlap.
- `outputs/cpt_weight_robustness_by_cohort.csv`: full cohort-by-scenario results.
- `outputs/cpt_dimension_coverage.csv`: per-dimension availability, ADME included.
- `outputs/covpdb_positive_set_enrichment_summary.csv`: fold enrichment with bootstrap intervals.
- `outputs/exact_site_discrimination_summary.csv`: gene-level versus exact-residue retention.
- `outputs_2026-08-03/`, `figures_2026-08-03/`: the superseded run, kept for comparison.
