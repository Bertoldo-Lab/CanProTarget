# Offline rebuild scripts

Run these from the repository root (`CanProTarget/`), not from this folder.

## Restore full DepMap 23Q4 CRISPR

Official Public 23Q4 CRISPR gene effect is 1,100 models × 18,443 genes (last
gene `ZZZ3`). Opening that CSV in Excel truncates it at `TNFRSF10C` (16,383
genes). Read it in R only.

Place these **gitignored** files under `data/`:

- `CRISPRGeneEffect_23Q4.csv` — DepMap Public 23Q4 gene effect
- `Model_v2.csv` — 23Q4 model metadata (do not use current Portal `Model.csv`)

```bash
Rscript docs/scripts/rebuild_depmap_23q4.R
```

Expect:

- `1100 models x 18443 genes; last=ZZZ3; TP53=TRUE; ACH-002509=FALSE`
- Acral Melanoma: 4 models in metadata, 3 in the CRISPR matrix

Writes `CRISPRGeneEffect_23Q4_clean.rds`, `cancer_model_data.rds`, and the
subtype lists. Those products are gitignored; do not commit the CSVs.

Optional local cache (also gitignored, ~1 hour):

```bash
Rscript docs/scripts/precompute_effectsizes_all_subtypes.R
```

`n_lines` in the log is metadata ∩ matrix. Acral CRISPR must be 3.

Do not remake `protein_binding_lookup_preprocessed.rds` for this rebuild.
Do not ingest 26Q1 CRISPR into the paper app.
