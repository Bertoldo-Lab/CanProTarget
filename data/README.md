# Bundled data

Runtime inputs and rebuild sources committed with the repo. Large DepMap
gene-effect matrices remain local / separately downloaded — see
`docs/DATA_PROVENANCE.md`.

Rebuild CRISPR 23Q4 from official CSV + `Model_v2.csv` (never Excel, never
current Portal `Model.csv`):

```bash
Rscript docs/scripts/rebuild_depmap_23q4.R
```

## Runtime (app reads these)

### `gene_index.rds`

Dependency summaries, engaged cysteines and probe rows for Discover and
Target. Built by `docs/scripts/build_gene_index.R` from the per-subtype
effect-size cache, the binding table, SwissADME and the atlas. The committed
file is what v1.0.0 runs. Rebuild it only if those inputs change.

### `cys_editing_atlas.rds`

the cysteine function layer atlas from
[`cravattlab/Cys_editing`](https://github.com/cravattlab/Cys_editing) at commit
`89bc6a268ea4385a687f2e193046b4486edd0185`, including conservation, ortholog,
and ClinVar fields. Rebuild:

```bash
Rscript docs/scripts/import_cys_editing_atlas.R /path/to/Cys_editing
```

Cite Li et al., *Nature Chemical Biology* (2023),
<https://doi.org/10.1038/s41589-023-01428-w>. License:
`docs/licenses/Cys_editing_LICENSE.txt`.

### `protein_binding_lookup_preprocessed.rds` / `swissadme_preprocessed.rds`

Chemistry and Chemistry tabs. Rebuild from `data/raw/` via the preprocess
scripts under `docs/scripts/`.

## Not tracked: `precomputed_effectsizes/`

Optional local cache of Discover effect-size tables (one TSV per
subtype × dataset). Gitignored on purpose: derived data, ~1 GB, and a truncated
copy fails silently rather than erroring. Build your own:

```bash
Rscript docs/scripts/precompute_effectsizes_all_subtypes.R
```

Without the cache the app recomputes each subtype in ~20-30 s. A cache you
did not build yourself should end at the last gene of the matrix (`ZZZ3` for
CRISPR 23Q4), not mid-alphabet.

## Rebuild sources

### `raw/`

Backus / SwissADME inputs (`*.xlsx` via Git LFS). See `raw/README.md`.
