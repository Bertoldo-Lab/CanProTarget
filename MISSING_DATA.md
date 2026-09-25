# Data notes

Discover, Target and Chemistry all run from data bundled in this repo. The
only inputs you still fetch yourself are the DepMap gene-effect matrices.

Data status is also listed in `data/data_versions.yaml` and surfaced in the
app's **About** tab.

## Bundled (in repo)

| Path | Tab / use | Notes |
|------|-----------|-------|
| `data/cys_editing_atlas.rds` | Discover: residue evidence; Target | Includes conservation, ortholog, ClinVar |
| `data/protein_binding_lookup_preprocessed.rds` | Discover: ligandability; Chemistry | Runtime; 10,588,541 rows: 998 probes, 8,425 genes, 37,774 cysteine sites, pooled from six studies across seven cell lines |
| `data/swissadme_preprocessed.rds` | Chemistry | Runtime; 1,000 probes x 53 descriptors; 996 join to the probe table |
| `data/raw/*` | Rebuild only | Reformatted chemoproteomics + SwissADME sources |

`data/raw/*.xlsx` is stored with Git LFS. After clone:

```bash
git lfs install
git lfs pull
```

Everything else in the repo is ordinary git, so a clone without Git LFS still
gives you a fully working app. Only the rebuild inputs will be placeholder
pointer files.

## Still local / download separately

DepMap inputs for the dependency layer:

- `CRISPRGeneEffect_23Q4_clean.rds` (or CSV → preprocess)
- `d2_gene_effect_headers_refined.rds`
- `cancer_model_data.rds` (from 23Q4 `Model_v2.csv`; do not use retracted Portal `Model.csv`)
- `cancer_subtypes_CRISPR.txt` / `cancer_subtypes_RNAi.txt`

## Precomputed subtype stats: not tracked, build locally

`data/precomputed_effectsizes/` is gitignored.
Without them, every subtype analysis is computed on demand, roughly **20-30
seconds per subtype**. The app is working correctly when it does this.

```bash
Rscript docs/scripts/precompute_effectsizes_all_subtypes.R
```

Writes one TSV per feasible subtype into `data/precomputed_effectsizes/`, which
the app prefers over recomputation. Budget about an hour for the full run.

**Do not commit these files, and do not accept them from another machine
without checking them first.** They are derived data, they are large (~5 MB per
subtype, ~1 GB total), and a partial file fails silently: the app reads it
without error and simply never sees the missing genes. Verify any candidate
cache with:

```r
# every file should end at the last gene of the matrix, not mid-alphabet
mat <- readRDS("data/CRISPRGeneEffect_23Q4_clean.rds")
n_expected <- ncol(mat)                       # 18,443 for CRISPR 23Q4
es <- readRDS("data/precomputed_effectsizes/CRISPR_Melanoma.rds")
nrow(es)                                      # should be close to n_expected
tail(es$gene_name, 1)                         # should be a Z gene, e.g. ZZZ3
```

Run that check after building a cache. A gene-effect matrix saved through
Excel is silently cut at its 16,384-column limit, and every statistic derived
from it is correct for the genes it kept, so nothing downstream complains.

## CRISPR 23Q4: rebuild in R, never Excel

Official DepMap Public 23Q4 `CRISPRGeneEffect.csv` is **1,100 x 18,443**, last
gene **`ZZZ3`**, with TP53 present. Excel's column limit is 16,384, so a copy
saved through Excel ends at **`TNFRSF10C`** and drops 2,508 genes including
TP53, WEE1, TOP2A, TYMS, VCP, WRN, XPO1 and YAP1.

Rebuild in R, never Excel:

```bash
# gitignored inputs under data/
#   CRISPRGeneEffect_23Q4.csv
#   Model_v2.csv   (23Q4 freeze; Portal Model.csv is retracted)
Rscript docs/scripts/rebuild_depmap_23q4.R
```

Expect `1100 models x 18443 genes; last=ZZZ3; TP53=TRUE; ACH-002509=FALSE` and
Acral Melanoma CRISPR n = 3.

`ACH-002509` (WM4235) is listed in 23Q4 `Model_v2.csv` (4 Acral models) but is
**not** in 23Q4 or 26Q1 CRISPR, so Acral Melanoma CRISPR is n = 3. RNAi
DEMETER2 has **zero** of these QIMR Acral lines.

Do **not** ingest current Portal `Model.csv` (2,154 models, 6 Acral) or 26Q1
CRISPR into the paper app. Pair 23Q4 CRISPR with 23Q4 `Model_v2.csv`.

## Coverage limits worth knowing

These are not gaps to be filled. They are properties of the assays the platform
is built on, and they should be stated plainly in the paper rather than
discovered by a reviewer.

### Absence means "not detected", never "not druggable"

Both the cysteine atlas (1,778 genes) and the probe table (8,425 genes) are
experimental coverage maps. 1,551 genes appear in both. A gene missing from
either was not observed in those screens; nothing was demonstrated about it.

### Mutation-created cysteines cannot appear at all

Both data types profile the **wild-type** cysteinome. A cysteine that exists only
because of a somatic mutation is structurally outside their scope.

KRAS is the clearest example, and it is worth understanding because KRAS/PDAC is
a headline case study:

| KRAS data | Status in CanProTarget |
|-----------|------------------------|
| Dependency (CRISPR + RNAi) | Present. `KRAS (3845)`; PDAC mean effect -2.03, selective, CPT 99.6 |
| Functional cysteines | Present. Two sites: **C118** functional but not ligandable (score 2.25), **C80** neither. Both fully conserved |
| Probe competition ratios | **Absent.** UniProt `P01116` is not in any of the six screens |
| Codon 12 (the G12C handle) | **Cannot exist here.** Wild-type KRAS has glycine at 12 |

The absence from the probe table is genuine non-detection, not a mapping bug: all
8,512 protein IDs in the table resolve through `id_mapping.tsv`, and the
homologues HRAS (`P01112`, 181 rows) and NRAS (`P01111`, 54 rows) are both
present. KRAS was simply never captured.

So CanProTarget will correctly report that KRAS has no probe evidence and no
ligandable cysteine, while KRAS G12C is a licensed drug target. That is the
system behaving as designed within its evidence base, not a defect.

### If you want broader cysteine coverage

**CysDB** is the natural superset. It aggregates nine chemoproteomic studies from
the Backus, Cravatt, Gygi, Wang and Yang groups into ~62,888 cysteines (about 24%
of the human cysteinome) with ligandability, hyperreactivity and disease
annotations. Six of those nine studies are already pooled in our probe table, so
CysDB would mainly add depth rather than a new kind of evidence.

- App: <https://backuslab.shinyapps.io/cysdb/>
- Boatner et al., *Cell Chemical Biology* (2023),
  <https://doi.org/10.1016/j.chembiol.2023.04.004>

CysDB is also wild-type, so it will **not** give you KRAS C12. Mutant-specific
covalent handles come from targeted medicinal chemistry and structural biology
(sotorasib, adagrasib) rather than cysteine chemoproteomics, which means ChEMBL
and the PDB, not a ligandability atlas. Adding that would be a different data
type and a scope decision, not a data refresh.

*Content in this section was rephrased from the cited sources for compliance with
licensing restrictions.*

## Raw chemoproteomics (`data/raw/`)

### table-s2.xlsx (Git LFS)

Not the raw Supplementary Table 1 from Backus *Nature* 2016. Reformatted with:

- "Ligandable Dataset" (wide CR matrix)
- "Compound Keys" (Compound_Name, Dataset, Cell_Line, SMILES)

### id_mapping.tsv

`From` (protein ID) → `To` (gene symbol).

### swissadme.csv

Batch export from swissadme.ch.

## Rebuild commands

```bash
# Chemoproteomics RDS from data/raw/ (see docs/DATA_PROVENANCE.md)
Rscript docs/scripts/preprocess_protein_binding.R

# DepMap 23Q4 CRISPR + Model_v2 (never Excel; do not use Portal Model.csv)
Rscript docs/scripts/rebuild_depmap_23q4.R

# Cys atlas (optional refresh from upstream Cys_editing)
git clone https://github.com/cravattlab/Cys_editing.git ../Cys_editing_reference
Rscript --vanilla docs/scripts/import_cys_editing_atlas.R ../Cys_editing_reference

# Precomputed subtype effect sizes (needs DepMap RDS + cancer_model_data.rds)
Rscript docs/scripts/precompute_effectsizes_all_subtypes.R
```

## Open question

The ADME dimension of the CPT Score is implemented but provisional, and carries
a default weight of 0 so that published scores are unchanged. Raising it above 0
changes every score, ranking and report, so the formula and its weight are a
methods decision rather than a setting. See `docs/CPT_SCORE.md` section 4.6.
