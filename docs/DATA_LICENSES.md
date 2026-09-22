# Third-party data licenses

CanProTarget source code is GNU Affero GPL v3.0 (see `LICENSE`).
The bundled and required datasets keep their own terms. This file is the
record the manuscript cites as `DATA_LICENSES.md`. It sits under `docs/`
because GitHub's licence detector reads every root-level file whose name
contains "license" and renders it as a licence for the repository: at the root
this file appeared as a second License tab beside the AGPL, which it is not.

| Dataset | File | License / terms | Notes |
|---------|------|-----------------|-------|
| DepMap CRISPR 23Q4 gene effect | `data/CRISPRGeneEffect_23Q4_clean.rds` (downloaded; not git-tracked) | [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) | Broad Institute DepMap. Cite DepMap and the 23Q4 release. |
| DepMap RNAi DEMETER2 v6 | `data/d2_gene_effect_headers_refined.rds` (downloaded; not git-tracked) | [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) | Same portal. |
| DepMap model metadata | `data/cancer_model_data.rds` | CC BY 4.0 | Derived from DepMap `Model.csv`. |
| Cys_editing functional-cysteine atlas | `data/cys_editing_atlas.rds` | MIT (copyright 2023 Jason Li) | Li et al., *Nat Chem Biol* 2023. Full text: `docs/licenses/Cys_editing_LICENSE.txt`. Source: [cravattlab/Cys_editing](https://github.com/cravattlab/Cys_editing). |
| Chemoproteomic competition-ratio table | `data/protein_binding_lookup_preprocessed.rds` | Follows source publications | Compiled from **six** CysDB-indexed ligandability studies, not a live CysDB dump: Kuljanin/Gygi, Backus/Cravatt, Vinogradova/Cravatt, Yang/Wang, Cao/Backus, Yan/Backus (see `docs/DATA_PROVENANCE.md` for probe/gene/cell-line counts). Where redistribution of a source table is restricted, rebuild from `data/raw/` with the scripts in `docs/scripts/`. |
| SwissADME descriptors | `data/swissadme_preprocessed.rds` | [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) | [swissadme.ch](http://www.swissadme.ch/). Daina et al., *Sci Rep* 2017. |
| CovPDB covalent-ligand structures | `analysis/software_benchmark/data_external/covpdb_cysteine_*` (benchmark only; not bundled in the app) | Terms of the source database — see [CovPDB](https://drug-discovery.vm.uni-freiburg.de/covpdb/) | Used solely as an independent positive set in the CPT Score benchmark. The retrieved index pages and the tables parsed from them are retained so the benchmark reproduces without re-querying the site. |
| UniProt ID mapping | `analysis/software_benchmark/data_external/uniprot_covpdb_gene_mapping.tsv` (benchmark only) | [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) | Retrieved through the UniProt REST API to map CovPDB accessions to gene symbols. |

Absence of a gene or cysteine from the atlas or probe table is missing experimental coverage, not evidence that the site is non-ligandable or non-functional.

Datasets marked *benchmark only* are inputs to the CPT Score benchmark in
`analysis/software_benchmark/` and are not loaded by the application.

Rebuild commands and coverage limits: `docs/DATA_PROVENANCE.md`.
