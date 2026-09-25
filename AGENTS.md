# CanProTarget - Agent Integration Guide

CanProTarget is a cancer protein target prioritization platform. It integrates DepMap gene dependency data with chemoproteomic ligandability and functional cysteine annotations to help researchers identify druggable cancer targets.

## How to use CanProTarget as an agent

CanProTarget is a domain-specific tool, not a general workbench. Use it when you need to answer questions about cancer target druggability, gene dependencies, or cysteine-targeted covalent chemistry. It does one thing well: prioritize cancer targets and tell you where the covalent handles are.

When integrating results into a broader workflow:
- Prefer `assess_target` / `rank_targets` for multi-step questions instead of chaining many calls
- Use structured JSON (including `provenance`, `methods_note`, `external_resources`) to chain into docking, literature, or clinical trial lookup elsewhere
- Generate reports via `generate_report` for self-contained, citable outputs
- Follow `external_resources` links (AlphaFold, UniProt, PDB, DepMap, OncoKB) rather than asking CanProTarget to render structures
- When using CanProTarget outputs in a paper or downstream tool, cite it — see [CITATION.cff](CITATION.cff).

## What this platform does

1. Identifies genes that are selectively essential in specific cancer types (using CRISPR/RNAi screens from DepMap)
2. Links dependency targets to covalent chemical probes (competition ratio data)
3. Annotates druggable cysteine residues using base-editing functional screens
4. Assesses drug-likeness of candidate probes (SwissADME properties)

## MCP Server

CanProTarget exposes an MCP server for programmatic access. The server uses a warm R subprocess that loads data once at startup, then responds to queries in milliseconds.

### Configuration

Add to your MCP client config:

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

Requires: Python 3.10+ with `mcp` package, R 4.x+ with `jsonlite`, data files in `data/`.

### Available tools

| Tool | Description | Key parameters |
|------|-------------|----------------|
| `list_subtypes` | Available cancer subtypes | `dataset`: "CRISPR" or "RNAi" |
| `list_genes` | All genes in a dataset | `dataset`: "CRISPR" or "RNAi" |
| `query_dependency` | Dependency score + stats for gene in a subtype | `gene`, `subtype`, `dataset` |
| `top_dependencies` | Ranked selective dependencies for a cancer type | `subtype`, `dataset`, `n` (default 20) |
| `gene_cysteines` | Functional/ligandable cysteine sites | `gene` |
| `cysteine_detail` | Full annotation for a specific cysteine | `site_id` (e.g. "EGFR_797") |
| `compare_subtypes` | Differential dependency between two cancers | `gene`, `subtype1`, `subtype2`, `dataset` |
| `canprotarget_score` | Composite CPT Score (0-100) with active/missing dimensions | `gene`, `subtype`, `dataset` |
| `rank_targets` | Top targets in a subtype by gene-level CPT Score | `subtype`, `dataset`, `n`, `require_ligandable` |
| `rank_site_targets` | Engaged cysteines ranked by evidence tier, then Site CPT | `subtype`, `dataset`, `n`, `max_targets`; manuscript filters: `effect_size_max`, `p_max`, `exclude_common_essentials` |
| `assess_target` | One-shot assessment (dependency + CPT + engaged sites/SMCL/tier + next steps) | `gene`, `subtype`, `dataset` |
| `pancancer_profile` | Dependency profile across all cancer types | `gene`, `dataset`, `top_n` (default 30) |
| `generate_report` | Branded HTML report | `report_type`, `gene`, `subtype` |
| `platform_info` | Version, data sources, citation | (none) |

### Example queries

```
"Is KRAS a dependency in pancreatic cancer?"
  -> query_dependency(gene="KRAS", subtype="Pancreatic Adenocarcinoma")

"What should I prioritize for covalent targeting in PDAC?"
  -> rank_targets(subtype="Pancreatic Adenocarcinoma", n=15)

"How good is KRAS as a target in pancreatic cancer?"
  -> assess_target(gene="KRAS", subtype="Pancreatic Adenocarcinoma")

"What are the top selective dependencies in lung cancer?"
  -> top_dependencies(subtype="Non-Small Cell Lung Cancer", n=10)

"What druggable cysteines does EGFR have?"
  -> gene_cysteines(gene="EGFR")

"Tell me about EGFR C797"
  -> cysteine_detail(site_id="EGFR_797")

"Is BRAF more essential in melanoma or colorectal cancer?"
  -> compare_subtypes(gene="BRAF", subtype1="Melanoma", subtype2="Colon Adenocarcinoma")

"What cancers can I query with CRISPR data?"
  -> list_subtypes(dataset="CRISPR")

"Show me where TP53 is essential across all cancers"
  -> pancancer_profile(gene="TP53")
```

### Understanding results

- **Gene effect scores**: Negative values indicate dependency. Below -0.5 is generally considered essential. Below -1.0 is strongly essential.
- **Effect size**: Difference between the cancer subtype mean and all other cell lines. More negative = more selectively essential in that cancer.
- **p-value**: From a one-sided t-test (cancer vs others). Significant selective dependencies have p < 0.05 and effect size < -0.2.
- **CPT Score dimensions**: Responses list which dimensions were active vs missing. Scores are not comparable across different numbers of active dimensions without checking `n_dimensions_used`.
- **`adme_druggability` is computed but weighted 0 by default.** It is scored for ~4,470 genes (mean fragment developability of covalent probes engaging the gene at CR >= 4) and reported in `dimensions`, but contributes nothing to `cpt_score` unless the weight is raised. The formula is provisional and pending methods review, so cite it as an indicator, not as a validated axis. `NA` means no probe engages the gene in the bundled screens, which is missing evidence and not a negative result.
- **Subtype names must match OncotreeSubtype exactly, and need >= 3 cell lines.** Near-miss names fail: "Colon Adenocarcinoma" has 45 CRISPR lines, but "Colorectal Adenocarcinoma" has only 1 and will error. Call `list_subtypes` first rather than guessing.
- **CPT Score methods**: Full description of weights, percentile ranking, composite formula, priority guards, and `rank_targets` pool logic is in [docs/CPT_SCORE.md](docs/CPT_SCORE.md). Implementation: `R/canprotarget_score.R`.
- **Priority labels**: Heuristic bands for ranking only. Non-dependencies (mean effect >= -0.5) are never labeled high priority.
- **Ligandability score**: 0-100 scale from the cysteine editing atlas. Higher = more evidence of ligandability from chemoproteomic competition.
- **Functional cysteines**: Sites where base-editing (ABE/CBE) causes fitness dropout, indicating the cysteine is important for protein function.
- **Absence is not evidence of absence.** Both the cysteine atlas and the probe table are experimental coverage maps. A gene or site that is missing was not detected in those screens; it has not been shown to be non-ligandable. Say "no data" rather than "not druggable".
- **Mutation-created cysteines are invisible to this platform.** The atlas and the probe screens profile the *wild-type* cysteinome. KRAS G12C is the clearest case: wild-type KRAS has glycine at codon 12, so no C12 site can exist in the atlas. CanProTarget lists only KRAS C80 and C118 (C118 functional, neither ligandable), and KRAS is absent from the probe table entirely while HRAS and NRAS are present. Do not read that as KRAS being undruggable. For mutant-specific covalent handles, follow `external_resources` out to UniProt/PDB and the medicinal-chemistry literature (sotorasib, adagrasib) instead.

## Data available

| Dataset | Coverage | Source |
|---------|----------|--------|
| CRISPR Gene Effect | 18,443 genes x 1,100 cell lines | DepMap 23Q4 |
| RNAi Gene Effect | 17,309 genes x 712 cell lines | DEMETER2 v6 |
| Cancer subtypes | 83 (CRISPR), 56 (RNAi) | OncotreeSubtype classification |
| Cysteine Editing Atlas | 13,872 sites, 1,778 genes | Li et al., Nat Chem Biol 2023 |
| Chemoproteomics (CR) | 998 probes x 8,425 genes x 37,774 cysteine sites | Six probe studies, 7 cell lines (see README) |
| SwissADME | 1,000 probes x 53 descriptors | swissadme.ch |

Wild-type coverage only, and per cysteine rather than per gene. 1,551 genes appear
in both the atlas and the probe table; 8,280 of the probe-table genes are also in
the CRISPR matrix.

## Project structure

```
R/api_functions.R          - Pure query functions (no Shiny deps)
R/mcp_worker.R             - Warm R subprocess for MCP server
R/functions.R              - Analysis functions (gene effect, probe linking)
R/canprotarget_score.R     - Composite scoring + annotate/radar/sensitivity
R/report_generator.R       - HTML reports
R/theme_canprotarget.R     - Publication ggplot2 theme
mcp/canprotarget_server.py - Python MCP server (FastMCP)
mcp/r_bridge.py            - R subprocess manager
app.R                      - Shiny app entry point
data/                      - Runtime RDS + data_versions.yaml (DepMap matrices gitignored)
data/raw/                  - Chemoproteomics rebuild inputs (xlsx via Git LFS)
tests/fixtures/            - Known-target recovery fixtures (YAML)
MISSING_DATA.md            - Data provenance, coverage limits, rebuild commands
```

## Performance note

There is no committed cache of subtype statistics, so `top_dependencies`,
`rank_targets` and the Shiny Dependencies tab compute limma plus t-tests over
~18,400 genes on demand. Budget **20-30 seconds per subtype** and do not treat a
slow first call as a failure. A local cache can be built with
`docs/scripts/precompute_effectsizes_all_subtypes.R`.

## Citation

See [CITATION.cff](CITATION.cff).