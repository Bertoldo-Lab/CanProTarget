# CanProTarget MCP Server

Model Context Protocol server for CanProTarget. Enables AI agents to query cancer dependency data, identify selective gene dependencies, and explore druggable cysteine sites.

## Architecture

```
AI Agent ↔ MCP Protocol (stdio)
                ↓
    Python MCP Server (FastMCP)
                ↓ JSON request via stdin
    R Subprocess (warm, persistent)
                ↓ loads data once at startup
    JSON response via stdout
```

The R subprocess loads large gene effect matrices (~18K genes × 1100 cell lines) once at startup, then handles queries in milliseconds. This avoids the 2-3 second cold start that would occur if R were launched per-query.

**Error handling:** Application errors (unknown gene, bad subtype, invalid dataset) are returned immediately without restarting R. Only transport failures (crash, timeout, broken pipe) restart the warm worker. On timeout the worker is killed before the next query so stdin/stdout cannot desync.

**Timeouts (optional env):**

| Variable | Default | Meaning |
|----------|---------|---------|
| `CPT_STARTUP_TIMEOUT` | 60 | Seconds to wait for R ready signal |
| `CPT_QUERY_TIMEOUT` | 120 | Seconds per tool call |

## Prerequisites

- **Python 3.10+** with `mcp` package
- **R 4.x+** with `jsonlite` package
- **Data files** in `data/` directory (see main project README)

## Setup

```bash
# From project root
cd /path/to/CanProTarget

# Install Python dependencies
pip install -r mcp/requirements.txt

# Verify R has jsonlite
Rscript --vanilla -e 'library(jsonlite); cat("OK\n")'
```

## Running

### Direct (stdio transport, for IDE/agent integration)

```bash
python mcp/canprotarget_server.py
```

### With MCP Inspector (for development/testing)

```bash
mcp dev mcp/canprotarget_server.py
```

### Configure in your MCP client

Add to your MCP configuration (e.g., `~/.kiro/settings/mcp.json`):

```json
{
  "mcpServers": {
    "canprotarget": {
      "command": "python",
      "args": ["mcp/canprotarget_server.py"],
      "cwd": "/path/to/CanProTarget",
      "env": {
        "CPT_PROJECT_ROOT": "/path/to/CanProTarget"
      }
    }
  }
}
```

## Available tools (14)

| Tool | Description |
|------|-------------|
| `list_subtypes` | Cancer subtypes for CRISPR or RNAi |
| `list_genes` | Genes in a dataset |
| `query_dependency` | Dependency score for a gene in a subtype |
| `top_dependencies` | Top selective dependencies for a cancer type |
| `gene_cysteines` | Atlas + engaged SMCL cysteines, with evidence tiers |
| `cysteine_detail` | Full annotation for a site (e.g. EGFR_797), including SMCL |
| `compare_subtypes` | Gene effect in subtype A vs B |
| `canprotarget_score` | CPT Score (0-100) with active/missing dimensions |
| `rank_targets` | Rank subtype genes by gene-level CPT Score |
| `rank_site_targets` | Rank engaged residues by evidence tier, then Site CPT |
| `assess_target` | One-shot dependency + CPT + engaged sites + next steps |
| `pancancer_profile` | Dependency across subtypes |
| `generate_report` | Branded HTML report |
| `platform_info` | Version, citation, credits |

Prefer `assess_target` / `rank_targets` for multi-step agent questions. Full examples: [AGENTS.md](../AGENTS.md).

## Example queries

```
# Top prioritization list
→ rank_targets(subtype="Pancreatic Adenocarcinoma", n=15)

# One-shot gene assessment
→ assess_target(gene="KRAS", subtype="Pancreatic Adenocarcinoma")

# Selective dependency stats
→ query_dependency(gene="KRAS", subtype="Pancreatic Adenocarcinoma")

# Cysteine handle
→ gene_cysteines(gene="EGFR")
→ cysteine_detail(site_id="EGFR_797")

# Lineage comparison (use Oncotree labels with enough cell lines)
→ compare_subtypes(gene="BRAF", subtype1="Melanoma", subtype2="Colon Adenocarcinoma")
```

## Data sources

- **DepMap CRISPR** (23Q4) and **RNAi** (DEMETER2 v6)
- **Cysteine Editing Atlas** (Li et al., Nat Chem Biol 2023)
- Chemoproteomics / SwissADME: not required for MCP dependency/cysteine tools; see [data/README.md](../data/README.md)

## Troubleshooting

**R worker fails to start:**
- Check that `Rscript` is on PATH: `which Rscript`
- Verify jsonlite: `Rscript --vanilla -e 'library(jsonlite)'`
- Check data files exist: `ls data/*.rds`

**Timeout on first query:**
- First query after startup loads large matrices (~3-5s). Subsequent queries are fast.
- Increase timeout via `CPT_QUERY_TIMEOUT` environment variable.

**Missing probe / SwissADME data:**
- Dependency, CPT, cysteine, and report tools work without those files.
- the Chemistry tab in the Shiny app needs the RDS files described in [data/README.md](../data/README.md).
