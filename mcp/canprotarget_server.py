"""
CanProTarget MCP Server

Exposes cancer protein target prioritization tools via the Model Context Protocol.
Communicates with a warm R subprocess to execute statistical analyses on DepMap
dependency data and the Cys_editing functional cysteine atlas.

Usage:
    # Direct (stdio transport for IDE/agent integration):
    python mcp/canprotarget_server.py

    # Via MCP inspector:
    mcp dev mcp/canprotarget_server.py
"""

import logging
import os
import sys
from pathlib import Path
from typing import Annotated, Optional

# Prevent the local mcp/ directory from shadowing the installed 'mcp' package.
# This script lives inside mcp/, so when run from the project root, Python would
# resolve 'import mcp' to this directory instead of the pip-installed package.
_this_dir = str(Path(__file__).parent)
_project_root = str(Path(__file__).parent.parent)
sys.path = [p for p in sys.path if os.path.abspath(p) not in (_this_dir, _project_root, "", ".")]
# Re-add project root AFTER the mcp package path so installed packages take priority
import importlib
import mcp as _mcp_pkg  # noqa: ensure installed mcp is found
sys.path.insert(0, _this_dir)  # for r_bridge import

from mcp.server.fastmcp import FastMCP
from r_bridge import (
    RBridge,
    RBridgeDomainError,
    RBridgeError,
    RBridgeTransportError,
)

# --- Logging setup ---
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(name)s] %(levelname)s: %(message)s",
    stream=sys.stderr,
)
logger = logging.getLogger("canprotarget_mcp")

# --- Server setup ---
mcp = FastMCP(
    "CanProTarget",
    instructions=(
        "Cancer Protein Target prioritization platform. "
        "Query DepMap CRISPR/RNAi dependency data, identify selective "
        "cancer dependencies, and explore functional cysteine sites for drug targeting."
    ),
)

# --- R Bridge (lazy initialization) ---
_bridge: RBridge | None = None


def get_bridge() -> RBridge:
    """Get or create the R bridge connection."""
    global _bridge
    if _bridge is None or not _bridge.is_running:
        project_root = os.environ.get(
            "CPT_PROJECT_ROOT",
            str(Path(__file__).parent.parent),
        )
        _bridge = RBridge(project_root=project_root)
        _bridge.start()
    return _bridge


def r_query(tool: str, params: dict | None = None):
    """Query R worker. Restart only on transport failures, never on domain errors.

    Domain errors (unknown gene, bad subtype, invalid dataset) are returned
    immediately so agents do not thrash a multi-second data reload.
    """
    bridge = get_bridge()
    try:
        return bridge.query(tool, params)
    except RBridgeDomainError:
        # Expected application failure — keep warm worker alive
        raise
    except RBridgeTransportError as e:
        logger.warning("R transport failure, restarting worker once: %s", e)
        try:
            bridge.restart()
        except RBridgeError as restart_err:
            logger.error("R worker restart failed: %s", restart_err)
            raise e from restart_err
        return bridge.query(tool, params)
    except RBridgeError as e:
        # Unknown subclass: restart once (legacy safety net)
        logger.warning("R query failed (generic), attempting restart: %s", e)
        bridge.restart()
        return bridge.query(tool, params)


# === MCP Tools ===


@mcp.tool()
def list_subtypes(
    dataset: Annotated[str, "Dataset to query: 'CRISPR' or 'RNAi'"] = "CRISPR",
) -> dict:
    """
    List all available cancer subtypes (OncotreeSubtype) for a given dataset.

    Use this to discover valid subtype names before querying dependencies.
    Returns the full list of cancer subtypes that have enough cell lines
    for statistical analysis in the specified dataset.
    """
    return r_query("list_subtypes", {"dataset": dataset})


@mcp.tool()
def list_genes(
    dataset: Annotated[str, "Dataset to query: 'CRISPR' or 'RNAi'"] = "CRISPR",
) -> dict:
    """
    List all genes available in a dataset.

    Returns gene symbols from the DepMap gene effect matrix.
    CRISPR (23Q4) has ~18,443 genes; RNAi (DEMETER2) has ~17,309 genes.
    """
    return r_query("list_genes", {"dataset": dataset})


@mcp.tool()
def query_dependency(
    gene: Annotated[str, "Gene symbol (e.g. 'KRAS', 'EGFR', 'TP53')"],
    subtype: Annotated[str, "Cancer subtype from OncotreeSubtype (e.g. 'Non-Small Cell Lung Cancer')"],
    dataset: Annotated[str, "Dataset: 'CRISPR' or 'RNAi'"] = "CRISPR",
) -> dict:
    """
    Query the dependency score for a specific gene in a cancer subtype.

    Returns whether the gene is a dependency (mean effect < -0.5) and whether
    it is SELECTIVE (significantly more essential in this subtype vs others).

    Statistical approach:
    - Computes mean gene effect in the subtype vs all other cell lines
    - Runs a one-sided t-test (alternative: cancer is more essential)
    - Reports effect size (difference in means) and p-value

    A more negative effect size means stronger selective dependency.
    """
    return r_query("query_dependency", {
        "gene": gene,
        "subtype": subtype,
        "dataset": dataset,
    })


@mcp.tool()
def top_dependencies(
    subtype: Annotated[str, "Cancer subtype (OncotreeSubtype)"],
    dataset: Annotated[str, "Dataset: 'CRISPR' or 'RNAi'"] = "CRISPR",
    n: Annotated[int, "Number of top genes to return (max 100)"] = 20,
) -> dict:
    """
    Get the top N most selectively essential genes for a cancer subtype.

    Ranks genes by effect size (cancer subtype mean minus other cell lines mean).
    The most negative effect sizes indicate genes that are selectively essential
    in the specified cancer type compared to all other cancers.

    This is the primary tool for target discovery — use it to find which genes
    a specific cancer type depends on that other cancers do not.
    """
    return r_query("top_dependencies", {
        "subtype": subtype,
        "dataset": dataset,
        "n": n,
    })


@mcp.tool()
def gene_cysteines(
    gene: Annotated[str, "Gene symbol (e.g. 'EGFR', 'KRAS', 'PIK3CA')"],
) -> dict:
    """
    Get cysteine sites for a gene, merging the Cys_editing atlas with
    chemoproteomic SMCL engagement (CR >= 4).

    Returns, per residue:
    - evidence_tier 1-4 (exact functional+atlas ligandable; functional;
      tested-non-functional; untested / not in atlas)
    - atlas functional / ligandable flags when the site was assayed
    - engaging SMCL records (probe name, CR, n_targets) when present

    Data sources: Li et al., Nat Chem Biol 2023; bundled CysDB-indexed
    chemoproteomic competition-ratio table (not a live CysDB query).
    """
    return r_query("gene_cysteines", {"gene": gene})


@mcp.tool()
def cysteine_detail(
    site_id: Annotated[str, "Cysteine site ID in format 'GENE_POSITION' (e.g. 'EGFR_797')"],
) -> dict:
    """
    Get full annotation for a specific cysteine site (e.g. EGFR_797).

    Returns atlas fields when the site was assayed (ABE/CBE dropout,
    ligandability, conservation, ClinVar) plus evidence_tier and engaging
    SMCL records (probe, CR, n_targets). Sites known only from chemoproteomics
    are returned as Tier 4 (untested), not as an error.

    Use gene_cysteines first to discover available site IDs for a gene.
    """
    return r_query("cysteine_detail", {"site_id": site_id})


@mcp.tool()
def compare_subtypes(
    gene: Annotated[str, "Gene symbol"],
    subtype1: Annotated[str, "First cancer subtype (OncotreeSubtype)"],
    subtype2: Annotated[str, "Second cancer subtype (OncotreeSubtype)"],
    dataset: Annotated[str, "Dataset: 'CRISPR' or 'RNAi'"] = "CRISPR",
) -> dict:
    """
    Compare the dependency of a gene between two cancer subtypes.

    Performs a two-sided t-test comparing gene effect scores between
    cell lines of subtype1 vs subtype2. Use this to understand if a
    gene is more essential in one cancer type than another.

    Useful for:
    - Validating selectivity of a target
    - Comparing related cancer types (e.g., lung adeno vs squamous)
    - Understanding tissue-specific dependencies
    """
    return r_query("compare_subtypes", {
        "gene": gene,
        "subtype1": subtype1,
        "subtype2": subtype2,
        "dataset": dataset,
    })


@mcp.tool()
def platform_info() -> dict:
    """
    Get CanProTarget platform metadata, version, data sources, and citation info.

    Returns platform version, available datasets with their sizes,
    data provenance, citation information, and list of available tools.
    """
    return r_query("platform_info")


@mcp.tool()
def canprotarget_score(
    gene: Annotated[str, "Gene symbol (e.g. 'KRAS', 'EGFR')"],
    subtype: Annotated[str, "Cancer subtype (OncotreeSubtype)"],
    dataset: Annotated[str, "Dataset: 'CRISPR' or 'RNAi'"] = "CRISPR",
) -> dict:
    """
    Compute the CanProTarget Score for a gene in a specific cancer subtype.

    The CPT Score is a composite metric (0-100) combining multiple evidence axes:
    - Dependency strength: how essential is this gene in the cancer?
    - Cancer selectivity: is it specific to this cancer or broadly essential?
    - Cysteine ligandability: can we drug it with a covalent inhibitor?
    - Evolutionary conservation: is the targetable site conserved?
    - Clinical evidence: ClinVar pathogenic annotations at the site?
    - ADME/drug-likeness: fragment developability of engaging probes (weight 0 by default)

    Higher scores indicate better druggable cancer targets. Scores above 70
    are high-priority, 40-70 moderate, below 40 low-priority.

    Returns the composite score plus individual dimension breakdowns.
    """
    return r_query("canprotarget_score", {
        "gene": gene,
        "subtype": subtype,
        "dataset": dataset,
    })


@mcp.tool()
def pancancer_profile(
    gene: Annotated[str, "Gene symbol (e.g. 'KRAS', 'EGFR', 'TP53')"],
    dataset: Annotated[str, "Dataset: 'CRISPR' or 'RNAi'"] = "CRISPR",
    top_n: Annotated[int, "Number of top subtypes to return (default 30)"] = 30,
) -> dict:
    """
    Get the pan-cancer dependency profile for a gene across all cancer subtypes.

    Shows how essential a gene is in every cancer type simultaneously.
    Answers the key question: "Is this target specific to my cancer or pan-essential?"

    Returns subtypes ranked by dependency strength (most dependent first),
    with mean gene effect, standard deviation, and cell line counts.

    A gene that is a dependency in >10 subtypes is likely a common essential
    (broadly required for cell survival, harder to target selectively).
    A gene dependent in 1-3 subtypes is a selective dependency (ideal target).
    """
    return r_query("pancancer_profile", {
        "gene": gene,
        "dataset": dataset,
        "top_n": top_n,
    })


@mcp.tool()
def rank_targets(
    subtype: Annotated[str, "Cancer subtype (OncotreeSubtype)"],
    dataset: Annotated[str, "Dataset: 'CRISPR' or 'RNAi'"] = "CRISPR",
    n: Annotated[int, "Number of ranked targets to return (max 100)"] = 20,
    pool: Annotated[int, "Selective-dependency pool size to re-score before cutting to n"] = 100,
    require_ligandable: Annotated[bool, "If true, only return genes with cysteine ligandability data"] = False,
) -> dict:
    """
    Rank genes in a cancer subtype by gene-level CanProTarget (CPT) Score.

    Takes the most selective dependencies, scores them on available CPT
    dimensions, and returns the top n by composite score. Each gene includes
    best_engaged_site (residue, evidence tier, top SMCL) when chemoproteomic
    records exist. Gene-level order is unchanged by Site CPT.

    For residue-resolved ranking that never promotes a worse evidence tier,
    use rank_site_targets.
    """
    return r_query("rank_targets", {
        "subtype": subtype,
        "dataset": dataset,
        "n": n,
        "pool": pool,
        "require_ligandable": require_ligandable,
    })


@mcp.tool()
def rank_site_targets(
    subtype: Annotated[str, "Cancer subtype (OncotreeSubtype)"],
    dataset: Annotated[str, "Dataset: 'CRISPR' or 'RNAi'"] = "CRISPR",
    n: Annotated[int, "Number of ranked sites to return (max 100)"] = 20,
    pool: Annotated[int, "Gene pool to expand into sites"] = 100,
    min_cr: Annotated[float, "Minimum competition ratio to count as engagement"] = 4,
    max_targets: Annotated[int, "Maximum n_targets for a prioritised SMCL (paper default 20)"] = 20,
    effect_size_max: Annotated[Optional[float], "If set, keep genes with effect size below this (paper RNAi: -0.1). Omit to use the CPT dependency pool."] = None,
    p_max: Annotated[Optional[float], "If set, keep genes with one-sided Welch p below this (paper: 0.05)"] = None,
    exclude_common_essentials: Annotated[bool, "If true, drop genes with whole-matrix mean gene effect < -0.5"] = False,
) -> dict:
    """
    Rank engaged cysteines in a subtype by evidence tier, then Site CPT.

    Default gene pool is the CPT dependency ranking. To match the manuscript
    RNAi discovery filters, set dataset='RNAi', effect_size_max=-0.1, p_max=0.05,
    exclude_common_essentials=True, min_cr=4, max_targets=20.

    Site CPT uses the gene's dependency/selectivity percentiles plus cysteine
    dimensions from that exact residue. Results are ordered by evidence tier
    first (1, then 2, 3, 4), so a high-scoring untested site cannot outrank a
    functional Tier 1 or 2 site.
    """
    params = {
        "subtype": subtype,
        "dataset": dataset,
        "n": n,
        "pool": pool,
        "min_cr": min_cr,
        "max_targets": max_targets,
        "exclude_common_essentials": exclude_common_essentials,
    }
    if effect_size_max is not None:
        params["effect_size_max"] = effect_size_max
    if p_max is not None:
        params["p_max"] = p_max
    return r_query("rank_site_targets", params)


@mcp.tool()
def assess_target(
    gene: Annotated[str, "Gene symbol (e.g. 'KRAS', 'EGFR')"],
    subtype: Annotated[str, "Cancer subtype (OncotreeSubtype)"],
    dataset: Annotated[str, "Dataset: 'CRISPR' or 'RNAi'"] = "CRISPR",
) -> dict:
    """
    One-shot target assessment for agent workflows.

    Bundles dependency stats, gene-level CPT Score, residue-resolved engaged
    sites (SMCL, evidence tier, Site CPT), cysteine summary, pan-cancer
    breadth, suggested next tool calls, and citation. Use this when a user
    asks "how good is gene X in cancer Y?" instead of chaining tools.
    """
    return r_query("assess_target", {
        "gene": gene,
        "subtype": subtype,
        "dataset": dataset,
    })


@mcp.tool()
def generate_report(
    report_type: Annotated[str, "Report type: 'gene_dependency' or 'cysteine_target'"],
    gene: Annotated[str, "Gene symbol (e.g. 'KRAS', 'EGFR')"],
    subtype: Annotated[str, "Cancer subtype (required for gene_dependency reports)"] = "",
    dataset: Annotated[str, "Dataset: 'CRISPR' or 'RNAi'"] = "CRISPR",
    site_id: Annotated[str, "Specific cysteine site ID (optional, for cysteine_target reports)"] = "",
) -> dict:
    """
    Generate a branded CanProTarget report as a self-contained HTML file.

    Report types:
    - gene_dependency: Full gene analysis including dependency stats,
      CPT Score breakdown, pan-cancer profile, and druggable cysteines.
      Requires: gene, subtype.
    - cysteine_target: Cysteine-level annotation report with ABE/CBE dropout
      plot, ligandability, conservation, and ClinVar data.
      Requires: gene. Optionally: site_id for a specific cysteine.

    The generated report includes CanProTarget branding, data provenance,
    citation information, and links to the source repository. It is designed
    to be shared as supplementary material or included in publications.

    Returns the file path to the generated HTML report and metadata.
    Always include the citation and source link when sharing results.
    """
    params = {"gene": gene, "dataset": dataset, "report_type": report_type}
    if subtype:
        params["subtype"] = subtype
    if site_id:
        params["site_id"] = site_id
    return r_query("generate_report", params)


# === Entry point ===

if __name__ == "__main__":
    transport = os.environ.get("MCP_TRANSPORT", "stdio")
    logger.info("Starting CanProTarget MCP server (transport=%s)", transport)
    mcp.run(transport=transport)
