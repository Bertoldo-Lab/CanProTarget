"""
Tests for the CanProTarget MCP server.

Run from project root:
    mcp/.venv/bin/pytest tests/test_mcp_server.py -v

Tests the full stack: Python MCP server → R bridge → R worker → data.
"""

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

# Project root
PROJECT_ROOT = str(Path(__file__).parent.parent)
os.environ["CPT_PROJECT_ROOT"] = PROJECT_ROOT

# Add mcp/ to path for imports
sys.path.insert(0, str(Path(PROJECT_ROOT) / "mcp"))


# === Fixtures ===


@pytest.fixture(scope="session")
def bridge():
    """Shared R bridge instance for all tests (avoids repeated data loading)."""
    from r_bridge import RBridge

    b = RBridge(project_root=PROJECT_ROOT)
    b.start()
    yield b
    b.stop()


# === R Bridge Tests ===


class TestRBridge:
    """Tests for the R subprocess bridge."""

    def test_bridge_starts(self, bridge):
        assert bridge.is_running

    def test_ping(self, bridge):
        result = bridge.query("ping")
        assert result["status"] == "ok"
        assert "timestamp" in result

    def test_ping_method(self, bridge):
        assert bridge.ping() is True

    def test_invalid_tool(self, bridge):
        from r_bridge import RBridgeDomainError

        with pytest.raises(RBridgeDomainError, match="Unknown tool"):
            bridge.query("nonexistent_tool")

    def test_missing_params(self, bridge):
        from r_bridge import RBridgeDomainError

        with pytest.raises(RBridgeDomainError):
            bridge.query("query_dependency", {"gene": "KRAS"})  # missing subtype


# === API Function Tests (via bridge) ===


class TestListSubtypes:
    """Tests for the list_subtypes tool."""

    def test_crispr_subtypes(self, bridge):
        result = bridge.query("list_subtypes", {"dataset": "CRISPR"})
        assert result["dataset"] == "CRISPR"
        assert result["count"] > 50
        assert isinstance(result["subtypes"], list)
        assert result["subtypes"] == sorted(result["subtypes"])

    def test_rnai_subtypes(self, bridge):
        result = bridge.query("list_subtypes", {"dataset": "RNAi"})
        assert result["dataset"] == "RNAi"
        assert result["count"] > 30

    def test_invalid_dataset(self, bridge):
        from r_bridge import RBridgeDomainError

        with pytest.raises(RBridgeDomainError, match="must be"):
            bridge.query("list_subtypes", {"dataset": "INVALID"})

    def test_invalid_dataset_does_not_kill_worker(self, bridge):
        """Domain errors must not require a restart (transport stays warm)."""
        from r_bridge import RBridgeDomainError

        pid_before = bridge._process.pid if bridge._process else None
        with pytest.raises(RBridgeDomainError):
            bridge.query("list_subtypes", {"dataset": "INVALID"})
        assert bridge.is_running
        assert bridge._process.pid == pid_before
        # Still healthy after domain error
        assert bridge.query("ping")["status"] == "ok"


class TestListGenes:
    """Tests for the list_genes tool."""

    def test_crispr_genes(self, bridge):
        result = bridge.query("list_genes", {"dataset": "CRISPR"})
        assert result["count"] > 18000
        assert "KRAS" in result["genes"]
        assert "EGFR" in result["genes"]

    def test_rnai_genes(self, bridge):
        result = bridge.query("list_genes", {"dataset": "RNAi"})
        assert result["count"] > 17000


class TestQueryDependency:
    """Tests for the query_dependency tool."""

    def test_kras_pancreatic(self, bridge):
        result = bridge.query("query_dependency", {
            "gene": "KRAS",
            "subtype": "Pancreatic Adenocarcinoma",
            "dataset": "CRISPR",
        })
        assert result["gene"] == "KRAS"
        assert result["is_selective_dependency"] is True
        assert result["effect_size"] < -1.0
        assert result["p_value"] < 0.001

    def test_non_dependency(self, bridge):
        result = bridge.query("query_dependency", {
            "gene": "A1BG",
            "subtype": "Melanoma",
            "dataset": "CRISPR",
        })
        assert result["is_dependency"] is False

    def test_invalid_gene(self, bridge):
        from r_bridge import RBridgeDomainError

        with pytest.raises(RBridgeDomainError, match="not found"):
            bridge.query("query_dependency", {
                "gene": "NOTREALGENE",
                "subtype": "Melanoma",
                "dataset": "CRISPR",
            })

    def test_invalid_subtype(self, bridge):
        from r_bridge import RBridgeDomainError

        with pytest.raises(RBridgeDomainError, match="not found"):
            bridge.query("query_dependency", {
                "gene": "KRAS",
                "subtype": "Fake Cancer Type",
                "dataset": "CRISPR",
            })

    def test_case_insensitive_gene(self, bridge):
        result = bridge.query("query_dependency", {
            "gene": "kras",
            "subtype": "Pancreatic Adenocarcinoma",
            "dataset": "CRISPR",
        })
        assert result["gene"] == "KRAS"


class TestTopDependencies:
    """Tests for the top_dependencies tool."""

    def test_melanoma_top_5(self, bridge):
        result = bridge.query("top_dependencies", {
            "subtype": "Melanoma",
            "dataset": "CRISPR",
            "n": 5,
        })
        assert result["subtype"] == "Melanoma"
        assert len(result["genes"]) == 5
        # Should be sorted by effect_size (ascending = most negative first)
        effects = [g["effect_size"] for g in result["genes"]]
        assert effects == sorted(effects)

    def test_known_dependency_in_top(self, bridge):
        result = bridge.query("top_dependencies", {
            "subtype": "Melanoma",
            "dataset": "CRISPR",
            "n": 10,
        })
        genes = [g["gene"] for g in result["genes"]]
        # SOX10 and BRAF are well-known melanoma dependencies
        assert any(g in genes for g in ["SOX10", "BRAF"])

    def test_max_cap(self, bridge):
        result = bridge.query("top_dependencies", {
            "subtype": "Melanoma",
            "dataset": "CRISPR",
            "n": 200,  # should be capped at 100
        })
        assert len(result["genes"]) <= 100
        assert result["top_n"] == 100

    def test_negative_n_clamped(self, bridge):
        result = bridge.query("top_dependencies", {
            "subtype": "Melanoma",
            "dataset": "CRISPR",
            "n": -5,
        })
        assert result["top_n"] >= 1
        assert len(result["genes"]) >= 1
        assert "is_dependency" in result["genes"][0]
        assert result.get("ranking") == "effect_size_mean_difference"


class TestGeneCysteines:
    """Tests for the gene_cysteines tool."""

    def test_egfr(self, bridge):
        result = bridge.query("gene_cysteines", {"gene": "EGFR"})
        assert result["gene"] == "EGFR"
        assert result["total_sites"] > 30
        assert result["n_functional"] > 0
        assert result["n_functional_ligandable"] > 0

    def test_site_fields(self, bridge):
        result = bridge.query("gene_cysteines", {"gene": "EGFR"})
        site = result["sites"][0]
        assert "site_id" in site
        assert "functional" in site
        assert "ligandable" in site
        assert "cysteine_position" in site
        assert "evidence_tier" in site
        assert site["evidence_tier"] in (1, 2, 3, 4)
        assert "smcls" in site

    def test_case_insensitive(self, bridge):
        result = bridge.query("gene_cysteines", {"gene": "egfr"})
        assert result["gene"] == "EGFR"

    def test_missing_gene(self, bridge):
        from r_bridge import RBridgeError

        with pytest.raises(RBridgeError, match="not found"):
            bridge.query("gene_cysteines", {"gene": "ZZZZNOTREAL"})


class TestCysteineDetail:
    """Tests for the cysteine_detail tool."""

    def test_egfr_797(self, bridge):
        result = bridge.query("cysteine_detail", {"site_id": "EGFR_797"})
        assert result["gene_symbol"] == "EGFR"
        assert result["cysteine_position"] == 797
        assert result["functional"] is True
        assert result["ligandable"] is True
        assert result["ligandability_score"] == 100  # manually annotated

    def test_has_all_fields(self, bridge):
        result = bridge.query("cysteine_detail", {"site_id": "EGFR_797"})
        expected_fields = [
            "site_id", "gene_symbol", "cysteine_position",
            "functional", "ligandable", "editor_support",
            "abe_mean_lfc", "cbe_mean_lfc",
        ]
        for field in expected_fields:
            assert field in result, f"Missing field: {field}"

    def test_missing_site(self, bridge):
        from r_bridge import RBridgeError

        with pytest.raises(RBridgeError, match="not found"):
            bridge.query("cysteine_detail", {"site_id": "FAKE_999"})


class TestCompareSubtypes:
    """Tests for the compare_subtypes tool."""

    def test_braf_melanoma_vs_nsclc(self, bridge):
        result = bridge.query("compare_subtypes", {
            "gene": "BRAF",
            "subtype1": "Melanoma",
            "subtype2": "Non-Small Cell Lung Cancer",
            "dataset": "CRISPR",
        })
        assert result["gene"] == "BRAF"
        # BRAF should be more essential in melanoma (more negative)
        assert result["subtype1"]["mean_effect"] < result["subtype2"]["mean_effect"]
        assert result["p_value"] < 0.01

    def test_small_subtype_error(self, bridge):
        from r_bridge import RBridgeError

        with pytest.raises(RBridgeError, match="<3"):
            bridge.query("compare_subtypes", {
                "gene": "BRAF",
                "subtype1": "Melanoma",
                "subtype2": "Adenosquamous Carcinoma of the Pancreas",
                "dataset": "CRISPR",
            })


class TestPlatformInfo:
    """Tests for the platform_info tool."""

    def test_returns_metadata(self, bridge):
        result = bridge.query("platform_info")
        assert result["name"] == "CanProTarget"
        assert "version" in result
        assert len(result["available_tools"]) == 14
        assert "assess_target" in result["available_tools"]
        assert "rank_targets" in result["available_tools"]
        assert "rank_site_targets" in result["available_tools"]
        assert len(result["data_sources"]) == 4


class TestCanprotargetScore:
    """CPT score plumbing via the R bridge."""

    def test_kras_pdac_not_median_dummy(self, bridge):
        result = bridge.query("canprotarget_score", {
            "gene": "KRAS",
            "subtype": "Pancreatic Adenocarcinoma",
            "dataset": "CRISPR",
        })
        assert result["gene"] == "KRAS"
        assert result["cpt_score"] is not None
        dep = result["dimensions"]["dependency_strength"]["score"]
        sel = result["dimensions"]["cancer_selectivity"]["score"]
        # Single-gene scoring used to hardcode 50; background rank should not.
        assert abs(dep - 50) > 5
        assert abs(sel - 50) > 5
        assert "adme_druggability" in result["dimensions_missing"]
        assert result["n_dimensions_used"] == len(result["dimensions_active"])


class TestRankTargets:
    def test_pdac_rank(self, bridge):
        result = bridge.query("rank_targets", {
            "subtype": "Pancreatic Adenocarcinoma",
            "dataset": "CRISPR",
            "n": 10,
            "pool": 60,
        })
        assert result["n_returned"] > 0
        genes = [g["gene"] for g in result["genes"]]
        assert "KRAS" in genes
        scores = [g["cpt_score"] for g in result["genes"]]
        assert scores == sorted(scores, reverse=True)


class TestAssessTarget:
    def test_kras_bundle(self, bridge):
        result = bridge.query("assess_target", {
            "gene": "KRAS",
            "subtype": "Pancreatic Adenocarcinoma",
            "dataset": "CRISPR",
        })
        assert result["gene"] == "KRAS"
        assert "summary" in result
        assert result["dependency"]["is_selective_dependency"] is True
        assert len(result["next_actions"]) >= 2
        assert "citation" in result


# === MCP Protocol Tests ===


class TestMCPProtocol:
    """Tests for the full MCP JSON-RPC protocol."""

    @pytest.fixture(scope="class")
    def mcp_process(self):
        """Start MCP server process for protocol testing."""
        proc = subprocess.Popen(
            [sys.executable, "mcp/canprotarget_server.py"],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            cwd=PROJECT_ROOT,
        )
        yield proc
        proc.stdin.close()
        proc.wait(timeout=10)

    def _send_recv(self, proc, message):
        proc.stdin.write(json.dumps(message) + "\n")
        proc.stdin.flush()
        line = proc.stdout.readline()
        return json.loads(line)

    def test_initialize(self, mcp_process):
        resp = self._send_recv(mcp_process, {
            "jsonrpc": "2.0",
            "id": 1,
            "method": "initialize",
            "params": {
                "protocolVersion": "2024-11-05",
                "capabilities": {},
                "clientInfo": {"name": "test", "version": "1.0"},
            },
        })
        assert "result" in resp
        assert resp["result"]["serverInfo"]["name"] == "CanProTarget"
        assert "tools" in resp["result"]["capabilities"]

    def test_list_tools(self, mcp_process):
        # Send initialized notification first
        mcp_process.stdin.write(
            json.dumps({"jsonrpc": "2.0", "method": "notifications/initialized"}) + "\n"
        )
        mcp_process.stdin.flush()

        import time
        time.sleep(0.3)

        resp = self._send_recv(mcp_process, {
            "jsonrpc": "2.0",
            "id": 2,
            "method": "tools/list",
            "params": {},
        })
        tools = resp["result"]["tools"]
        tool_names = [t["name"] for t in tools]
        assert len(tools) == 14
        assert "query_dependency" in tool_names
        assert "gene_cysteines" in tool_names
        assert "rank_site_targets" in tool_names

    def test_call_tool(self, mcp_process):
        resp = self._send_recv(mcp_process, {
            "jsonrpc": "2.0",
            "id": 3,
            "method": "tools/call",
            "params": {
                "name": "platform_info",
                "arguments": {},
            },
        })
        content = resp["result"]["content"]
        assert len(content) > 0
        data = json.loads(content[0]["text"])
        assert data["name"] == "CanProTarget"
