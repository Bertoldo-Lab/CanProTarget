"""
R subprocess bridge for CanProTarget MCP server.

Manages a warm R process that loads data once at startup and handles
JSON-line requests via stdin/stdout. This avoids the 2-3s cold start
penalty of launching Rscript for each query.
"""

import json
import logging
import os
import subprocess
import threading
import time
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)


class RBridgeError(Exception):
    """Base error from the R subprocess bridge."""
    pass


class RBridgeDomainError(RBridgeError):
    """Application-level R error (unknown gene/subtype, bad args).

    These are expected query failures. Callers must NOT restart the worker.
    """
    pass


class RBridgeTransportError(RBridgeError):
    """Process crash, timeout, broken pipe, or protocol desync.

    The worker may be dead or stdin/stdout may be out of sync; restart is safe.
    """
    pass


class RBridge:
    """Manages a persistent R subprocess for MCP tool execution."""

    def __init__(
        self,
        project_root: str | None = None,
        r_executable: str = "Rscript",
        startup_timeout: float | None = None,
        query_timeout: float | None = None,
    ):
        self.project_root = Path(project_root or self._find_project_root())
        self.r_executable = r_executable
        # Env overrides (documented in mcp/README.md)
        self.startup_timeout = float(
            os.environ.get("CPT_STARTUP_TIMEOUT", startup_timeout if startup_timeout is not None else 60.0)
        )
        self.query_timeout = float(
            os.environ.get("CPT_QUERY_TIMEOUT", query_timeout if query_timeout is not None else 120.0)
        )
        self._process: subprocess.Popen | None = None
        self._lock = threading.Lock()
        self._stderr_thread: threading.Thread | None = None
        self._stderr_lines: list[str] = []
        self._reader_alive = False  # True if a timed-out reader may still hold stdout

    @staticmethod
    def _find_project_root() -> str:
        """Find project root by looking for R/mcp_worker.R relative to this file."""
        # This file is at mcp/r_bridge.py, project root is one level up
        return str(Path(__file__).parent.parent)

    @property
    def is_running(self) -> bool:
        """Check if the R process is alive."""
        return self._process is not None and self._process.poll() is None

    def start(self) -> None:
        """Start the R worker process and wait for the ready signal."""
        if self.is_running:
            logger.info("R worker already running (pid=%d)", self._process.pid)
            return

        worker_script = self.project_root / "R" / "mcp_worker.R"
        if not worker_script.exists():
            raise RBridgeTransportError(f"R worker script not found: {worker_script}")

        env = os.environ.copy()
        env["CPT_PROJECT_ROOT"] = str(self.project_root)

        logger.info("Starting R worker: %s --vanilla %s", self.r_executable, worker_script)

        self._process = subprocess.Popen(
            [self.r_executable, "--vanilla", str(worker_script)],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            cwd=str(self.project_root),
            env=env,
            bufsize=1,  # Line-buffered
        )
        self._reader_alive = False

        # Start stderr reader thread (for logging, not protocol)
        self._stderr_lines = []
        self._stderr_thread = threading.Thread(
            target=self._read_stderr, daemon=True
        )
        self._stderr_thread.start()

        # Wait for ready signal
        ready = self._read_response(timeout=self.startup_timeout)
        if ready is None:
            self.stop()
            stderr_output = "\n".join(self._stderr_lines[-20:])
            raise RBridgeTransportError(
                f"R worker did not send ready signal within {self.startup_timeout}s.\n"
                f"Stderr:\n{stderr_output}"
            )

        if ready.get("status") != "ready":
            self.stop()
            raise RBridgeTransportError(f"Unexpected startup response: {ready}")

        logger.info(
            "R worker ready (pid=%d)", ready.get("pid", self._process.pid)
        )

    def stop(self) -> None:
        """Stop the R worker process."""
        if self._process is None:
            return

        logger.info("Stopping R worker (pid=%d)", self._process.pid)

        try:
            if self._process.stdin and not self._process.stdin.closed:
                self._process.stdin.close()
        except Exception:
            pass

        try:
            self._process.wait(timeout=5.0)
        except subprocess.TimeoutExpired:
            logger.warning("R worker did not exit gracefully, terminating")
            self._process.terminate()
            try:
                self._process.wait(timeout=3.0)
            except subprocess.TimeoutExpired:
                logger.warning("R worker did not terminate, killing")
                self._process.kill()

        self._process = None

    def restart(self) -> None:
        """Restart the R worker (e.g., after a crash)."""
        logger.info("Restarting R worker...")
        self.stop()
        time.sleep(0.5)
        self.start()

    def query(self, tool: str, params: dict[str, Any] | None = None) -> Any:
        """
        Send a query to the R worker and return the result.

        Raises:
            RBridgeDomainError: gene/subtype/arg errors from R (do not restart)
            RBridgeTransportError: crash, timeout, pipe, bad protocol (restart ok)
        """
        with self._lock:
            if not self.is_running:
                raise RBridgeTransportError(
                    "R worker is not running. Call start() first."
                )

            # Never reuse a pipe after a timed-out reader may still hold stdout
            if self._reader_alive:
                logger.warning("Discarding desynced R worker before next query")
                self.stop()
                raise RBridgeTransportError(
                    "R worker was desynced after a prior timeout; restart required."
                )

            request = {"tool": tool}
            if params:
                request["params"] = params

            request_json = json.dumps(request)
            logger.debug("Sending to R: %s", request_json[:200])

            try:
                self._process.stdin.write(request_json + "\n")
                self._process.stdin.flush()
            except (BrokenPipeError, OSError) as e:
                self.stop()
                raise RBridgeTransportError(
                    f"Failed to write to R process: {e}"
                ) from e

            response = self._read_response(timeout=self.query_timeout)

            if response is None:
                # Timeout or unreadable line: kill worker so the next start()
                # cannot read a late line as the wrong response.
                stderr_tail = "\n".join(self._stderr_lines[-10:])
                crashed = not self.is_running
                self.stop()
                if crashed:
                    raise RBridgeTransportError(
                        f"R worker crashed during query '{tool}'.\n"
                        f"Stderr:\n{stderr_tail}"
                    )
                raise RBridgeTransportError(
                    f"R worker timed out after {self.query_timeout}s on tool '{tool}'"
                )

            if response.get("error"):
                raise RBridgeDomainError(
                    response.get("message", "Unknown error from R worker")
                )

            return response.get("result")

    def ping(self) -> bool:
        """Health check — returns True if R worker responds."""
        try:
            result = self.query("ping")
            return result.get("status") == "ok"
        except RBridgeError:
            return False

    def _read_response(self, timeout: float) -> dict | None:
        """Read a single JSON response line from R stdout with timeout."""
        if not self._process or not self._process.stdout:
            return None

        result_holder: list[str | None] = [None]

        def _reader():
            try:
                line = self._process.stdout.readline()
                result_holder[0] = line
            except Exception:
                result_holder[0] = None

        reader_thread = threading.Thread(target=_reader, daemon=True)
        reader_thread.start()
        reader_thread.join(timeout=timeout)

        if reader_thread.is_alive():
            # Timeout — thread still blocked on readline; mark desync and
            # let caller stop() the process so this thread dies with the pipe.
            self._reader_alive = True
            return None

        self._reader_alive = False
        line = result_holder[0]
        if not line:
            return None

        line = line.strip()
        if not line:
            return None

        try:
            return json.loads(line)
        except json.JSONDecodeError as e:
            logger.warning("Invalid JSON from R: %s (error: %s)", line[:100], e)
            # Protocol corruption: treat as transport failure
            self._reader_alive = True
            return None

    def _read_stderr(self) -> None:
        """Background thread: read stderr for logging."""
        try:
            while self._process and self._process.stderr:
                line = self._process.stderr.readline()
                if not line:
                    break
                line = line.rstrip("\n")
                self._stderr_lines.append(line)
                # Keep last 100 lines
                if len(self._stderr_lines) > 100:
                    self._stderr_lines = self._stderr_lines[-50:]
                logger.debug("[R] %s", line)
        except Exception:
            pass

    def __enter__(self):
        self.start()
        return self

    def __exit__(self, *args):
        self.stop()

    def __del__(self):
        try:
            self.stop()
        except Exception:
            pass
