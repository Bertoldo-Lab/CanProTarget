# Usage and diagnostics logging

## Is file-based usage logging a good idea?

**For local development or a server you control:** yes, append-only logs under `logs/usage/` are simple and keep the app code small. You can aggregate JSON Lines with standard tools (`jq`, R, Python).

**For shinyapps.io:** the container filesystem is **ephemeral**. Anything written under the app directory during a session is **not** a durable audit trail unless you copy it elsewhere. For production analytics you typically use:

- shinyapps.io **Metrics** and **Logs** in the admin UI  
- An external service (Plausible, Google Analytics with privacy settings, your own API endpoint, etc.)

**Privacy:** logging free-text queries (e.g. protein names) can be **personally or commercially sensitive**. The default logger only records **counts and parameters** unless you extend it. If you add raw search strings, document retention and comply with your institution’s policy.

## Optional file logging in this repo

1. Set the environment variable **`CPT_LOG_USAGE=1`** before starting the app (or in shinyapps.io: Settings → Environment).

2. Logs are appended as **one JSON object per line** under `logs/usage/usage_YYYY-MM-DD.jsonl`.

3. Log files are **gitignored** (`logs/usage/*.jsonl`) so they do not bloat the repository or deployment bundle.

4. Toggle in code: `app_config$logging$enabled <- TRUE` in `R/app_config.R` (still respects `CPT_LOG_USAGE` unless you set `enabled = TRUE` for local-only testing; see `R/usage_logger.R`).

## Events currently logged

| Event | Module | Fields (examples) |
|-------|--------|-------------------|
| `deps_run_analysis` | Discover | `dataset`, `subtype`, `n_subtype` |
| `deps_load_probes` | Discover | `trigger` (`layer` or `cr_below_4`) |
| `swiss_load_binding` | Chemistry | `probe`, `cr_cutoff` |

## Crashes (e.g. OOM)

R **Out-of-memory** kills usually **do not** flush a custom log line. Rely on **shinyapps.io Logs** right after a disconnect, and reduce peak memory (subset matrices, lazy load, `gc()` after large objects) as described in deployment docs.
