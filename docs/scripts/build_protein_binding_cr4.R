# ============================================================
# Script:   build_protein_binding_cr4.R
# Purpose:  Engagement-only subset of the chemoproteomics binding table.
#
#   Usage:  Rscript docs/scripts/build_protein_binding_cr4.R
#
# The full table is 10,588,541 probe x cysteine measurements, of which 30,219
# (0.3%) reach the CR >= 4 engagement threshold the app and the manuscript use.
# The other 99.7% are measurements showing no displacement -- necessary to know
# what was tested, but never shown unless someone lowers the CR floor below 4.
#
# Holding all of it in memory to serve 0.3% costs 555 MB and a visible wait on
# first use. This writes the engaged rows to their own file so the default
# views load from something small, and the full table is read only when the CR
# floor actually goes below 4.
# ============================================================

suppressMessages({
  library(dplyr)
})

data_dir <- "data"
src <- file.path(data_dir, "protein_binding_lookup_factored.rds")
if (!file.exists(src)) src <- file.path(data_dir, "protein_binding_lookup_preprocessed.rds")
stopifnot(file.exists(src))
out_path <- file.path(data_dir, "protein_binding_cr4.rds")

msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "-", ..., "\n")

msg("Reading", basename(src), "...")
t0 <- Sys.time()
pb <- readRDS(src)
msg(sprintf("  %s rows in %.1f s", format(nrow(pb), big.mark = ","),
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))

cr <- suppressWarnings(as.numeric(as.character(pb$CR)))
keep <- !is.na(cr) & cr >= 4
engaged <- pb[keep, , drop = FALSE]

# Factors carried over from the parent keep every unused level, which costs
# more than the rows themselves at this size.
fac <- vapply(engaged, is.factor, logical(1))
engaged[fac] <- lapply(engaged[fac], droplevels)

msg(sprintf("  keeping %s engaged rows (%.2f%%)",
            format(nrow(engaged), big.mark = ","), 100 * mean(keep)))
saveRDS(engaged, out_path)

msg(sprintf("Wrote %s: %.1f MB on disk, %.1f MB in memory (parent: %.1f MB / %.1f MB)",
            basename(out_path),
            file.size(out_path) / 1e6,
            as.numeric(object.size(engaged)) / 1e6,
            file.size(src) / 1e6,
            as.numeric(object.size(pb)) / 1e6))

t1 <- Sys.time(); invisible(readRDS(out_path))
msg(sprintf("Read back in %.2f s", as.numeric(difftime(Sys.time(), t1, units = "secs"))))
