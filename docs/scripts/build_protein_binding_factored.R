# ============================================================
# Script:   build_protein_binding_factored.R
# Purpose:  Store the protein-binding table with its character columns as
#           factors, which the app loads in place of the original.
#
#   Usage:  Rscript docs/scripts/build_protein_binding_factored.R
#
# Why this exists
#
#   The table is 10,588,541 rows. Every character column costs about 85 MB
#   whatever it contains, because R keeps one pointer per element — the
#   "ligandable" column has two distinct values and still costs 85 MB. As
#   factors the same data is an integer vector plus a short level table.
#
#     original   936 MB in RAM,  22.8 s to load
#     factored   555 MB in RAM,   1.8 s to load
#
#   Behaviour is unchanged: comparisons, %in%, match() and split() treat
#   factors the same, and the two places that do string work
#   (simplify_probe_name and the probe index) call as.character() first.
#
#   The app prefers the factored file and falls back to converting the
#   original on the fly, so this script is an optimisation, not a dependency.
# ============================================================

src <- file.path("data", "protein_binding_lookup_preprocessed.rds")
dst <- file.path("data", "protein_binding_lookup_factored.rds")
if (!file.exists(src)) stop("Missing ", src, " — see docs/DATA_PROVENANCE.md")

message("Reading ", src, " ...")
x <- readRDS(src)
before <- as.numeric(object.size(x)) / 1e6

chr <- vapply(x, is.character, logical(1))
message("Converting ", sum(chr), " character columns to factors ...")
if (any(chr)) x[chr] <- lapply(x[chr], as.factor)
after <- as.numeric(object.size(x)) / 1e6

message("Writing ", dst, " (xz; this takes a minute) ...")
# xz, deliberately. Unlike the DepMap matrices, reading this table is bound by
# deserialising ten million rows rather than by decompression: gzip saved under
# two seconds and cost more than a hundred megabytes in the deploy bundle. The
# app reads the factored copy at about a second either way.
saveRDS(x, dst, compress = "xz")

message(sprintf("Done. %.0f MB -> %.0f MB in memory (%.0f%% smaller); %.0f MB on disk.",
                before, after, 100 * (1 - after / before), file.size(dst) / 1e6))
