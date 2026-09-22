# ============================================================
# Script: preprocess_protein_binding.R
# Purpose: Offline preprocessing for protein lookup
#          - Calculate n_targets for each probe (CR >= 4)
#          - Map proteinid to gene_name
#          - Rename ACRYL/CL probes according to legacy logic
#          - Output lookup RDS for Shiny (CR == 0 dropped; gene_name_key added)
#
# IMPORTANT:
#   Lives under docs/scripts/ so it is not auto-sourced by Shiny (only `R/*.R` are).
#   Run only via docs/scripts/preprocess_data.R or directly with Rscript.
# ============================================================

library(readxl)
library(dplyr)
library(tidyr)

# Configuration
# Rebuild inputs live under data/raw/ (see data/raw/README.md). These paths
# pointed at data/ itself, where the files have not been since they were moved,
# so the script could not run at all.
INPUT_FILE <- "data/raw/table-s2.xlsx"
ID_MAPPING_FILE <- "data/raw/id_mapping.tsv"
OUTPUT_FILE <- "data/protein_binding_lookup_preprocessed.rds"
CR_THRESHOLD <- 4  # Define what counts as a "target"

cat("Starting protein binding preprocessing...\n")

# ============================================================
# 1. Load ID mapping (proteinid -> gene_name)
# ============================================================
cat("Loading ID mapping...\n")
id_mapping <- read.delim(ID_MAPPING_FILE, header = TRUE, stringsAsFactors = FALSE) %>%
  rename(proteinid = From, gene_name = To)

# ============================================================
# 2. Load Ligandable Dataset
# ============================================================
cat("Loading Ligandable Dataset...\n")
ligandable_data <- read_excel(INPUT_FILE, sheet = "Ligandable Dataset")

# ============================================================
# 3. Load Compound Keys
# ============================================================
cat("Loading Compound Keys...\n")
compound_keys <- read_excel(INPUT_FILE, sheet = "Compound Keys")

# ============================================================
# 4. Transform ligandable data from wide to long format
# ============================================================
cat("Reshaping data to long format...\n")

# Identify probe columns (exclude metadata columns)
metadata_cols <- c("proteinid", "cysteineid", "resid", "ligandable")
probe_cols <- setdiff(colnames(ligandable_data), metadata_cols)

# Convert all probe columns to character first (handles mixed types)
ligandable_data <- ligandable_data %>%
  mutate(across(all_of(probe_cols), as.character))

# Pivot to long format
protein_binding_long <- ligandable_data %>%
  pivot_longer(
    cols = all_of(probe_cols),
    names_to = "probe_name_raw",
    values_to = "CR"
  ) %>%
  # Remove missing or invalid CR values
  filter(!is.na(CR), CR != "--", CR != "", trimws(CR) != "") %>%
  # Convert to numeric after filtering out non-numeric values
  mutate(CR = suppressWarnings(as.numeric(CR))) %>%
  # Remove any rows where conversion to numeric failed
  filter(!is.na(CR))

# ============================================================
# 5. Apply ACRYL/CL renaming logic (from old functions.R)
# ============================================================
cat("Renaming probes (ACRYL/CL logic)...\n")

probe_name_raw_lookup <- NULL   # filled below, used by the collision guard

protein_binding_long <- protein_binding_long %>%
  mutate(
    probe_num = suppressWarnings(as.numeric(sub("ACRYL_|CL_|OTHER_", "", probe_name_raw))),
    probe_type = case_when(
      startsWith(probe_name_raw, "ACRYL_") ~ "ACRYL",
      startsWith(probe_name_raw, "CL_") ~ "CL",
      startsWith(probe_name_raw, "OTHER_") ~ "OTHER",
      TRUE ~ "OTHER"
    ),
    probe_name = case_when(
      # ACRYL logic
      #
      # ACRYL_6 is named explicitly because the band below would send it to
      # ACRYL_4, which ACRYL_4 itself already occupies. Checked against the
      # Compound Keys sheet: ACRYL_4 is publication compound 14 and ACRYL_6 is
      # compound 15, so that band merged two different structures and pooled
      # their competition ratios. Compound 15 is SwissADME's AC5, where its
      # other replicate ACRYL_7 already lands, so ACRYL_6 belongs beside it --
      # one name per compound, with Cell_Line telling the replicates apart,
      # exactly as CL_19 and CL_20 are handled.
      probe_type == "ACRYL" & probe_num == 6 ~ "ACRYL_5",
      probe_type == "ACRYL" & probe_num >= 1 & probe_num <= 3 ~ paste0("ACRYL_", probe_num - 1),
      probe_type == "ACRYL" & probe_num >= 5 & probe_num <= 11 ~ paste0("ACRYL_", probe_num - 2),
      probe_type == "ACRYL" & probe_num == 12 ~ "OTHER_6",
      probe_type == "ACRYL" & probe_num == 13 ~ "OTHER_7",
      probe_type == "ACRYL" & probe_num >= 14 ~ paste0("ACRYL_", probe_num - 4),

      # CL logic
      probe_type == "CL" & probe_num >= 1 & probe_num <= 19 ~ probe_name_raw,
      probe_type == "CL" & probe_num >= 20 ~ paste0("CL_", probe_num - 1),

      # Default fallback (no changes)
      TRUE ~ probe_name_raw
    )
  )

probe_name_raw_lookup <- unique(
  protein_binding_long[, c("probe_name_raw", "probe_name")]
)
names(probe_name_raw_lookup) <- c("raw", "out")

protein_binding_long <- protein_binding_long %>%
  select(-probe_num, -probe_type, -probe_name_raw)

# ============================================================
# 5b. A merge is only allowed between replicates of one compound
# ============================================================
# Two source columns may legitimately share an output name when they are the
# same compound measured in different cell lines -- CL_19 and CL_20 are one
# structure in MDA-MB-231 and Ramos, and Cell_Line still tells them apart.
# What must never happen is two different structures sharing a name, which
# pools the competition ratios of unrelated molecules. Compound Keys carries
# the SMILES, so the difference is checkable rather than a matter of trust.
local({
  key <- unique(data.frame(
    raw = probe_name_raw_lookup$raw,
    out = probe_name_raw_lookup$out,
    stringsAsFactors = FALSE
  ))
  key$smiles <- compound_keys$SMILES[match(key$raw, compound_keys$Compound_Name)]
  merged <- key[key$out %in% key$out[duplicated(key$out)], , drop = FALSE]
  if (nrow(merged)) {
    n_struct <- tapply(merged$smiles, merged$out, function(x) length(unique(x)))
    bad <- names(n_struct)[n_struct > 1]
    if (length(bad)) {
      clash <- merged[merged$out %in% bad, , drop = FALSE]
      clash <- clash[order(clash$out, clash$raw), , drop = FALSE]
      stop("Probe renaming merges different structures:\n",
           paste(sprintf("  %s -> %s", clash$raw, clash$out), collapse = "\n"),
           "\nCheck these against the Compound Keys sheet before rebuilding.",
           call. = FALSE)
    }
    message("  note: ", length(unique(merged$out)),
            " name(s) carry replicates of one compound across cell lines")
  }
})

# ============================================================
# 6. Map proteinid to gene_name
# ============================================================
cat("Mapping proteinid to gene_name...\n")

protein_binding_long <- protein_binding_long %>%
  left_join(id_mapping, by = "proteinid") %>%
  # Remove rows without gene_name mapping
  filter(!is.na(gene_name))

# ============================================================
# 7. Join with compound keys to get Dataset & Cell_Line
# ============================================================
cat("Joining with compound keys...\n")

# Create lookup from compound keys
compound_lookup <- compound_keys %>%
  select(Compound_Name, Dataset, Cell_Line, SMILES) %>%
  rename(probe_name = Compound_Name)

protein_binding_long <- protein_binding_long %>%
  left_join(compound_lookup, by = "probe_name")

# ============================================================
# 8. Calculate n_targets for each probe (CR >= threshold)
# ============================================================
cat("Calculating n_targets per probe (CR >= ", CR_THRESHOLD, ")...\n")

# Count distinct proteins per probe where CR >= 4
target_counts <- protein_binding_long %>%
  filter(CR >= CR_THRESHOLD) %>%
  group_by(probe_name) %>%
  summarise(n_targets = n_distinct(proteinid), .groups = "drop")

# Add n_targets to all rows
protein_binding_final <- protein_binding_long %>%
  left_join(target_counts, by = "probe_name") %>%
  # If a probe has no targets with CR >= 4, set n_targets = 0
  mutate(n_targets = ifelse(is.na(n_targets), 0, n_targets))

# ============================================================
# 9. Select and order final columns
# ============================================================
cat("Finalizing data structure...\n")

protein_binding_final <- protein_binding_final %>%
  select(
    probe_name,
    proteinid,
    gene_name,
    CR,
    n_targets,
    cysteineid,
    ligandable,
    Dataset,
    Cell_Line,
    SMILES
  ) %>%
  arrange(probe_name, gene_name)

# ============================================================
# 9b. Optional: align SMILES with SwissADME (reduces work in Shiny)
#     Build data/swissadme_preprocessed.rds first (SwissADME preprocessing).
# ============================================================
swiss_rds <- file.path(dirname(OUTPUT_FILE), "swissadme_preprocessed.rds")
if (file.exists(swiss_rds)) {
  functions_r <- normalizePath(
    file.path(dirname(OUTPUT_FILE), "..", "R", "functions.R"),
    mustWork = FALSE
  )
  if (!is.na(functions_r) && file.exists(functions_r)) {
    source(functions_r, local = FALSE)
    protein_binding_final <- overlay_canonical_smiles_from_swissadme(
      protein_binding_final,
      readRDS(swiss_rds)
    )
    cat("Aligned SMILES with SwissADME (swissadme_preprocessed.rds).\n")
  }
}

# ============================================================
# 10. Drop CR == 0 (not used in app), add gene_name_key for lookups
# ============================================================
n_before <- nrow(protein_binding_final)
protein_binding_final <- protein_binding_final %>%
  dplyr::filter(is.na(.data$CR) | .data$CR != 0) %>%
  dplyr::mutate(gene_name_key = tolower(as.character(.data$gene_name)))
cat(
  "Dropped CR == 0 rows:",
  n_before - nrow(protein_binding_final),
  "of",
  n_before,
  "\n"
)

# ============================================================
# 11. Save to RDS (fast loading for Shiny)
# ============================================================
cat("Saving to RDS...\n")
# xz, deliberately. Unlike the DepMap matrices, reading this table is bound by
# deserialising ten million rows rather than by decompression: gzip saved under
# two seconds and cost more than a hundred megabytes in the deploy bundle. The
# app reads the factored copy at about a second either way.
saveRDS(protein_binding_final, OUTPUT_FILE, compress = "xz")

# ============================================================
# 12. Summary statistics
# ============================================================
cat("\n=== Preprocessing Complete ===\n")
cat("Output file:", OUTPUT_FILE, "\n")
cat("Total rows:", nrow(protein_binding_final), "\n")
cat("Unique probes:", n_distinct(protein_binding_final$probe_name), "\n")
cat("Unique proteins:", n_distinct(protein_binding_final$proteinid), "\n")
cat("Unique genes:", n_distinct(protein_binding_final$gene_name), "\n")
cat("\nProbes with most targets (CR >= 4):\n")
print(
  protein_binding_final %>%
    select(probe_name, n_targets) %>%
    distinct() %>%
    arrange(desc(n_targets)) %>%
    head(10)
)

cat("\nSample of preprocessed data:\n")
print(head(protein_binding_final, 20))

