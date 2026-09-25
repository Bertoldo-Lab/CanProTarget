# Why these files are archived

## backus_swiss_link.xlsx

Archived as part of the v1.0.1 fix (protein_binding SMILES-annotation drift).

Historically this workbook carried an `updated_compound_names` column that
was used by the SwissADME preprocess flow (`docs/scripts/reference_swissadme_from_csv.R`)
to attach CysDB Compound_Names to each SwissADME `Molecule N`. From v1.0.1
onwards the pipeline reads CysDB `Compound Keys` (`data/raw/table-s2.xlsx`)
directly and does not need this link file at all — the app's probe renaming
already lives in `docs/scripts/preprocess_protein_binding.R` (section 5).

Do not re-introduce this file into the build pipeline. It is kept here only
for provenance of pre-v1.0.1 rebuilds.
