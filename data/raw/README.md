# Raw chemoproteomics rebuild inputs

These files rebuild the bundled Chemistry RDS outputs
via `docs/scripts/preprocess_protein_binding.R` and
`docs/scripts/preprocess_data.R`.

| File | Role |
|------|------|
| `table-s2.xlsx` | Backus ligandable matrix + Compound Keys, restructured (Git LFS) — see below |
| `id_mapping.tsv` | Protein ID → gene symbol |
| `swissadme.csv` | SwissADME batch export |
| `backus_swiss_link.xlsx` | Probe name mapping for SwissADME (optional) |

## What `table-s2.xlsx` has to contain

`docs/scripts/preprocess_protein_binding.R` reads two sheets by name:

- **`Ligandable Dataset`** — one row per cysteine, in wide form: the metadata
  columns `proteinid`, `cysteineid`, `resid` and `ligandable`, then one column
  per probe holding that probe's competition ratio. Every other column is
  treated as a probe. Cells may be blank or `--`; those are dropped, and the
  remainder are coerced to numeric. The script pivots this to one row per
  cysteine and probe.
- **`Compound Keys`** — the probe-to-structure key, which supplies SMILES.

This layout is not the one the study was published in: the supplementary
table was restructured into these two sheets before it was committed here.
**The steps taken to restructure it are not recorded**, so the file cannot
currently be regenerated from the publication.

One consequence was found here and has been fixed in the script; the shipped
RDS still carries it until the chemoproteomics table is rebuilt.

The 1,000 probe columns are one per compound and cell line. `Compound Keys`
gives each its publication identifier and SMILES, and the preprocessing script
renumbers the columns to the scheme SwissADME uses. Checking every probe by
heavy-atom formula against the SwissADME table, **997 of 998 agreed**; the
numbering is sound. The exception was a band boundary that sent `ACRYL_6` onto
`ACRYL_4`, a name `ACRYL_4` already held:

| Source | Publication compound | Structure |
|---|---|---|
| `ACRYL_4` | 14 | C11H7F6NO — SwissADME `AC4` |
| `ACRYL_6` | 15 | C22H17F3N2O2 — SwissADME `AC5` |

Two different molecules under one name, pooling 4,096 records. Compound 15's
other replicate `ACRYL_7` already mapped to `AC5`, so `ACRYL_6` is now named
explicitly to join it — one name per compound, `Cell_Line` telling replicates
apart, the same way `CL_19` and `CL_20` are handled. With that rule in place
all **998 of 998** probes agree with SwissADME by formula.

`CL_19` and `CL_20` share a name legitimately: one structure, publication
compound 13, in MDA-MB-231 and Ramos. The script's guard allows a merge only
when the source columns share a SMILES and stops the build otherwise.

**Status: corrected here, not upstream.** The reassignment of `ACRYL_6` to
`AC5` is a decision taken by the CanProTarget authors on the evidence in the
`Compound Keys` sheet. The authors of the source study were not consulted and
have not reviewed it, and `table-s2.xlsx` is unchanged. The evidence is that
`Compound Keys` gives `ACRYL_4` and `ACRYL_6` different structures,
publication compounds 14 and 15, and that compound 15's other replicate
`ACRYL_7` already mapped to `AC5`; with `ACRYL_6` beside it all 998 probes
agree with the SwissADME table by heavy-atom formula, against 997 before.
`AC4` now holds 2,859 records, all compound 14, and `AC5` holds 2,119, being
compound 15's two cell-line replicates; before the change `AC4` held 4,096
spanning both compounds. Anyone rebuilding from the published supplementary
table should be aware this differs from the numbering as distributed.

**The rebuild is required for the fix to reach the app.** `ACRYL_4` and
`ACRYL_6` are both MDA-MB-231, so the pooled rows in the existing
`protein_binding_lookup_preprocessed.rds` cannot be told apart after the fact
and cannot be patched in place.

Both `.xlsx` files are stored with **Git LFS** (`table-s2.xlsx` is ~124 MB).
After clone:

```bash
git lfs install
git lfs pull
```

These are rebuild inputs only — the app never reads them at runtime, so a clone
without Git LFS still works. `docs/scripts/fetch_data_assets.sh` offers a
release-based fallback via the `data-assets-v1` release (`chemoproteomics_raw.tar.gz`
only). Provenance is in `docs/DATA_PROVENANCE.md`.
