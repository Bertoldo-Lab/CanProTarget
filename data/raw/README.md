# Raw chemoproteomics rebuild inputs

These files rebuild the bundled Chemistry RDS outputs via
`docs/scripts/preprocess_protein_binding.R` and
`docs/scripts/preprocess_data.R`.

| File | Role |
|------|------|
| `table-s2.xlsx` | CysDB ligandable matrix + Compound Keys (Git LFS, ~124 MB) |
| `id_mapping.tsv` | Protein ID → gene symbol |
| `swissadme.csv` | SwissADME batch export |
| `_archive/backus_swiss_link.xlsx` | **Archived**. Not used by the app or the current build. Kept for provenance. See `_archive/WHY_ARCHIVED.md`. |

## What `table-s2.xlsx` has to contain

`docs/scripts/preprocess_protein_binding.R` reads two sheets by name:

- **`Ligandable Dataset`** — one row per cysteine, in wide form: metadata
  columns `proteinid`, `cysteineid`, `resid`, `ligandable`, then one column
  per probe holding that probe's competition ratio. Cells may be blank or
  `--`; those are dropped, and the remainder are coerced to numeric.
- **`Compound Keys`** — the probe-to-structure key. Supplies SMILES,
  Dataset and Cell_Line for every probe.

## How ACRYL and CL probes are renumbered

The `Compound_Name` column in `table-s2.xlsx` uses an indexing that does
**not** match what the public **CysDB web app**
(<https://backuslab.shinyapps.io/cysdb/>) shows for the same molecule.
CanProTarget follows the web app's naming so users cross-referencing our
probe labels against CysDB get the intended compound. The mapping is
applied in `docs/scripts/preprocess_protein_binding.R` §5.

| S2 raw `Compound_Name` | app name | Notes |
|---|---|---|
| `CL_1` .. `CL_19`        | `CL1` .. `CL19` | unchanged |
| `CL_20+`                 | `CL(N−1)`       | shift −1 (e.g. S2 `CL_175` → app `CL174`) |
| `ACRYL_1..3`             | `AC0..AC2`      | shift −1 |
| `ACRYL_4`                | `AC4`           | unchanged (band boundary) |
| `ACRYL_5, 7..11`         | `AC(N−2)`       | shift −2 |
| `ACRYL_6`                | `AC5`           | explicit override so it lands beside its replicate `ACRYL_7` |
| `ACRYL_12`               | `OTHER_6`       | carve-out |
| `ACRYL_13`               | `OTHER_7`       | carve-out |
| `ACRYL_14+`              | `AC(N−4)`       | shift −4 |

**Worked example.** The molecule the CysDB web app labels `CL_174` is
stored in the S2 spreadsheet as `Compound_Name = CL_175`, with SMILES
`ClCC(=O)N(C1CC1)C2=CCCCC2`. The app displays it as `CL174`. To find
our `CL174` in `table-s2.xlsx`, search for `CL_175`.

**Replicate names.** When the same molecule was tested in two cell lines
and CysDB assigned it two adjacent `Compound_Name` values (e.g.
`ACRYL_1` in MDA-MB-231 and `ACRYL_2` in Ramos), the app carries them as
two adjacent probe names (`AC0` and `AC1`). Both share the same
Canonical.SMILES and ADME properties; `Cell_Line` distinguishes the
replicates. This is visible on the Chemistry Explorer as two adjacent
dropdown entries with identical structure images — by design.

**Collision guard.** The rename script refuses to build if two rows with
different SMILES would collapse to the same app name. Future
contributors cannot silently reintroduce a merge.

## Git LFS

`table-s2.xlsx` is stored with **Git LFS**. After clone:

```bash
git lfs install
git lfs pull
```

The app never reads these files at runtime, so a clone without Git LFS
still runs. `docs/scripts/fetch_data_assets.sh` offers a release-based
fallback via the `data-assets-v1` release. See `MISSING_DATA.md`.
