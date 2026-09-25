# Raw chemoproteomics rebuild inputs

These files rebuild the bundled Chemistry RDS outputs via
`docs/scripts/preprocess_protein_binding.R` and
`docs/scripts/preprocess_data.R`.

| File | Role |
|------|------|
| `table-s2.xlsx` | CysDB ligandable matrix + Compound Keys (Git LFS, ~124 MB) |
| `id_mapping.tsv` | Protein ID → gene symbol |
| `swissadme.csv` | SwissADME batch export |

## What `table-s2.xlsx` has to contain

`docs/scripts/preprocess_protein_binding.R` reads two sheets by name:

- **`Ligandable Dataset`** — one row per cysteine, in wide form: metadata
  columns `proteinid`, `cysteineid`, `resid`, `ligandable`, then one column
  per probe holding that probe's competition ratio. Cells may be blank or
  `--`; those are dropped, and the remainder are coerced to numeric.
- **`Compound Keys`** — the probe-to-structure key. Supplies SMILES,
  Dataset and Cell_Line for every probe.

## Probe naming

The `Compound_Name` column in `table-s2.xlsx` uses an indexing that does
not match what the public CysDB web app shows for the same molecule.
Rebuilders do not need to worry about this — the preprocess script
(`docs/scripts/preprocess_protein_binding.R` §5) applies the mapping and a
collision guard that refuses to build if two different structures would
collapse to the same app name. The full rename tables, collision-guard
rationale and a worked example live in
[`docs/DATA_PROVENANCE.md`](../../docs/DATA_PROVENANCE.md).

## Git LFS

`table-s2.xlsx` is stored with **Git LFS**. After clone:

```bash
git lfs install
git lfs pull
```

The app never reads these files at runtime, so a clone without Git LFS
still runs. `docs/scripts/fetch_data_assets.sh` offers a release-based
fallback via the `data-assets-v1` release. See `REBUILD.md`.
