# Release checklist — CanProTarget v1.0.0

Tag `v1.0.0` is already on `main`. This file is only the remaining archive
step. Delete it once the Zenodo DOIs are written into the repository.

## Waiting on Zenodo

Zenodo mints the DOI when a GitHub Release is published, and only if this
repository was already connected. A DOI cannot be reserved beforehand on
that path.

1. Make this repository public: Settings → General → Danger zone → Change
   visibility. Leave `Bertoldo-Lab/CanProTarget_dev` private.
2. Sign in to [zenodo.org](https://zenodo.org) with GitHub. Open
   [GitHub settings](https://zenodo.org/account/settings/github/) and enable
   `Bertoldo-Lab/CanProTarget`. If the repository is missing from that list,
   a Bertoldo-Lab owner approves the Zenodo GitHub application, then
   refreshes the page.
3. Publish the GitHub Release from the existing tag `v1.0.0`. Do not create
   another tag and do not attach files. Do not use `data-assets-v1`. That
   release is only `chemoproteomics_raw.tar.gz`, the chemoproteomics rebuild
   inputs.
4. Zenodo then issues two DOIs:
   - **Version DOI** for `v1.0.0`. Cite this in the manuscript.
   - **Concept DOI**, which always resolves to the newest version. Put this
     in `CITATION.cff`, the README citation, and the About page in `app.R`.
5. The `v1.0.0` archive will still say "DOI pending release", because the
   DOI did not exist when the tag was cut. After those three files carry the
   concept DOI, tag `v1.0.1` and publish that release the same way. Then
   delete this file.

## Already settled

- `main` is one commit. Development history stays in
  `Bertoldo-Lab/CanProTarget_dev`.
- Version `1.0.0` is set in `R/api_functions.R`, `data/data_versions.yaml`
  and `CITATION.cff`.
- There is no `renv.lock`.
- `data/raw/*.xlsx` is Git LFS and marked `export-ignore`, so a release zip
  does not contain stub workbooks. The app never reads those files.
- The manuscript states that this archive carries the code and that the
  dependency matrices come from DepMap. What a clone still fetches is in
  `data/README.md`.
- How to reproduce the benchmark, and which inputs sit outside the
  repository, is in the README section "Reproducing the manuscript
  benchmark".
