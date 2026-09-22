#!/usr/bin/env bash
# Fallback for clones without Git LFS: download the raw chemoproteomics
# rebuild inputs from a GitHub Release and extract into data/raw/.
#
# Usage (from project root):
#   docs/scripts/fetch_data_assets.sh [tag]
#
# Default tag: data-assets-v1
#
# The release exists: Bertoldo-Lab/CanProTarget @ data-assets-v1.
# It now carries only chemoproteomics_raw.tar.gz.
# precomputed_stats.tar.gz was deleted 20 Aug 2026 (Excel-truncated CRISPR cache).
# Do not re-upload it. Build a local cache with:
#   Rscript docs/scripts/precompute_effectsizes_all_subtypes.R
#
# Prefer `git lfs pull` when Git LFS is available.
#
# !! This script downloads ONLY chemoproteomics_raw.tar.gz. !!

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"

TAG="${1:-data-assets-v1}"
REPO="${GITHUB_REPO:-Bertoldo-Lab/CanProTarget}"

need() { command -v "$1" >/dev/null 2>&1 || { echo "Missing dependency: $1" >&2; exit 1; }; }
need gh
need tar

if ! gh release view "$TAG" -R "$REPO" >/dev/null 2>&1; then
  echo "Release '${TAG}' not found in ${REPO}." >&2
  echo "Use 'git lfs pull' instead, or create the release (see header)." >&2
  exit 1
fi

mkdir -p data/raw tmp_data_assets
cd tmp_data_assets

echo "Downloading release assets for ${REPO}@${TAG} ..."
gh release download "$TAG" -R "$REPO" --clobber -p 'chemoproteomics_raw.tar.gz'

echo "Extracting chemoproteomics raw sources ..."
tar -xzf chemoproteomics_raw.tar.gz -C "$ROOT/data/raw"

cd "$ROOT"
rm -rf tmp_data_assets

echo "Done."
echo "  data/raw/ -> $(ls data/raw | tr '\n' ' ')"
