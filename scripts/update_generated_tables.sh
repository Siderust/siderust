#!/usr/bin/env bash
# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Vallés Puig, Ramon
#
# update_generated_tables.sh
#
# Refreshes the committed Rust tables in src/generated/ by re-running the
# full download + codegen pipeline for VSOP87, ELP2000, and IERS EOP.
#
# Usage:
#   ./scripts/update_generated_tables.sh [--datasets-dir DIR]
#
# Options:
#   --datasets-dir DIR   Where to download/cache raw source files.
#                        Defaults to $SIDERUST_DATASETS_DIR, or
#                        <repo_root>/.siderust_datasets if unset.
#
# Environment variables:
#   SIDERUST_DATASETS_DIR   Alternative to --datasets-dir.
#
# After this script completes, run:
#   git diff dev-deps/siderust/src/generated/
# to review changes, then commit them.

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

DATASETS_DIR=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --datasets-dir)
      DATASETS_DIR="$2"
      shift 2
      ;;
    *)
      echo "Unknown option: $1" >&2
      exit 1
      ;;
  esac
done

# Resolve datasets directory.
if [[ -z "$DATASETS_DIR" ]]; then
  DATASETS_DIR="${SIDERUST_DATASETS_DIR:-"$ROOT_DIR/.siderust_datasets"}"
fi

echo "==> Downloading source datasets into: $DATASETS_DIR"
SIDERUST_DATASETS_DIR="$DATASETS_DIR" \
  bash "$ROOT_DIR/scripts/prefetch_datasets.sh" --minimal

echo ""
echo "==> Regenerating Rust tables from downloaded datasets..."
echo "    Output: $ROOT_DIR/src/generated/"

SIDERUST_REGEN=1 \
SIDERUST_DATASETS_DIR="$DATASETS_DIR" \
  cargo build --manifest-path "$ROOT_DIR/Cargo.toml" 2>&1

echo ""
echo "==> Writing datasets.lock.json..."

GENERATED_DIR="$ROOT_DIR/src/generated"
mkdir -p "$GENERATED_DIR"

TIMESTAMP="$(date -u +%Y-%m-%dT%H:%M:%SZ)"

# Compute SHA-256 of each output file.
sha256_of() {
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$1" | awk '{print $1}'
  else
    shasum -a 256 "$1" | awk '{print $1}'
  fi
}

VSOP87A_SHA="$(sha256_of "$GENERATED_DIR/vsop87a.rs")"
VSOP87E_SHA="$(sha256_of "$GENERATED_DIR/vsop87e.rs")"
ELP_SHA="$(sha256_of "$GENERATED_DIR/elp_data.rs")"
IERS_SHA="$(sha256_of "$GENERATED_DIR/iers_eop_data.rs")"

cat > "$GENERATED_DIR/datasets.lock.json" <<EOF
{
  "generated_at": "${TIMESTAMP}",
  "sources": {
    "vsop87": {
      "url": "https://ftp.imcce.fr/pub/ephem/planets/vsop87/",
      "output": "vsop87a.rs",
      "sha256": "${VSOP87A_SHA}"
    },
    "vsop87e": {
      "url": "https://ftp.imcce.fr/pub/ephem/planets/vsop87/",
      "output": "vsop87e.rs",
      "sha256": "${VSOP87E_SHA}"
    },
    "elp2000": {
      "url": "https://cdsarc.cds.unistra.fr/ftp/VI/79/",
      "output": "elp_data.rs",
      "sha256": "${ELP_SHA}"
    },
    "iers_eop": {
      "url": "https://datacenter.iers.org/products/eop/rapid/standard/finals2000A.all",
      "output": "iers_eop_data.rs",
      "sha256": "${IERS_SHA}"
    }
  }
}
EOF

echo ""
echo "==> Done. Updated files:"
ls -lh "$GENERATED_DIR"
echo ""
echo "Review changes with:"
echo "  git diff dev-deps/siderust/src/generated/"
echo ""
echo "Then commit:"
echo "  git add dev-deps/siderust/src/generated/"
echo "  git commit -m 'chore: refresh generated dataset tables'"
