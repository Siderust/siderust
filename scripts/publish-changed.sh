#!/usr/bin/env bash
# Publish workspace crates whose versions are not yet on crates.io.
#
# Usage (from repository root):
#   CARGO_REGISTRY_TOKEN=... bash scripts/publish-changed.sh
#   CARGO_REGISTRY_TOKEN=... bash scripts/publish-changed.sh --dry-run
#   CARGO_REGISTRY_TOKEN=... bash scripts/publish-changed.sh --confirm-ffi
#
# Tag-triggered releases (e.g. v0.9.0) should set CARGO_REGISTRY_TOKEN in CI.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

DRY_RUN=0
CONFIRM_FFI=0
for arg in "$@"; do
  case "$arg" in
    --dry-run) DRY_RUN=1 ;;
    --confirm-ffi) CONFIRM_FFI=1 ;;
    -h | --help)
      sed -n '2,12p' "$0"
      exit 0
      ;;
    *)
      echo "error: unknown argument: $arg" >&2
      exit 2
      ;;
  esac
done

if [[ -z "${CARGO_REGISTRY_TOKEN:-}" ]]; then
  echo "error: CARGO_REGISTRY_TOKEN must be set to publish to crates.io" >&2
  exit 1
fi

export CARGO_REGISTRY_TOKEN

is_publishable() {
  local manifest="$1"
  local publish
  publish="$(cargo metadata --no-deps --format-version 1 --manifest-path "$manifest" \
    | python3 -c 'import json,sys; m=json.load(sys.stdin)["packages"][0]; print(str(m.get("publish", True)).lower())')"
  [[ "$publish" != "false" ]]
}

crate_name() {
  cargo metadata --no-deps --format-version 1 --manifest-path "$1" \
    | python3 -c 'import json,sys; print(json.load(sys.stdin)["packages"][0]["name"])'
}

crate_version() {
  cargo metadata --no-deps --format-version 1 --manifest-path "$1" \
    | python3 -c 'import json,sys; print(json.load(sys.stdin)["packages"][0]["version"])'
}

already_on_crates_io() {
  local name="$1"
  local version="$2"
  if cargo search "$name" --limit 1 2>/dev/null | grep -q "^$name = \"$version\""; then
    return 0
  fi
  return 1
}

publish_manifest() {
  local manifest="$1"
  local name version
  name="$(crate_name "$manifest")"
  version="$(crate_version "$manifest")"

  if ! is_publishable "$manifest"; then
    echo "skip $name ($manifest): publish = false"
    return 0
  fi

  if already_on_crates_io "$name" "$version"; then
    echo "skip $name $version: already on crates.io"
    return 0
  fi

  echo "publish $name $version from $manifest"
  if [[ "$DRY_RUN" -eq 1 ]]; then
    cargo publish --manifest-path "$manifest" --dry-run --allow-dirty
  else
    cargo publish --manifest-path "$manifest" --allow-dirty
  fi
}

# Primary crate (tag releases such as v0.9.0).
publish_manifest "$ROOT/Cargo.toml"

if [[ "$CONFIRM_FFI" -eq 1 ]]; then
  ffi_manifest="$ROOT/siderust-ffi/Cargo.toml"
  if ! is_publishable "$ffi_manifest"; then
    echo "error: --confirm-ffi requested but siderust-ffi has publish = false" >&2
    echo "Set publish = true in siderust-ffi/Cargo.toml before releasing the FFI crate." >&2
    exit 1
  fi
  publish_manifest "$ffi_manifest"
fi

echo "publish-changed: done"
