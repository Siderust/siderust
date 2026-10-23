#!/usr/bin/env bash
# Publish workspace crates whose versions are not yet on crates.io.
#
# Usage (from repository root):
#   CARGO_REGISTRY_TOKEN=... bash scripts/publish-changed.sh
#   bash scripts/publish-changed.sh --dry-run
#   CARGO_REGISTRY_TOKEN=... bash scripts/publish-changed.sh --confirm-ffi
#
# Tag-triggered releases (e.g. v0.15.0) should set CARGO_REGISTRY_TOKEN in CI.
# Dry runs do not require a registry token.
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

if [[ "$DRY_RUN" -eq 0 && -z "${CARGO_REGISTRY_TOKEN:-}" ]]; then
  echo "error: CARGO_REGISTRY_TOKEN must be set to publish to crates.io" >&2
  exit 1
fi

if [[ -n "${CARGO_REGISTRY_TOKEN:-}" ]]; then
  export CARGO_REGISTRY_TOKEN
fi

RELEASE_BIN="$ROOT/target/debug/siderust-release"
ensure_release_tool() {
  if [[ ! -x "$RELEASE_BIN" ]]; then
    cargo build -q -p siderust-release
  fi
}

# Prefer the crates.io version API over `cargo search`. Search only reports the
# newest crate and misses already-published versions during reruns.
already_on_crates_io() {
  local name="$1"
  local version="$2"
  ensure_release_tool
  if "$RELEASE_BIN" already-published "$name" "$version"; then
    return 0
  fi
  return 1
}

publish_manifest() {
  local manifest="$1"
  local name version status
  ensure_release_tool
  read -r name version status < <("$RELEASE_BIN" package-meta "$manifest")

  if [[ "$status" != "publishable" ]]; then
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

# Primary crate (tag releases such as v0.15.0).
publish_manifest "$ROOT/Cargo.toml"

if [[ "$CONFIRM_FFI" -eq 1 ]]; then
  ffi_manifest="$ROOT/siderust-ffi/Cargo.toml"
  ensure_release_tool
  read -r _name _version status < <("$RELEASE_BIN" package-meta "$ffi_manifest")
  if [[ "$status" != "publishable" ]]; then
    echo "error: --confirm-ffi requested but siderust-ffi has publish = false" >&2
    echo "Set publish = true in siderust-ffi/Cargo.toml before releasing the FFI crate." >&2
    exit 1
  fi
  publish_manifest "$ffi_manifest"
fi

echo "publish-changed: done"
