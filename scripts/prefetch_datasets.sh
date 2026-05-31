#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
prefetch_datasets.sh

Downloads Siderust build-time datasets and opt-in runtime SPKs into
$SIDERUST_DATASETS_DIR.

Usage:
  ./scripts/prefetch_datasets.sh --minimal
  ./scripts/prefetch_datasets.sh --all
  ./scripts/prefetch_datasets.sh --de440
  ./scripts/prefetch_datasets.sh --de441
  ./scripts/prefetch_datasets.sh --planet-centers
  ./scripts/prefetch_datasets.sh --mar099
  ./scripts/prefetch_datasets.sh --jup365
  ./scripts/prefetch_datasets.sh --sat441l
  ./scripts/prefetch_datasets.sh --ura184
  ./scripts/prefetch_datasets.sh --nep097

Environment:
  SIDERUST_DATASETS_DIR  Destination directory (recommended).
                        If unset, defaults to ./.siderust_datasets
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATASETS_DIR="${SIDERUST_DATASETS_DIR:-"$ROOT_DIR/.siderust_datasets"}"

mkdir -p "$DATASETS_DIR"

downloader=""
if command -v curl >/dev/null 2>&1; then
  downloader="curl"
elif command -v wget >/dev/null 2>&1; then
  downloader="wget"
else
  echo "error: need curl or wget" >&2
  exit 1
fi

download() {
  local url="$1"
  local out="$2"
  local tmp="${out}.download"
  mkdir -p "$(dirname "$out")"
  if [[ -s "$out" ]]; then
    return 0
  fi
  echo "downloading: $url"
  rm -f "$tmp"
  if [[ "$downloader" == "curl" ]]; then
    curl -fSL --retry 3 --retry-delay 2 -o "$tmp" "$url" || {
      rm -f "$tmp"
      return 1
    }
  else
    wget -q --tries=3 --waitretry=2 -O "$tmp" "$url" || {
      rm -f "$tmp"
      return 1
    }
  fi
  mv "$tmp" "$out"
}

prefetch_vsop87() {
  local dest="$DATASETS_DIR/vsop87_dataset"
  mkdir -p "$dest"

  local base="https://ftp.imcce.fr/pub/ephem/planets/vsop87/"
  local listing
  if [[ "$downloader" == "curl" ]]; then
    listing="$(curl -fsSL "$base")"
  else
    listing="$(wget -qO- "$base")"
  fi

  mapfile -t files < <(printf '%s' "$listing" | grep -Eo 'VSOP87[A-Z]\.[A-Za-z]{3}' | sort -u)
  if [[ "${#files[@]}" -eq 0 ]]; then
    echo "error: could not parse VSOP87 directory listing" >&2
    exit 1
  fi

  for f in "${files[@]}"; do
    download "${base}${f}" "$dest/$f"
  done
}

prefetch_elp2000() {
  local dest="$DATASETS_DIR/elp2000_dataset"
  mkdir -p "$dest"
  local base="https://cdsarc.cds.unistra.fr/ftp/VI/79"
  for n in $(seq 1 36); do
    local name="ELP${n}"
    download "${base}/${name}" "$dest/$name"
  done
}

prefetch_iers() {
  local dest="$DATASETS_DIR/iers_dataset"
  mkdir -p "$dest"

  local out="$dest/finals2000A.all"
  if [[ -s "$out" ]]; then
    return 0
  fi

  local primary="https://datacenter.iers.org/products/eop/rapid/standard/finals2000A.all"
  local fallback="https://maia.usno.navy.mil/ser7/finals2000A.all"

  if download "$primary" "$out"; then
    return 0
  fi
  rm -f "$out"
  download "$fallback" "$out"
}

prefetch_de440() {
  local dest="$DATASETS_DIR/de440_dataset"
  mkdir -p "$dest"
  local out="$dest/de440.bsp"
  local url="https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440.bsp"
  download "$url" "$out"
}

prefetch_de441() {
  local dest="$DATASETS_DIR/de441_dataset"
  mkdir -p "$dest"
  local out="$dest/de441_part-2.bsp"
  local url="https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de441_part-2.bsp"
  download "$url" "$out"
}

prefetch_mar099() {
  download \
    "https://ssd.jpl.nasa.gov/ftp/eph/satellites/bsp/mar099.bsp" \
    "$DATASETS_DIR/mar099.bsp"
}

prefetch_jup365() {
  download \
    "https://ssd.jpl.nasa.gov/ftp/eph/satellites/bsp/jup365.bsp" \
    "$DATASETS_DIR/jup365.bsp"
}

prefetch_sat441l() {
  download \
    "https://ssd.jpl.nasa.gov/ftp/eph/satellites/bsp/sat441l.bsp" \
    "$DATASETS_DIR/sat441l.bsp"
}

prefetch_ura184() {
  download \
    "https://ssd.jpl.nasa.gov/ftp/eph/satellites/bsp/ura184.bsp" \
    "$DATASETS_DIR/ura184.bsp"
}

prefetch_nep097() {
  download \
    "https://ssd.jpl.nasa.gov/ftp/eph/satellites/bsp/nep097.bsp" \
    "$DATASETS_DIR/nep097.bsp"
}

want_minimal=0
want_all=0
want_de440=0
want_de441=0
want_planet_centers=0
want_mar099=0
want_jup365=0
want_sat441l=0
want_ura184=0
want_nep097=0

for arg in "$@"; do
  case "$arg" in
    --minimal) want_minimal=1 ;;
    --all) want_all=1 ;;
    --de440) want_de440=1 ;;
    --de441) want_de441=1 ;;
    --planet-centers) want_planet_centers=1 ;;
    --mar099) want_mar099=1 ;;
    --jup365) want_jup365=1 ;;
    --sat441l) want_sat441l=1 ;;
    --ura184) want_ura184=1 ;;
    --nep097) want_nep097=1 ;;
    *)
      echo "error: unknown arg: $arg" >&2
      usage >&2
      exit 2
      ;;
  esac
done

if [[ "$#" -eq 0 ]]; then
  want_minimal=1
fi

if [[ "$want_all" -eq 1 ]]; then
  want_minimal=1
  want_de440=1
  want_de441=1
fi

if [[ "$want_planet_centers" -eq 1 ]]; then
  want_de440=1
  want_mar099=1
  want_jup365=1
  want_sat441l=1
  want_ura184=1
  want_nep097=1
fi

if [[ "$want_minimal" -eq 1 ]]; then
  prefetch_vsop87
  prefetch_elp2000
fi

if [[ "$want_de440" -eq 1 ]]; then
  prefetch_de440
fi

if [[ "$want_de441" -eq 1 ]]; then
  prefetch_de441
fi

if [[ "$want_mar099" -eq 1 ]]; then
  prefetch_mar099
fi

if [[ "$want_jup365" -eq 1 ]]; then
  prefetch_jup365
fi

if [[ "$want_sat441l" -eq 1 ]]; then
  prefetch_sat441l
fi

if [[ "$want_ura184" -eq 1 ]]; then
  prefetch_ura184
fi

if [[ "$want_nep097" -eq 1 ]]; then
  prefetch_nep097
fi

echo
echo "done. datasets dir: $DATASETS_DIR"
echo "tip: export SIDERUST_DATASETS_DIR=\"$DATASETS_DIR\""
