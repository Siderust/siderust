# DE441 Build Pipeline

This module provides build-time integration of JPL's DE441 ephemeris using
NAIF's `de441_part-2.bsp` kernel.

## Overview

DE441 is a long-span numerical integration ephemeris from NASA JPL. This
pipeline mirrors the DE440 workflow:

1. Acquire the DE441 SPK (download or reuse a cached copy)
2. Parse the DAF (Double Precision Array File) container
3. Extract SPK Type 2 Chebyshev segments for Sun, EMB, and Moon
4. Generate statically-embedded Rust data and accessors

## Files

| File | Purpose |
|------|---------|
| `mod.rs` | Thin configuration wrapper calling `../jpl/pipeline.rs` |
| `../jpl/pipeline.rs` | Shared build pipeline (download, parse, codegen) |
| `../jpl/daf.rs` | Shared DAF parser |
| `../jpl/spk.rs` | Shared SPK Type 2 segment reader |

## Data Source

NAIF HTTPS download:

```text
https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de441_part-2.bsp
```

## Build Process

When building with `--features de441`, the build script:

1. Reuses `de441_part-2.bsp` if already present in the dataset directory (see `doc/datasets.md`)
2. Downloads from NAIF if missing
3. Parses DAF summaries and extracts Sun/EMB/Moon Type-2 segments
4. Writes binary blobs (`de441_{sun,emb,moon}.bin`)
5. Generates `de441_data.rs` with 8-byte aligned `include_bytes!` accessors

## Skipping the Download (Typecheck Only)

To make local `cargo check --all-features` fast without downloading the BSP,
set `SIDERUST_JPL_STUB` to `de441`:

```bash
SIDERUST_JPL_STUB=de441 cargo check --all-features
```

For tests, you can use the same variable:

```bash
SIDERUST_JPL_STUB=de441 cargo test --features de441
```

To stub all JPL datasets (DE440 + DE441), use:

```bash
SIDERUST_JPL_STUB=all cargo check --all-features
```

This generates a small stub `de441_data.rs` so the crate compiles, but any
direct use of low-level DE441 JPL data will panic.

For the public ephemeris backend (`calculus::ephemeris::De441Ephemeris`), the
same setting enables a mock backend that falls back to VSOP87/ELP2000 so tests
can execute without downloading DE441.

## Manual Download

```bash
export SIDERUST_DATASETS_DIR="$PWD/.siderust_datasets"
mkdir -p "$SIDERUST_DATASETS_DIR/de441_dataset"
curl -fSL --max-time 5400 \
  https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de441_part-2.bsp \
  -o "$SIDERUST_DATASETS_DIR/de441_dataset/de441_part-2.bsp"
```

Or with `wget`:

```bash
export SIDERUST_DATASETS_DIR="$PWD/.siderust_datasets"
mkdir -p "$SIDERUST_DATASETS_DIR/de441_dataset"
wget --timeout=5400 \
  https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de441_part-2.bsp \
  -O "$SIDERUST_DATASETS_DIR/de441_dataset/de441_part-2.bsp"
```

## References

- NAIF generic planetary kernels: <https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/>
- SPICE toolkit docs: <https://naif.jpl.nasa.gov/naif/toolkit.html>
- DAF required reading: <https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/daf.html>
- SPK required reading: <https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html>
