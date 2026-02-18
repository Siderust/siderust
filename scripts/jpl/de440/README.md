# DE440 Build Pipeline

This module provides build-time integration of JPL's DE440 planetary ephemeris.

## Overview

DE440 is a high-precision numerical integration ephemeris from NASA JPL, providing positions and velocities for major solar system bodies. This build pipeline:

1. Downloads the DE440 Binary SPK (SPICE Kernel) file from NAIF
2. Parses the DAF (Double Precision Array File) format
3. Extracts Chebyshev polynomial segments
4. Generates statically-embedded Rust data structures

## Files

| File | Purpose |
|------|---------|
| `mod.rs` | Thin configuration wrapper calling `../jpl/pipeline.rs` |
| `../jpl/pipeline.rs` | Shared build pipeline (download, parse, codegen) |
| `../jpl/daf.rs` | Shared DAF file format parser (NAIF Double Precision Array File) |
| `../jpl/spk.rs` | Shared SPK Type 2 segment reader (Chebyshev polynomials) |

## Data Source

NAIF HTTPS download:

```
https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440.bsp
```

## Build Process

When building with `--features de440`, the build script:

1. **Checks cache**: Reuses `de440.bsp` if already present in the dataset directory (see `doc/datasets.md`)
2. **Downloads**: Uses `curl`/`wget` to download from NAIF (900s timeout) if missing
3. **Parses DAF**: Reads file records, summary records, name records
4. **Extracts segments**: Type 2 Chebyshev segments for Sun, Earth-Moon barycenter, Moon
5. **Generates code**: Creates `de440_data.rs` with embedded coefficients aligned to 8 bytes

## Skipping the Download (Typecheck Only)

To make local `cargo check --all-features` fast without downloading BSP files,
set `SIDERUST_JPL_STUB`:

```bash
SIDERUST_JPL_STUB=de440 cargo check --all-features
```

To stub all JPL datasets (DE440 + DE441), use:

```bash
SIDERUST_JPL_STUB=all cargo check --all-features
```

This generates a small stub `de440_data.rs` so the crate compiles, but any
attempt to use DE440 at runtime will panic.

## Manual Download

If you need to manually download the file:

```bash
export SIDERUST_DATASETS_DIR="$PWD/.siderust_datasets"
mkdir -p "$SIDERUST_DATASETS_DIR/de440_dataset"
curl -fSL --max-time 900 \
  https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440.bsp \
  -o "$SIDERUST_DATASETS_DIR/de440_dataset/de440.bsp"
```

Or with `wget`:

```bash
export SIDERUST_DATASETS_DIR="$PWD/.siderust_datasets"
mkdir -p "$SIDERUST_DATASETS_DIR/de440_dataset"
wget --timeout=900 \
  https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440.bsp \
  -O "$SIDERUST_DATASETS_DIR/de440_dataset/de440.bsp"
```

## References

- **DE440**: [JPL Planetary and Lunar Ephemerides](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/aareadme_de440_de441.txt)
- **NAIF SPICE**: [SPICE Toolkit Documentation](https://naif.jpl.nasa.gov/naif/toolkit.html)
- **DAF Format**: [Required Reading - DAF](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/daf.html)
- **SPK Format**: [Required Reading - SPK](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html)
