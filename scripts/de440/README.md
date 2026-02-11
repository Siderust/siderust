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
| `mod.rs` | Build pipeline entry point, orchestrates the process |
| `fetch.rs` | Downloads `de440.bsp` from NAIF or copies from Git LFS backup |
| `daf.rs` | DAF file format parser (NAIF Double Precision Array File) |
| `spk.rs` | SPK Type 2 segment reader (Chebyshev polynomials) |
| `dataset/de440.bsp` | **Git LFS backup** of the 115 MB ephemeris file |

## Data Source

**Primary**: Git LFS backup in `dataset/de440.bsp` (115 MB)  
**Fallback**: NAIF HTTPS download

```
https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440.bsp
```

## Build Process

When building with `--features de440`, the build script:

1. **Checks cache**: If `$OUT_DIR/de440_dataset/de440.bsp` exists with valid size (>100 MB), uses it
2. **Tries Git LFS**: Copies from `scripts/de440/dataset/de440.bsp` if present
3. **Downloads fallback**: Uses `curl`/`wget` to download from NAIF (900s timeout)
4. **Parses DAF**: Reads file records, summary records, name records
5. **Extracts segments**: Type 2 Chebyshev segments for Sun, Earth-Moon barycenter, Moon
6. **Generates code**: Creates `de440_data.rs` with embedded coefficients aligned to 8 bytes

## Git LFS

The `de440.bsp` file is tracked via Git LFS to avoid:
- Slow CI builds (115 MB HTTPS download on every clean build)
- Network failures interrupting builds
- NAIF server rate limiting

Configuration in `.gitattributes`:
```gitattributes
scripts/**/dataset/** filter=lfs diff=lfs merge=lfs -text
```

This pattern automatically tracks all files in `scripts/{theory}/dataset/` directories, including:
- `scripts/vsop87/dataset/` (VSOP87 theory data, ~27 MB)
- `scripts/elp2000/dataset/` (ELP2000-82B lunar theory, ~2.5 MB)
- `scripts/de440/dataset/` (DE440 ephemeris, ~115 MB)

## Manual Download

If you need to manually download the file:

```bash
curl -fSL --max-time 900 \
  https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440.bsp \
  -o scripts/de440/dataset/de440.bsp
```

Or with `wget`:

```bash
wget --timeout=900 \
  https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440.bsp \
  -O scripts/de440/dataset/de440.bsp
```

## References

- **DE440**: [JPL Planetary and Lunar Ephemerides](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/aareadme_de440_de441.txt)
- **NAIF SPICE**: [SPICE Toolkit Documentation](https://naif.jpl.nasa.gov/naif/toolkit.html)
- **DAF Format**: [Required Reading - DAF](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/daf.html)
- **SPK Format**: [Required Reading - SPK](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html)
