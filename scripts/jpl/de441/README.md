# DE441 Build Pipeline

This module provides build-time integration of JPL's DE441 ephemeris using
NAIF's `de441_part-2.bsp` kernel.

## Overview

DE441 is a long-span numerical integration ephemeris from NASA JPL. This
pipeline mirrors the DE440 workflow:

1. Acquire the DE441 SPK (from Git LFS backup or NAIF fallback)
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
| `dataset/de441_part-2.bsp` | Git LFS backup of DE441 part-2 kernel |

## Data Source

Preferred: Git LFS backup at `dataset/de441_part-2.bsp`.
Fallback: NAIF HTTPS download.

```text
https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de441_part-2.bsp
```

## Build Process

When building with `--features de441`, the build script:

1. Reuses `$OUT_DIR/de441_dataset/de441_part-2.bsp` when valid
2. Tries `scripts/de441/dataset/de441_part-2.bsp` (Git LFS)
3. Downloads from NAIF as fallback
4. Parses DAF summaries and extracts Sun/EMB/Moon Type-2 segments
5. Writes binary blobs (`de441_{sun,emb,moon}.bin`)
6. Generates `de441_data.rs` with 8-byte aligned `include_bytes!` accessors

## Git LFS

The kernel is tracked by the repository-wide dataset rules in `.gitattributes`:

```gitattributes
scripts/**/dataset/** filter=lfs diff=lfs merge=lfs -text
```

## Manual Download

```bash
curl -fSL --max-time 5400 \
  https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de441_part-2.bsp \
  -o scripts/de441/dataset/de441_part-2.bsp
```

Or with `wget`:

```bash
wget --timeout=5400 \
  https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de441_part-2.bsp \
  -O scripts/de441/dataset/de441_part-2.bsp
```

## References

- NAIF generic planetary kernels: <https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/>
- SPICE toolkit docs: <https://naif.jpl.nasa.gov/naif/toolkit.html>
- DAF required reading: <https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/daf.html>
- SPK required reading: <https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html>
