# Datasets, Codegen, and JPL Backends (no Git LFS)

Siderust relies on several datasets that are generated/embedded at build time.
This document explains what is embedded, when downloads happen, and how to do
fast/offline development loops.

This repository previously stored large datasets via **Git LFS**, but GitHub
applies LFS bandwidth limits. The project does not strictly need LFS because
the build pipeline can fetch the raw datasets on demand. As a result, datasets
are **not committed** anymore.

## Always-generated datasets (build-time)

- **VSOP87** (planetary analytical series)
- **ELP2000-82B** (lunar analytical series)
- **IERS EOP** (`finals2000A.all`) for Earth-orientation parameters

These are generated into `OUT_DIR` by `build.rs` and included by the crate at
compile time.

## Optional JPL DE4xx datasets (feature-gated)

When `de440` and/or `de441` features are enabled, the build scripts may download
large NAIF `.bsp` kernels and extract the segments required by the library.

Typical sizes:
- DE440: ~120 MB
- DE441 part-2: ~1.65 GB

## Where downloads go (caching)

By default, downloads are stored under Cargo’s `OUT_DIR` (so they may be wiped
by `cargo clean`).

To reuse datasets across clean builds and to support offline builds after a
one-time download, set a persistent cache directory:

```bash
export SIDERUST_DATASETS_DIR="$PWD/.siderust_datasets"
```

The build will then store datasets under that directory (e.g.
`$SIDERUST_DATASETS_DIR/de440_dataset/de440.bsp`).

## Manual prefetch (recommended)

Use the helper script to download datasets into `SIDERUST_DATASETS_DIR`:

```bash
export SIDERUST_DATASETS_DIR="$PWD/.siderust_datasets"
./scripts/prefetch_datasets.sh --minimal
```

Optional (large):

```bash
./scripts/prefetch_datasets.sh --de440
./scripts/prefetch_datasets.sh --de441
```

## Manual download (direct URLs)

If you prefer downloading by hand, place files here:

- VSOP87 → `$SIDERUST_DATASETS_DIR/vsop87_dataset/`
  - Source listing: `https://ftp.imcce.fr/pub/ephem/planets/vsop87/`
- ELP2000 → `$SIDERUST_DATASETS_DIR/elp2000_dataset/` (files `ELP1`…`ELP36`)
  - Base URL: `https://cdsarc.cds.unistra.fr/ftp/VI/79/`
- IERS EOP → `$SIDERUST_DATASETS_DIR/iers_dataset/finals2000A.all`
  - Primary: `https://datacenter.iers.org/products/eop/rapid/standard/finals2000A.all`
  - Mirror:  `https://maia.usno.navy.mil/ser7/finals2000A.all`
- DE440 → `$SIDERUST_DATASETS_DIR/de440_dataset/de440.bsp`
  - `https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440.bsp`
- DE441 → `$SIDERUST_DATASETS_DIR/de441_dataset/de441_part-2.bsp`
  - `https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de441_part-2.bsp`

## Stubbing JPL datasets (fast/offline)

To compile with JPL features enabled while skipping downloads/codegen:

```bash
SIDERUST_JPL_STUB=all cargo test --all-features
```

Supported values:
- `de440` (stub DE440 only)
- `de441` (stub/mock DE441 only)
- `de440,de441` or `all` (stub both)

See `README.md` for the recommended local override file pattern
(`.cargo/config.local.toml`).
