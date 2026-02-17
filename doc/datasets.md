# Datasets, Codegen, and JPL Backends

Siderust relies on several datasets that are generated/embedded at build time.
This document explains what is embedded, when downloads happen, and how to do
fast/offline development loops.

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
