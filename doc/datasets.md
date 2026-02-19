# Datasets, Codegen, and JPL Backends (offline-by-default)

Siderust relies on several datasets that are embedded at compile time.
This document explains what is embedded, how the offline-first strategy works,
and how to refresh or extend the pre-generated data.

## Offline-by-default: committed generated tables

The three core datasets are **pre-generated and committed** under
`src/generated/`.  Normal builds — including docs.rs — do **not** download
anything:

| File | Source dataset | Update frequency |
|------|---------------|-----------------|
| `src/generated/vsop87a.rs` | VSOP87A (planetary theory) | ~stable |
| `src/generated/vsop87e.rs` | VSOP87E (planetary theory) | ~stable |
| `src/generated/elp_data.rs` | ELP2000-82B (lunar theory) | ~stable |
| `src/generated/iers_eop_data.rs` | IERS finals2000A.all (EOP) | weekly |

The library includes these files directly:
```rust
include!(concat!(env!("CARGO_MANIFEST_DIR"), "/src/generated/vsop87a.rs"));
```

`src/generated/datasets.lock.json` records the generation timestamp and
SHA-256 checksums of each file.

## Optional JPL DE4xx datasets (feature-gated)

When `de440` and/or `de441` features are enabled, the build scripts **still**
download large NAIF `.bsp` kernels (these cannot be committed due to size).

Typical sizes:
- DE440: ~120 MB
- DE441 part-2: ~1.65 GB

## How to regenerate locally

Use the helper script to download source datasets and regenerate all committed
Rust tables in one step:

```bash
# From the siderust crate root:
./scripts/update_generated_tables.sh
```

This does three things:
1. Downloads the source datasets via `scripts/prefetch_datasets.sh --minimal`.
2. Runs `SIDERUST_REGEN=1 cargo build` which overwrites `src/generated/`.
3. Writes an updated `src/generated/datasets.lock.json`.

Review with `git diff src/generated/`, then commit:

```bash
git add src/generated/
git commit -m "chore: refresh generated dataset tables"
```

Alternatively, you can drive the regen directly without the script:

```bash
export SIDERUST_DATASETS_DIR="$PWD/.siderust_datasets"
./scripts/prefetch_datasets.sh --minimal
SIDERUST_REGEN=1 cargo build
```

## How `SIDERUST_REGEN=1` works

`build.rs` checks the `SIDERUST_REGEN` environment variable at compile time:

- **Not set (default):** The build script is a no-op for VSOP87/ELP2000/IERS.
  The library reads the committed files from `src/generated/`.
- **`SIDERUST_REGEN=1`:** The build script downloads (or reuses cached) source
  datasets, parses them, and overwrites `src/generated/*.rs`.

## Automated CI refresh

A GitHub Actions workflow (`.github/workflows/update-datasets.yml`) runs every
Monday and on `workflow_dispatch`.  If any generated file changes (typically
`iers_eop_data.rs` is updated weekly with new IERS data), it opens a pull
request automatically.

## Where downloads go (caching for SIDERUST_REGEN and JPL)

By default, downloads are stored under Cargo's `OUT_DIR` (wiped by
`cargo clean`).

To reuse datasets across clean builds, set a persistent cache directory:

```bash
export SIDERUST_DATASETS_DIR="$PWD/.siderust_datasets"
```

## Manual download (direct URLs)

If you prefer downloading by hand, place files in the cache directory:

- VSOP87 → `$SIDERUST_DATASETS_DIR/vsop87_dataset/`
  - Source: `https://ftp.imcce.fr/pub/ephem/planets/vsop87/`
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
