# Datasets, Codegen, and JPL Backends (offline-by-default)

Dataset cataloging, runtime download/cache managers, and IERS time-data providers
live in the [`siderust-archive`](../../archive/crates/siderust-archive) crate.
Import `siderust_archive::datasets` and `siderust_archive::runtime` directly;
`siderust` no longer re-exports those modules.

Siderust relies on several datasets that are embedded at compile time.
This document explains what is embedded, how the offline-first strategy works,
and how to refresh or extend the pre-generated data.

## Offline-by-default: committed generated tables

The core Siderust ephemeris datasets are **pre-generated and committed** under
`src/generated/`. Normal builds, including docs.rs, do **not** download
anything:

| File / asset | Embedded table | Official reference | Canonical upstream | Update frequency |
|------|---------------|--------------------|--------------------|-----------------|
| `src/generated/vsop87a.rs` | VSOP87A rectangular heliocentric series | Bretagnon, P., & Francou, G. (1988), *Astronomy and Astrophysics* **202**, 309, "Planetary theories in rectangular and spherical variables. VSOP87 solutions" | CDS VI/81 / IMCCE VSOP87 distribution | ~stable |
| `src/generated/vsop87e.rs` | VSOP87E barycentric spherical series | Bretagnon, P., & Francou, G. (1988), *Astronomy and Astrophysics* **202**, 309, "Planetary theories in rectangular and spherical variables. VSOP87 solutions" | CDS VI/81 / IMCCE VSOP87 distribution | ~stable |
| `src/generated/elp_data.rs` | ELP2000-82B lunar series (`ELP1`…`ELP36`) | Chapront-Touzé, M., & Chapront, J. (1983), *Astronomy and Astrophysics* **124**, 50, "The Lunar Ephemeris ELP 2000"; and Chapront-Touzé, M., & Chapront, J. (1988), *Astronomy and Astrophysics* **190**, 342, "ELP 2000-85: a semi-analytical lunar ephemeris adequate for historical times" | CDS VI/79 | ~stable |

Earth-orientation data is owned by `tempoch`; `siderust::astro::IersEop`
consumes `tempoch`'s bundled EOP tables instead of regenerating a second active
IERS table in this crate. The legacy `siderust::astro::iers_data` module is
kept for compatibility, but new code should use `tempoch`/`IersEop`.

Other bundled scientific tables that are shipped directly as data files:

| File / asset | Embedded table | Official reference | Canonical upstream |
|------|---------------|--------------------|--------------------|
| `data/o3trans.dat` | Ozone transmittance vs wavelength | Patat, F., et al. (2008), *A&A* **481**, 575, "An Atlas of the Sky Background Spectrum over Cerro Paranal" | Bundled from the NSB / `darknsb` lineage; canonical copy in `siderust::atmosphere::ozone` |
| `data/passbands/bessell1990/U.dat` | Johnson U passband | Bessell, M. S. (1990), *PASP* **102**, 1181, "UBVRI Passbands" | SVO Filter Profile Service `Generic/Bessell.U` |
| `data/passbands/bessell1990/B.dat` | Johnson B passband | Bessell, M. S. (1990), *PASP* **102**, 1181, "UBVRI Passbands" | SVO Filter Profile Service `Generic/Bessell.B` |
| `data/passbands/bessell1990/V.dat` | Johnson V passband | Bessell, M. S. (1990), *PASP* **102**, 1181, "UBVRI Passbands" | SVO Filter Profile Service `Generic/Bessell.V` |
| `data/passbands/bessell1990/R.dat` | Cousins R passband | Bessell, M. S. (1990), *PASP* **102**, 1181, "UBVRI Passbands" | SVO Filter Profile Service `Generic/Bessell.R` |
| `data/passbands/bessell1990/I.dat` | Cousins I passband | Bessell, M. S. (1990), *PASP* **102**, 1181, "UBVRI Passbands" | SVO Filter Profile Service `Generic/Bessell.I` |

For the bundled UTC-TAI, ΔT, and IERS EOP tables consumed through `tempoch`,
see the corresponding reference inventory in `tempoch/README.md`.

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

Use the archive helper script to prefetch source datasets (VSOP87, ELP2000, JPL
kernels). Table regeneration is owned by `siderust-archive`; see
[`archive/README.md`](../../archive/README.md).

```bash
# From the siderust crate root (monorepo layout):
../archive/scripts/prefetch_datasets.sh --minimal
```

For JPL kernels only:

```bash
export SIDERUST_DATASETS_DIR="$PWD/.siderust_datasets"
../archive/scripts/prefetch_datasets.sh --de440
```

To regenerate Sun-Earth Lagrange SCK kernels, use the archive tool:

```bash
cd ../archive
cargo run -p generate-lagrange-cheby -- \
    --source vsop87 --out src/lagrange/vsop87
```

Legacy note: `siderust/scripts/update_generated_tables.sh` and
`siderust/scripts/prefetch_datasets.sh` were removed; use the archive paths
above.

## How `SIDERUST_REGEN=1` works

`build.rs` checks the `SIDERUST_REGEN` environment variable at compile time:

- **Not set (default):** The build script is a no-op for VSOP87/ELP2000.
  The library reads the committed files from `src/generated/`.
- **`SIDERUST_REGEN=1`:** The build script downloads (or reuses cached) source
  datasets, parses them, and overwrites `src/generated/*.rs`.

## Automated CI refresh

A GitHub Actions workflow (`.github/workflows/update-datasets.yml`) runs on
`workflow_dispatch` for the active Siderust tables. EOP refreshes belong in
the `tempoch` data pipeline.

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
