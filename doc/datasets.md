# Datasets and runtime acquisition

Canonical scientific datasets (VSOP87, ELP2000, nutation, gravity, Lagrange
SCK kernels, optional JPL DE4xx) live in the
[`siderust-archive`](https://crates.io/crates/siderust-archive) crate. Siderust
consumes archive-backed static snapshots at compile time; it does not download,
regenerate, or mutate datasets locally.

IERS leap-second, ΔT, and EOP data are owned by [`tempoch`](https://crates.io/crates/tempoch)
and loaded at runtime when needed. See `tempoch/README.md` for the time-data
bundle contract.

Import `siderust_archive::datasets` and `siderust_archive::runtime` directly for
catalog and optional fetch helpers. Enable Siderust's `runtime-data` feature for
`DatasetManager`-style ephemeris download.

## Regeneration (archive workspace)

Dataset maintenance and generators belong in the
[`siderust-archive`](https://github.com/Siderust/archive) repository, not in
this crate. To regenerate Sun–Earth Lagrange SCK kernels:

```bash
cd path/to/archive
cargo run -p generate-lagrange-cheby -- \
    --source vsop87 --out src/lagrange/vsop87
```

For prefetch and validation scripts, see `archive/README.md` in the archive
checkout.

## Legacy note

Former `siderust/scripts/`, `src/generated/`, `src/data/`, and local `build.rs`
dataset generation were removed. This crate is orchestration-only for canonical
data.

## Runtime ephemeris cache

When using `runtime-data`, downloaded kernels are cached under a configurable
directory (see `siderust_archive::runtime` and example `12_runtime_ephemeris`).

For offline JPL testing without large downloads:

```bash
SIDERUST_JPL_STUB=all cargo test --all-features
```

See `README.md` for the recommended local override file pattern
(`.cargo/config.local.toml`).
