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

## Runtime planet-center SPKs

Generic DE planetary kernels provide the outer planet-system barycenters. Exact
Mars through Neptune center chains need the corresponding satellite SPK offset
loaded together with a DE kernel at runtime. Cache and acquisition are handled
by `siderust_archive::runtime`; manual setups should keep the NAIF file identity
visible in run provenance.

| Center | Runtime SPK |
|---|---|
| Mars `499` | `mar099.bsp` |
| Jupiter `599` | `jup365.bsp` |
| Saturn `699` | `sat441l.bsp` |
| Uranus `799` | `ura184.bsp` |
| Neptune `899` | `nep097.bsp` |

These files are runtime data; they are not embedded into the crate.

## Offline / CI testing

`cargo test --all-features` does not download JPL BSP files. CI sets
`SIDERUST_JPL_STUB=all` for deterministic workspace checks.

For real-kernel checks, provide a local BSP:

```bash
SIDERUST_BSP_PATH=/path/to/de440.bsp cargo test --test test_jpl_real_backend
```

See `README.md` for `runtime-data` download and cache layout.
