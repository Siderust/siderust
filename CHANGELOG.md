# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Add typed tabulated Cartesian ephemerides with cubic Hermite interpolation,
  multi-body lookup, coverage queries, and TDB CCSDS OEM loading.
- Add a reusable POD helper for one-way inter-satellite range with a light-time
  correction.
- Add a LISA POD example using the tabulated ephemeris provider and range helper.
- Add typed HEALPix primitives for reusable full-sky maps: validated `Nside`,
  `HealpixIndex`, `HealpixOrdering`, `HealpixGrid`, and frame-typed
  `HealpixMap<F, T>`.
- Add frame-neutral HEALPix pixelization from typed Cartesian `Direction<F>`
  values, with a documented private RING kernel using HEALPix `theta`, `phi`,
  and `z = cos(theta)` terminology.
- Add integrated-starlight map APIs for Galactic HEALPix surface-brightness maps,
  including catalogue records, apparent-magnitude validation, S10 photometry,
  provenance metadata, deterministic CSV serialization, and scientific
  validators for generated map assets.

### Changed

- Update `affn` to 0.8, `optica` to 0.2, `parquet` to 59, and the FFI dependency
  to the crates.io `tempoch-ffi` 0.6.6 release.
- Make stable-toolchain FFI builds reuse the checked-in C header while retaining
  header generation on nightly toolchains.
- Document the new HEALPix and starlight modules with scientific scope,
  technical scope, and primary references, following the workspace
  `missing_docs = "deny"` policy.

## [0.10.1] - 2026-06-20

### Added

* Added fixed ICRS ↔ Galactic direction-frame transforms.
* Added composed EquatorialMeanJ2000 ↔ Galactic and EclipticMeanJ2000 ↔ Galactic direction-frame transforms.
* Added spherical `Direction` frame-transform coverage for Galactic coordinates through the existing Cartesian round-trip blanket implementation.
* Added Galactic frame regression tests covering round-trip stability, SOFA reference matrix values, and the IAU Galactic north-pole direction.

### Changed

* Isolated the IAU/Hipparcos Galactic realization in `coordinates::transform::frames::galactic`, keeping the J2000 `bias` module focused on frame-bias and obliquity helpers.

### Fixed

* Avoided requiring nightly/cbindgen header generation for regular `siderust-ffi` builds when a checked-in FFI header is available.

## [0.10.0] - 2026-06-07

### Added

- Internal solar daily predictor for solar altitude threshold events.
- Internal Chebyshev-first generic crossing engine using `cheby 0.4`.
- Internal local scan+Brent fallback baseline.
- `bench-internals` feature for Criterion baseline comparisons.

### Changed

- Public altitude period API standardized on [`altitude_ranges`], [`above_threshold`], and [`below_threshold`], plus event-level [`crossings`] and [`culminations`].
- [`SearchOpts`] follows Option A and contains only `time_tolerance`.
- [`AltitudeProvider`] now only represents single-point altitude evaluation.
- Solar altitude events use an internal daily predictor with precise Brent validation and local fallback.
- Lunar altitude events use `MoonAltitudeContext` + Chebyshev-first crossing discovery.
- FFI [`SiderustSearchOpts`] mirrors public [`SearchOpts`] (time tolerance only).

### Removed

- `unstable-event-search` feature and experimental algorithm/tuning hooks.
- FFI `_v2` altitude/crossing functions and v2 tuning structs.
- Public algorithm/tuning structs and public scan-step override.
- `AltitudeQuery`, `AltitudePeriodsProvider`, `altitude_periods`, `compute_altitude_periods`, and target-specific legacy period wrappers.
- Local Chebyshev root code superseded by `cheby`.

## [0.9.1] - 2026-06-06

### Changed

- Speed up solar altitude, daylight, twilight, and night-period searches by using directed threshold crossings that avoid extra midpoint classification probes.
- Route default apparent-altitude threshold searches through provider-optimized paths when custom scan options are not requested.
- Update solar altitude benchmarks to cover the optimized night-event search path.

### Added

- Add directed above-threshold and in-range interval search helpers, with regression coverage against the existing generic scan implementation.

## [0.9.0] - 2026-06-01

### Changed

- Siderust is now a clean orchestration crate: canonical scientific datasets live in the [`siderust-archive`](https://crates.io/crates/siderust-archive)
  dependency; time-scale conversion, leap seconds, ΔT, and EOP ownership belong to [`tempoch`](https://crates.io/crates/tempoch).
- Removed `siderust::datasets`, `siderust::archive`, and the former `src/data/` tree. Dataset catalog and optional runtime acquisition are
  provided by `siderust-archive` (`siderust_archive::datasets`, `siderust_archive::runtime`).
- Dropped legacy Cargo features `archive-data`, `embedded-data`, `generated-tables`, and `external-data`. The `runtime-data` feature now
  enables `siderust-archive/fetch` for on-demand kernel download.
- The `lagrange-centers` feature is consumer-only: it enables archive-backed Sun–Earth Lagrange SCK kernels via `siderust-archive/lagrange`. Generation
  lives in `siderust-archive/tools/generate-lagrange-cheby`; Siderust no longer ships a local generator binary or build-time data mutation.
- Earth-rotation helpers delegate UTC/EOP indexing to `tempoch`'s active bundle. `try_jd_utc_from_tt` and `try_gmst_with_eop` fail when no runtime EOP
  is loaded; `gmst_default` uses the ΔT model only. `jd_utc_from_tt` remains as a ΔT-based compatibility wrapper.
- VSOP87/ELP coefficient tables are consumed from `siderust-archive` static snapshots; Siderust does not regenerate them at build time.
- Generic grid/spectrum/math ownership delegated to `optica`; orbital mechanics integration/STM/covariance delegated to `principia`.
- Reorganized internal namespaces: former `calculus::*` responsibilities moved into `event::*`, `ephemeris::*`, and `astro::*`; `spectra` renamed to  `photometry`.
- Force-model and propagation APIs use typed `qtty` quantities (`GravitationalParameter`, `Second`, …) instead of raw scalars.
- Default features are now `serde` only; POD lives behind the explicit `pod` feature family (`pod-parquet`, `pod-doris`).
- Extend runtime ephemeris loading to distinguish SPK data types and resolve Type 2/Type 3 segments consistently.
- Update runtime ephemeris internals to use stacked SPK segments for Sun, EMB, Moon, and optional Earth chains.
- Document required runtime satellite SPKs for exact Mars-through-Neptune planet-center chains.
- Clarify ERA implementation against the ERFA/SOFA `era00` split-date algorithm and avoid fractional-day double counting.

### Added

- Dependency on `siderust-archive` for embedded ephemeris tables, nutation, gravity, atmosphere, Pluto, and optional runtime fetch.
- `spice` feature with `siderust::spice` / `formats::spice` for SPICE kernels and `SpiceContext`.
- `pod` feature folding the former `siderust-pod` crate into `siderust::pod` (force models, propagation, estimation, QC, products, I/O).
- Spacecraft dynamics under `astro::dynamics`, SGP4/TLE/OMM support, and mission/format modules (CCSDS, IGS, ILRS, RINEX, VLBI).
- Lagrange-center support via archive-backed Chebyshev kernels.
- FFI dynamics bindings and tests.
- Add SPK Type 3 segment support, including velocity coefficient handling and Type 3-aware position/velocity evaluation.
- Add runtime SPK kernel management via `SpkKernelSet`, enabling dynamic loading and target-center geometric state resolution.
- Support multiple SPK segments per target/center pair through segment stacks, improving handling of overlapping kernel coverage.
- Add runtime planet-center SPK support for Mars, Jupiter, Saturn, Uranus, and Neptune satellite offset kernels.
- Add `de440` feature flag and register `de441`, `spk_perf`, and `sidereal` Criterion benchmarks.
- Add SPK parsing/evaluation/runtime benchmarks and dedicated sidereal GMST/GAST benchmarks.

### Removed

- Local canonical dataset regeneration (`regen-data`, `build.rs` archive generation, `scripts/generate-lagrange-cheby.rs`, and related build deps).
