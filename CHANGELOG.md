# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [unreleased]

### Changed

- Removed `siderust::datasets` and `siderust::archive`. Dataset catalog and
  runtime acquisition now come from the `siderust-archive` dependency
  (`siderust_archive::datasets`, `siderust_archive::runtime`).
- Dropped `archive-data`, `embedded-data`, `generated-tables`, and
  `external-data` Cargo features; `build.rs` no longer emits
  `archive_registry.rs`.
- `runtime-data` now enables `siderust-archive/runtime` instead of bundling
  `ureq`/`sha2` directly in this crate.
- Sun-Earth Lagrange SCK regeneration moved to
  `siderust-archive/tools/generate-lagrange-cheby` (archive workspace tool).
  The `generate-lagrange-cheby` binary and `scripts/` shims were removed from
  this crate.
- Dataset maintenance shell scripts now live under `siderust-archive/scripts/`
  only; the former `siderust/scripts/` redirect shims were deleted.

### Removed

- Deleted `src/datasets/` and `src/archive.rs`.
- Removed `scripts/generate-lagrange-cheby.rs`, `scripts/jpl/`, and other
  archive redirect shims from this crate.
- Removed the `archive/` git submodule. The reusable archive crate (now named
  `siderust-archive`) is published on crates.io; downstream consumers that
  want archive-backed data should clone https://github.com/Siderust/archive
  separately or set `SIDERUST_ARCHIVE_ROOT` to point at a checkout. The
  `build.rs` lagrange layout helper resolves the archive root via, in order:
  `SIDERUST_ARCHIVE_ROOT`, `./archive/`, `../archive/`.
- The `lagrange-centers` feature now requires the archive to be reachable at
  build time (same resolution rules). The build emits a `compile_error!` with
  an actionable message when the SCK kernels cannot be found, instead of the
  former opaque `include_bytes!` failure.

### Added

- Dependency on `siderust-archive` for dataset catalog and optional runtime
  ephemeris download.
- Added a reusable `siderust-archive-data` Rust crate inside the archive
  submodule (`archive/crates/siderust-archive-data`) so every repository in the
  Siderust family can share one data-access layer via a path dependency into
  the submodule. It provides: TOML manifest parsing (schema v1, top-level
  registry + per-family manifests), SHA-256 checksum verification, provenance,
  and a runtime IERS time-data downloader (`time::TimeDataManager`, behind the
  `fetch` feature) that downloads, verifies, caches, and records provenance for
  leap seconds, ΔT, and Earth Orientation Parameters. The crate carries no
  embedded data and declares an empty `[workspace]` so it is not absorbed into a
  consumer's cargo workspace. Added `archive/time/eop/` describing the IERS
  EOP/ΔT sources (fetched at runtime, not committed).
- Established the `archive/` git submodule as the canonical store for
  scientific datasets, kernels, manifests, generators, and validation
  reports. Added `.gitmodules`, `archive/README.md`, `archive/MANIFEST.toml`
  (top-level registry), per-family `manifest.toml` skeletons for VSOP87,
  IAU 2000A nutation, ELP2000-82B, Meeus 1998 Pluto, leap seconds, frames,
  body constants, and the generated Sun-Earth Lagrange Chebyshev kernels.
  All archive metadata is TOML; JSON is no longer used in the archive.
- Added `archive/schema/archive-manifest-v1.md` (manifest contract) and
  `archive/schema/sck-v1.md` (Siderust Chebyshev Kernel binary spec).
- Added directory skeletons for `archive/generators/lagrange/`,
  `archive/tools/validate/`, and `archive/tools/convert/` with README stubs
  describing planned standalone Cargo crates.
- Added `plans/archive-plan.md` documenting the 12-phase data-architecture
  migration.
- Added Cargo features `archive-data`, `embedded-data`, `generated-tables`,
  and `external-data` to opt in to the archive-backed dataset stack.
- Added `siderust::data::archive` with the build-time generated
  `ARCHIVE_ENTRIES` table and `lookup_family` helper. The registry is
  populated from `archive/MANIFEST.toml` when the `archive-data` feature
  is enabled and the submodule is checked out; otherwise it is empty.
- Added `siderust::formats::sck` reader for the Siderust Chebyshev Kernel
  v1 binary format used by the Lagrange archive payloads.
- Extended `build.rs` to detect the `archive/` submodule, emit an
  `archive_registry.rs` artifact into `OUT_DIR`, and watch the relevant
  manifest paths via `cargo:rerun-if-changed`. The build never regenerates
  scientific datasets and falls back to an empty registry when the
  submodule is missing.
- The Lagrange generator (`scripts/generate-lagrange-cheby.rs`) now
  accepts `--block-days` to control the Chebyshev block size, writes its
  outputs directly into `archive/lagrange/<source>/`, and emits an
  archive-compatible `manifest.toml` next to the `.sck` files.

- Added the `spice` feature and `siderust::spice` / `formats::spice` stack for SPICE text, DAF/SPK, FK, PCK, CK, SCLK, IK, kernel metadata, and high-level `SpiceContext`.
- Added `photometry` support for Johnson–Cousins UBVRI passbands and delegated generic sampled-spectrum/grid interpolation to `optica`.
- Added `pod` as an explicit opt-in feature, folding the former `siderust-pod` functionality into `siderust::pod` with force-model configuration, propagation, process noise, observations, estimation, QC, products, and SPICE/I/O adapters.
- Added spacecraft dynamics under `astro::dynamics`: typed orbital states, perturbation context, drag, SRP, geopotential, third-body, relativity, empirical accelerations, covariance/STM helpers, and propagator wrappers.
- Added SGP4/TLE/OMM support, CCSDS/IGS/ILRS/RINEX/VLBI format modules, and related examples/tests.
- Added aircraft identity/state and ICAO ISA atmosphere support under `bodies::aircraft`.
- Added Lagrange-center support, generated Chebyshev tooling/placeholders, and center-shift providers.
- Added FFI dynamics bindings and tests.

### Changed

- Reorganized major internal namespaces:
  - `calculus::*` functionality moved into `event::*`, `ephemeris::*`, and `astro::*`.
  - `archive` / `provenance` data moved under `data::*`.
  - generic grid/table functionality moved to `optica::grid`.
- Replaced the old local generic dynamics implementation with the extracted `principia` mechanics crate for integration, STM, covariance transport, gravity kernels, and local trajectory frames.
- Updated force-model and propagation APIs to use typed `qtty` quantities such as `GravitationalParameter` and `Second` instead of raw scalars.
- Renamed `siderust::spectra` to `siderust::photometry`; use `optica::spectrum` for generic spectrum infrastructure.
- Made `pod` non-default; default features now include `serde` only.
- Marked `siderust-ffi` as `publish = false`.
- Corrected and modernized examples, benches, docs, and generated bindings for the new module layout.

### Removed

- **Dissolved `src/data/` module** — the `data` namespace no longer exists.
  Sub-modules were redistributed:
  - `data::checksum` → `siderust::checksum` (top-level)
  - `data::archive` → `siderust::archive` (top-level)
  - `data::{DatasetId, DatasetMeta, DatasetError, …}` → `siderust::datasets`
  - `data::runtime` → `siderust::datasets::runtime`
  - `data::compiled::jpl` → `siderust::ephemeris::jpl` (co-located)
  - `data::provenance` / `data::DataSource` replaced by `optica::data::{Provenance, DataSource}` directly
  - Static data files (`o3trans.dat`, Bessell passbands) moved next to their only consumers
- Removed public `siderust::numeric`; event-search internals are now crate-private under `event::search`.
- Removed `siderust::calculus` as a public namespace after moving its responsibilities into `event`, `ephemeris`, and `astro`.
- Removed `siderust::archive` and old data/provenance layout in favor of `siderust::data`.
- Removed `siderust::tables` and the `tables` feature; use `optica::grid`.
- Removed `siderust::spectra`; use `siderust::photometry` and `optica::spectrum`.
- Removed the standalone `siderust-pod` crate and `pod-lisa` feature; POD now lives behind `siderust::pod`, and LISA-specific logic lives in examples.
- Removed the parallel `pod::dynamics` force-model wrapper layer in favor of canonical `AccelerationModel` / `CompositeModel`.
- Removed root-level `siderust::instruments`, `siderust::mission_geometry`, and
  `siderust::mission_context`; mission APIs now live under `siderust::mission`.
- **Removed `src/data/compiled/lagrange/`** — the 296 K-line generated Rust
  array is replaced by five compact SCK binary kernels committed to the
  `archive/` submodule (`archive/lagrange/vsop87/l{1..5}.sck`, each ~475 KB,
  2283 records, 1900–2100 coverage).  The runtime loader uses `include_bytes!`
  and a lazy parser instead of a giant static array.
- **Removed `src/data/compiled/{vsop87a,vsop87e,elp_data,nut00a_tables,pluto_tables}.rs`**
  from the `data::compiled` module; each file now lives co-located with its
  sole consumer (`src/ephemeris/vsop87/`, `src/ephemeris/elp2000/`,
  `src/astro/nutation/`, `src/ephemeris/`).
- **Removed `src/data/source/lagrange/`** — metadata superseded by
  `archive/lagrange/vsop87/manifest.toml`.
- `data::compiled` now contains only the feature-gated JPL DE440/DE441 binary
  sub-modules; all other compiled tables are gone.

### Fixed

- Fixed misspelled `bodies::satelite` path to `bodies::satellite`.
- Made RINEX NAV parsing strict by default, returning located `FormatError`s instead of silently zero-filling malformed records.
- Fixed ILRS CPF/CRD doctest import paths.
- Removed unreachable non-finite epoch tests now covered by stricter `JulianDate` construction.


## [0.8.0] - 2026-05-18

### Added

* `siderust::JulianDate` and `siderust::ModifiedJulianDate` re-exported at
  crate root so user code can write `siderust::JulianDate::new(…)` without
  an additional `use siderust::time::…` import.

### Changed

* **Dependency**: `qtty` upgraded from `0.7` to `0.8`.
* **Dependency**: `tempoch` upgraded from `0.4` to `0.6`.
  `siderust-ffi` dependencies updated to `qtty 0.8`, `tempoch 0.6`, and
  `tempoch-ffi 0.6`.
* `JulianDate::new(f64)` and `ModifiedJulianDate::new(f64)` are now the
  canonical constructors throughout the codebase; `from_raw_unchecked` and
  the fallible `try_new` variants are no longer used internally.
* Internal time-period durations are now computed via `.end()` / `.start()`
  accessors rather than raw-value arithmetic.
* Internal time calculations use the `julian_centuries()` method for
  consistency.
* All public code examples (doctests) updated to use the unambiguous
  `JulianDate::new(2451545.0)` / `ModifiedJulianDate::new(60000.0)` form.

### Removed

* **`siderust-ffi` dynamics module** (`siderust_ffi::dynamics`) and its
  entire C ABI surface have been removed.  The orbital/satellite-propagation
  FFI layer is no longer part of the crate.

### Fixed

* Resolved 17 doctest compilation errors caused by passing
  `qtty::Day::new(x)` (a typed quantity) where the `f64` argument of
  `JulianDate::new` / `ModifiedJulianDate::new` is expected.  Affected
  modules: `calculus::altitude`, `calculus::azimuth`,
  `calculus::lunar::phase`, `calculus::lunar::moon_equations`,
  `calculus::solar`, `coordinates::transform`, `astro::orbit`, and
  `targets`.

## [0.7.0] - 2026-05-07

This release contains **breaking API changes** throughout the public surface.
Callers must migrate before upgrading.

### Added

* **New dimensionless `qtty` newtypes** in `crate::qtty` (re-exported from
  `siderust::qtty`):
  `OpticalDepth` / `OpticalDepths`, `Airmass` / `Airmasses`,
  `Albedo` / `Albedos`, `IlluminationFraction` / `IlluminationFractions`,
  `Refractivity` / `Refractivities`, `CipCoordinate` / `CipCoordinates`.
  All carry dimensional safety at compile time — accidental raw-`f64` mixing
  is now a type error.
* **`Asteroid::with_albedo(albedo: Albedos) -> Self`** const builder.
  The optional `albedo: Option<Albedos>` field is now part of `Asteroid`,
  matching the existing pattern on `Planet` and `Satellite`.
* **`DEFAULT_SCALE_HEIGHT: Kilometers`** constant in
  `siderust::atmosphere::rayleigh` replaces the deprecated
  `DEFAULT_SCALE_HEIGHT_KM: f64` (which remains with `#[deprecated]` for one
  cycle).
* **Crate-wide academic documentation sweep**: every module now has structured
  `## Scientific scope`, `## Technical scope`, and `## References`
  (BibTeX-style citations) sections.
* **Authoring-convention reference** at `doc/conventions.md` formalises the
  documentation and API style rules now followed throughout the crate.

### Changed

* **`EarthOrientationModel` runtime enum removed** from
  `coordinates::transform::context` and from the `pub use` in
  `coordinates::transform`.
  Model selection is now fully static via phantom-typed dispatch:
  * `to_frame_as::<F2, Nut: NutationModel>(jd)` on `DirectionAstroExt`,
    `SphericalDirectionAstroExt`, `VectorAstroExt`, and `PositionAstroExt`.
  * `to_as::<C2, F2, Nut>(jd)` on `PositionAstroExt` for a simultaneous
    center + frame transform.
  * The default nutation marker (`DefaultNutationModel`) is `Iau2006A`.
    No runtime branch, no heap allocation — all dispatch is monomorphised.

* **`rayleigh_optical_depth_bodhaine99` signature** (in
  `siderust::atmosphere::rayleigh`) now takes typed quantities:
  * `pressure_hpa: f64` → `pressure: Hectopascals`
  * `scale_height_km: f64` → `scale_height: Kilometers`

* **`AtmosphereProfile.rayleigh_scale_height_km: f64`** renamed and retyped to
  `rayleigh_scale_height: Kilometers`.

* **`astro::cio::CipCio`** fields `x` and `y` retyped from `f64` to
  `CipCoordinates`. `cip_xy()` now returns `(CipCoordinates, CipCoordinates)`.
  `cio_locator_s` accepts `CipCoordinates` arguments.

* **`astro::iers_data::lookup`** now accepts `Days` instead of `f64` for the
  Julian-date offset argument, ensuring callers cannot silently pass a raw
  number in the wrong time scale.

* **`astro::orbit::PreparedOrbit.mean_motion`** is now typed as
  `AngularRate<Radian, Day>`.

* **`astro::conic::MeanMotionOrbit.mean_motion`** is now typed as
  `AngularRate<Degree, Day>`.

* **`calculus::lunar::phase::MoonPhaseGeometry.illuminated_fraction`**
  retyped from `f64` to `IlluminationFractions`.
  `Moon::illumination_above`, `illumination_below`, and `illumination_range`
  now accept `IlluminationFractions` thresholds instead of bare `f64`.

### Removed

* `EarthOrientationModel` enum and all `to_*_model(...)` helper methods
  (added in 0.6.1). Replaced by the phantom-typed `to_frame_as` / `to_as`
  API described above.

### Migration

#### 1 — Frame transforms (nutation model)

```rust
// 0.6.x (runtime enum, removed)
use siderust::coordinates::transform::EarthOrientationModel;
let equatorial = ecliptic.to_frame(jd, EarthOrientationModel::Iau2006A);

// 0.7.0 (compile-time phantom type)
use siderust::astro::nutation::Iau2006A;
let equatorial = ecliptic.to_frame_as::<EquatorialFrame, Iau2006A>(jd);

// Using the default (Iau2006A) via DefaultNutationModel:
use siderust::coordinates::transform::context::DefaultNutationModel;
let equatorial = ecliptic.to_frame_as::<EquatorialFrame, DefaultNutationModel>(jd);
```

#### 2 — Rayleigh optical depth

```rust
// 0.6.x (raw f64 parameters, removed)
let tau = rayleigh_optical_depth_bodhaine99(wavelength, co2_ppm, 744.0, latitude, 8.0);

// 0.7.0 (typed quantities)
use siderust::qtty::{Hectopascals, Kilometers};
let tau = rayleigh_optical_depth_bodhaine99(
    wavelength,
    co2_ppm,
    Hectopascals::new(744.0),
    latitude,
    Kilometers::new(8.0),
);
```

#### 3 — Rayleigh scale height constant

```rust
// 0.6.x (deprecated f64 constant)
use siderust::atmosphere::rayleigh::DEFAULT_SCALE_HEIGHT_KM;
let h: f64 = DEFAULT_SCALE_HEIGHT_KM;

// 0.7.0 (typed constant)
use siderust::atmosphere::rayleigh::DEFAULT_SCALE_HEIGHT;
let h: Kilometers = DEFAULT_SCALE_HEIGHT;
```

#### 4 — AtmosphereProfile field

```rust
// 0.6.x
profile.rayleigh_scale_height_km   // f64

// 0.7.0
profile.rayleigh_scale_height      // Kilometers
```

#### 5 — Illuminated fraction

```rust
// 0.6.x
let fraction: f64 = geometry.illuminated_fraction;

// 0.7.0
use siderust::qtty::IlluminationFractions;
let fraction: IlluminationFractions = geometry.illuminated_fraction;
let raw: f64 = fraction.value();
```

#### 6 — IERS lookup

```rust
// 0.6.x
let entry = iers_data::lookup(jd_value_f64);

// 0.7.0
use siderust::qtty::Days;
let entry = iers_data::lookup(Days::new(jd_value_f64));
```

### Fixed

* Resolved all `cargo doc --no-deps` intra-doc link warnings (64 warnings → 0):
  * Removed link brackets from private submodule references in `calculus::altitude`,
    `calculus::azimuth`, `calculus::solar`, `calculus::lunar`, and `astro::nutation`.
  * Removed link brackets from private function references in `calculus::solar`,
    `calculus::lunar`, `calculus::math_core::root_finding`, `astro::earth_rotation_provider`,
    and `astro::proper_motion`.
  * Qualified previously-unresolved links: `[crate::qtty]`, `[crate::targets::Target]`,
    `[crate::coordinates::frames::ICRF]`, `[crate::coordinates::frames::EclipticMeanJ2000]`,
    `[crate::astro::orientation::IauRotationParams]`, and `[crate::astro::precession::…]`.
  * Replaced non-existent `[IersEop::from_entries]` / `[IersEop::from_file]` references
    with accurate prose pointing to `IersEop::new`.
  * Converted feature-gated links (`atmosphere`, `runtime-data`, `spectra`, `tables`)
    to plain text to avoid unresolved-link warnings in the default doc build.
  * Fixed eight redundant explicit link targets (e.g. `[Foo](path::Foo)` → `[path::Foo]`).
  * Replaced unresolvable generic-parameter links (e.g. `` [`Direction<Horizontal>`] ``)
    with code-formatted backtick spans.
* Removed unused imports that triggered `dead_code` / `unused-imports` compiler warnings:
  `MJD` (seven files in `calculus/`), `CoordinateScale` in `calculus/jpl/eval.rs`, and
  `Second` in `coordinates/observation/observer_state.rs`.
* Replaced deprecated API calls throughout:
  * `EncodedTime::value()` → `.jd_value()` (six call sites).
  * `Rotation3::from_matrix` → `from_matrix_unchecked` (three const statics).
  * `SphericalPosition::new_raw` → `new_unchecked` (thirteen call sites).
  * `SphericalPosition::new_raw_with_params` → `new_unchecked_with_params` (four call sites).
  * `SphericalDirection::new_raw` → `new_unchecked` (nine call sites).
* Updated altitude benchmarks (`benches/solar_altitude.rs`, `moon_altitude.rs`,
  `star_altitude.rs`, `altitude_comparison.rs`) to use `ModifiedJulianDate::from_chrono`
  (replaces the removed `from_utc`) and `Period<ModifiedJulianDate>` (replaces `Period<MJD>`).
* Removed erroneous `Eq` derive from `ConicError`: the `OutOfRange { value: f64 }` variant
  contains an `f64` which does not implement `Eq`.
* Fixed `build.rs` to eliminate duplicate `jpl_daf` / `jpl_pipeline` / `jpl_spk` module
  declarations that caused `E0428` errors when both `de440` and `de441` features are enabled.
* Added `de441 = []` to `[features]` in `Cargo.toml` and wired full `de441` build support
  in `build.rs` (matching the existing `de440` pipeline), resolving `unexpected cfg condition
  value: de441` warnings and `siderust_mock_de441` cfg check-cfg errors.

### Added
* **Generic 1D and 2D gridded tables** under the new optional `tables`
  feature (`siderust::tables`). Provides `Grid1D<X, V>` and `Grid2D<X, Y, V>`
  with strict-monotonic axis validation, per-axis `OutOfRange` policy
  (clamp / zero / error), and untyped `algo::{linear_1d, bilinear,
  bilinear_unit}` `f64` kernels for callers needing bit-for-bit parity
  with `numpy`-style table lookups (e.g. NSB's Leinert zodiacal lookup).
* **Atmospheric optics primitives** under the new optional `atmosphere`
  feature (`siderust::atmosphere`). Provides:
  * `airmass(zenith: Radians, formula)` with the named
    `AirmassFormula::{PlaneParallel, Young1994, Rozenberg1966,
    KrisciunasSchaefer1991}` variants.
  * Rayleigh optical depth `rayleigh_optical_depth_bodhaine99`
    (Bodhaine et al. 1999) with parameterizable scale height, plus
    `rayleigh_phase`.
  * Mie / aerosol optical depth via the Patat 2011 power law
    (`MieParams { tau0, alpha, lambda_ref }`) with a
    `MieParams::PARANAL` preset.
  * Beer-Lambert `transmission(tau, airmass)`.
* **Shared `OutOfRange` and `Provenance` types** at the crate root
  (`siderust::interp::OutOfRange`, `siderust::provenance::{Provenance,
  DataSource}`) so the `spectra` and `tables` features can share
  validation, boundary, and dataset-metadata vocabulary. The
  `siderust::spectra` re-exports preserve the existing public API.
* **Johnson B/V passband convenience accessors** in
  `siderust::spectra::passbands`: `johnson_b()` and `johnson_v()`
  alias the Bessell (1990) realization (PASP 102, 1181), which
  Bessell & Murphy (2012, PASP 124, 140) recommend as the canonical
  Johnson–Cousins reference. Lets downstream consumers (e.g. NSB's
  zodiacal component) replace nearest-grid-point B/V wavelength hacks
  with real filter integrals without committing to a specific dataset
  module name.
* **Generic typed sampled spectra** under the new optional `spectra` feature
  (`siderust::spectra`). Provides `SampledSpectrum<X: Unit, Y: Unit, S = f64>`
  with strict-monotonic validation, configurable `Interpolation` (linear) and
  `OutOfRange` (clamp / zero / error) policies, trapezoidal `integrate` /
  `integrate_range` / `integrate_weighted` (returning unit-correct
  `Quantity<Prod<Y, X>>`), `Provenance` / `DataSource` metadata, and a
  generic ASCII two-column loader (`spectra::loaders::ascii::two_column`).
  Untyped `f64` numerical kernels are exposed under `spectra::algo` for
  callers that need bit-for-bit `numpy.interp`/trapezoid parity.
* **Conic orbit API** for propagation models that are not plain elliptic Kepler elements.
  * `astro::conic::{ConicOrbit, MeanMotionOrbit, ConicError}` plus re-exported `ConicKind`
  * `calculus::conic_equations::{calculate_conic_position, calculate_mean_motion_position}`
  * `ConicOrbit::position_at(...)` for elliptic and hyperbolic propagation (`ParabolicUnsupported` for `e == 1`)
* `PreparedOrbit` for repeated elliptic propagation with cached mean motion, precomputed orientation trig, and a fast `position_at(...)` path.
* `astro::units::{GaussianYear, GaussianYears, GAUSSIAN_YEAR}` for AU-day / Gaussian-year orbital period work.
* Crate-root re-exports for `KeplerianOrbit`, `PreparedOrbit`, `ConicOrbit`, `MeanMotionOrbit`, `ConicError`, `ConicKind`, `twilight`, and `Twilight`.
* **FFI orbit extensions** in `siderust-ffi`.
  * New C structs: `siderust_mean_motion_orbit_t` and `siderust_conic_orbit_t`
  * New propagation APIs: `siderust_kepler_position_ex`, `siderust_mean_motion_position`, and `siderust_conic_position`
  * New prepared-orbit handle lifecycle: `siderust_prepared_orbit_create`, `siderust_prepared_orbit_position`, and `siderust_prepared_orbit_destroy`
* New type-level IAU nutation model markers in `astro::nutation`: `Iau2000A`, `Iau2000B`, `Iau2006`, and `Iau2006A`, plus shared runtime identifiers via `NutationModelId`.
* Runtime Earth-orientation model selection via `EarthOrientationModel`, including ergonomic `to_*_model(...)` coordinate transform methods for model-specific rotation paths.
* New nutation model showcase example `examples/14_nutation_models.rs`, covering default transforms and custom `AstroContext::with_model::<...>()` usage.
* Checked-in FFI bindings matrix at `doc/ffi_bindings_matrix.md`, tracking canonical Rust concepts across `qtty-ffi`, `tempoch-ffi`, `siderust-ffi`, and the C++/Python/JS adapters.
* **FFI transform contexts and richer target payloads** in `siderust-ffi`.
  * New opaque `SiderustContext` handle plus `SiderustEarthOrientationModel`
    for model-sensitive frame and horizontal transforms.
  * New context-aware transform entry points and explicit TT/UT1 horizontal APIs.
  * New canonical `siderust_generic_target_create(...)` constructor that accepts
    tagged spherical-direction, spherical-position, or cartesian-position payloads.

### Changed
* `astro::orbit::Orbit` has been renamed to `KeplerianOrbit`, making the elliptic-only semantics explicit across public builders, body constants, examples, and `BodycentricParams`.
* Orbit construction now exposes validated `try_new(...)` entry points for `KeplerianOrbit`, `MeanMotionOrbit`, and `ConicOrbit`, while unchecked `new(...)` remains available for trusted compile-time constants.
* Orbital shape and orientation data now flow through `affn` conic primitives; `siderust` keeps the astronomy-specific epoch, anomaly, and propagation semantics on top.
* Terminology now consistently uses `argument_of_periapsis` / `arg_periapsis_deg` instead of perihelion-specific naming in both Rust and FFI orbit types.
* `siderust_kepler_position` is now documented as a legacy FFI entry point; `siderust_kepler_position_ex` should be preferred when the orbit reference center matters.
* Transform defaults now use `Iau2006A` as `DefaultNutationModel`, while keeping opt-in support for `Iau2000A`, `Iau2000B`, and `Iau2006` through typed model contexts.
* `siderust-ffi` is now aligned with `tempoch-ffi 0.4.x`, bumps its ABI version to `400`, and generates its C header into both Cargo `OUT_DIR` and `include/` for consistent adapter consumption.
* `siderust-ffi` now follows the `tempoch-ffi 0.4` scalar carrier contract for JD/MJD values and extends `SiderustSphericalPos` with an explicit `length_unit`, so target payloads preserve distance-unit semantics at the ABI boundary.

### Removed
* Legacy compatibility namespaces `coordinates::types::direction` and `coordinates::types::position`; use the explicit aliases or the `coordinates::spherical::{direction, position}` modules instead.
* `calculus::horizontal::equatorial_to_horizontal_with_ctx`; use `equatorial_to_horizontal_true_of_date_with_ctx` instead.

### Fixed
* FFI orbit propagation can now preserve the declared orbit reference center in output metadata via `siderust_kepler_position_ex`, instead of always tagging results as heliocentric.
* Implemented `serde::Serialize` / `serde::Deserialize` for `time::JulianDate`
  and `time::ModifiedJulianDate` under the `serde` feature flag; these types
  are now serialisable when derived in downstream structs (`ConicOrbit`,
  `MeanMotionOrbit`, `CoordinateWithPM`, etc.).
* Replaced manual `>= &&<=` range checks in `atmosphere::mie` tests with
  `RangeInclusive::contains` (`clippy::manual_range_contains`).
* Replaced `.max(0).min(35)` clamp pattern in `tables::grid2d` with
  `clamp(0, 35)` (`clippy::manual_clamp`).
* Removed unused `use serde_json;` import from `tests/test_serde.rs`
  (`clippy::unused_imports` with `--all-features`).

## [0.6.1] - 2026-05-09

## [0.6.0] - 08/03/2026

### Added
* **Moon phase module** (`calculus::lunar::phase`) with continuous photometric geometry and discrete 8-label classification
  * `moon_phase_geocentric<E>()`, phase angle, illuminated fraction, elongation, waxing flag (generic over `Ephemeris` backend)
  * `moon_phase_topocentric<E>()`, same with observer parallax correction
  * `MoonPhaseGeometry` struct with `.label()` → `MoonPhaseLabel` (8 classical names)
  * `find_phase_events<E>()`, locate New Moon, First Quarter, Full Moon, Last Quarter instants via scan + Brent root-finding
  * `MoonPhaseSeries<E>::sample()` / `sample_topocentric()`, batch sampling for plotting
  * `PhaseThresholds` for customizable label bin widths
  * `PhaseKind`, `PhaseEvent`, `PhaseSearchOpts` types
* Re-exports of the full phase API at crate root for ergonomic access
* Integration test suite `tests/test_moon_phase.rs` with layered validation (L1–L7)

* **Unified Azimuth API** (`calculus::azimuth`) providing azimuth-based
  event finding and interval queries analogous to the existing altitude API
  (crossings, extrema, and in-range periods).
  * `AzimuthProvider` trait: `azimuth_at(...)` and `azimuth_periods(&AzimuthQuery)`.
  * Free functions: `azimuth_periods()`, `azimuth_crossings()`, `azimuth_extrema()`,
    `azimuth_ranges()`, `in_azimuth_range()`, `outside_azimuth_range()`.
  * Types: `AzimuthQuery`, `AzimuthCrossingEvent`, `AzimuthExtremum`,
    `AzimuthExtremumKind`, and re-use of `CrossingDirection` semantics.
  * Body implementations for `Sun`, `Moon`, `Star<'_>`, and `direction::ICRS`.
  * Per-body scalar engines: `calculus::solar::sun_azimuth_rad`,
    `calculus::lunar::moon_azimuth_rad`, and
    `calculus::stellar::fixed_star_azimuth_rad` (delegates to
    `calculus::horizontal::star_horizontal`).
  * Robust handling of the 0°/360° discontinuity:
    - Crossing detection uses `sin(az − bearing)` to avoid wrap artifacts.
    - Range queries use a midpoint-cosine transform `cos(az − mid)`.
    - Extrema detection uses a stateful unwrapping closure that accumulates
      ±2π offsets and wraps results back into `[0°,360°)`.
  * Crate-root re-exports for ergonomic access and a new integration test
    suite `tests/test_azimuth_api.rs` covering trait and free-function APIs.

* **C ABI / FFI layer**: new workspace crate `siderust-ffi` exposing a flat C API (generated header `siderust-ffi/include/siderust_ffi.h` via `cbindgen`) for coordinate transforms, ephemeris queries, altitude/azimuth events and periods, moon-phase sampling, observatories, and target tracking.
* VSOP87 planets (Mercury–Neptune) now implement `AltitudePeriodsProvider` and `AzimuthProvider`, enabling altitude/azimuth queries and event finding for planets via the unified APIs.

### Changed
* Center-shift transformations are unified under `coordinates::transform::centers::TransformCenter` (`pos.to_center(...)` / `pos.to_center_with(...)`), replacing the previous per-center extension traits and legacy `to_*centric` modules.
* The topocentric expert API is now exposed as `coordinates::transform::to_topocentric_with_ctx(...)` (free function) instead of `ToTopocentricExt::to_topocentric_with_ctx(...)`.
* `Transform` blanket impl for `Position` is restricted to standard centers (`Params = ()`); parameterised centers (Topocentric/Bodycentric) require explicit parameters via `to_center(...)`.
* Example suite reorganized and expanded: renamed/renumbered examples under `examples/` and updated `examples/README.md`.
* Workspace and dependency updates: `siderust` is now a workspace with `siderust-ffi` member; dependency bumps (notably `qtty` 0.4 and `affn` 0.3.3); local development now patches `cheby` and `tempoch` to `dev-deps/` via `[patch.crates-io]`.
* Core ephemeris API refactor: `Ephemeris` methods `sun_barycentric`,
  `earth_barycentric` and `earth_heliocentric` now return a plain
  `Position<...>` instead of `CoordinateWithPM<Position<...>>`. This
  simplifies consumers that only needed instantaneous positions and
  eliminates unnecessary proper-motion wrappers in the coordinate
  transform pipeline.
* VSOP87 API changes: `vsop87a` and `vsop87e` generated helpers now
  return `Position<...>` directly (and `vsop87*_pos_vel` return
  `(Position, Velocity)`). Call sites across transforms, tests and
  examples were updated to use the direct `Position` return values.
* JPL body helpers under `calculus::jpl::bodies` now return
  `Position<...>` directly (instead of `CoordinateWithPM::new_static(...)`).
* `Trackable` / solar-system unit integration: VSOP87-backed
  `Trackable` implementations now wrap the returned `Position` into a
  `CoordinateWithPM::new_static(..., jd)` when a `CoordinateWithPM` is
  still required by the `Trackable` trait, keeping the public
  `Trackable` contract intact while removing the wrapper from the
  ephemeris primitives.
* Tests and examples: updated numerous `.get_position()` usages,
  dereferences and formatting to match the new return types for
  `vsop87a`/`vsop87e` and JPL helpers.

### Fixed
* Compilation and doc examples updated to reflect the API change
  (e.g. printing or inspecting a `Position` no longer accesses
  `.position` field). All tests run locally after the refactor.

## [0.5.2] - 19/02/2026

### Added
* IERS Earth Orientation Parameters integration via `astro::eop` with `IersEop`, `EopProvider`, `EopValues`, and `NullEop`
* Build-time IERS `finals2000A.all` ingestion pipeline under `scripts/iers/`, with generated embedded EOP tables
* `astro::earth_rotation` helpers for TT→UT1 conversion and GMST wrappers (`jd_ut1_from_tt`, `gmst_from_tt`, `gmst_from_tt_eop`)
* New Earth-orientation modules in `astro`: `cio`, `era`, `polar_motion`, and `light_deflection`
* IAU integration coverage in `tests/test_iau_compliance.rs` for the full GCRS→ITRS chain and model consistency checks
* Horizontal coordinate transforms (`EquatorialTrueOfDate` ↔ `Horizontal`) via `coordinates::transform::horizontal::{ToHorizontal, FromHorizontal}` with UT1+TT support
* Ecliptic-of-date direction transforms via `coordinates::transform::ecliptic_of_date::{ToEclipticTrueOfDate, FromEclipticTrueOfDate}`
* `SphericalDirectionAstroExt` for time-dependent frame transformations on `spherical::Direction<F>`
* New conversion examples: `examples/all_center_conversions.rs` and `examples/all_frame_conversions.rs`
* Regression suite `tests/test_high_precision_earth_rotation_regression.rs` with SOFA reference vectors for true-of-date horizontal and topocentric site-vector paths
* ECEF-first coordinate aliases in `coordinates::types`: `EcefSphericalDir` and `EcefPos`

### Changed
* Migrated precession to IAU 2006 and nutation to IAU 2000B across core astronomy and transform pipelines
* Updated sidereal-time computation to ERA-based IAU 2006 functions (`gmst_iau2006`, `gast_iau2006`) with explicit UT1/TT handling
* `AstroContext` now defaults to `IersEop` (`DefaultEop`) for EOP-aware transformations
* Horizontal/topocentric, lunar, and stellar calculations now consume the updated IAU/EOP-based Earth-rotation flow
* Geodetic/topocentric APIs now use `Geodetic<ECEF>` directly (including `to_topocentric*`, horizontal transforms, altitude providers, and observatory constants)
* Earth-fixed → equatorial rotation is now centralized in `astro::earth_rotation_provider` and shared across topocentric/observer pipelines
* `calculus::horizontal` true-of-date conversion now uses `GAST` (via `gast_iau2006`) instead of `GMST` for `EquatorialTrueOfDate` inputs
* New context-aware APIs: `ToTopocentricExt::to_topocentric_with_ctx(...)`
* New context-aware APIs: `ObserverState::topocentric_with_ctx(...)`
* New context-aware APIs: `calculus::horizontal::{geocentric_j2000_to_apparent_topocentric_with_ctx, equatorial_to_horizontal_true_of_date_with_ctx}`
* Topocentric translation and diurnal-velocity rotation now default to the full IAU 2006 terrestrial→celestial chain (`W · R3(ERA) · Q`) with EOP (`dut1`, `xp`, `yp`, `dx`, `dy`) from `AstroContext::default()`
* Fast/approx behavior remains selectable by supplying a context with `NullEop`
* Coordinate extension traits (`DirectionAstroExt`, `PositionAstroExt`, `VectorAstroExt`) now provide context-free IAU defaults, plus `_with` variants and `.using(&AstroContext)` for overrides
* Coordinate transformation APIs and examples migrated to `EclipticMeanJ2000` naming
* Coordinate transformation examples now primarily use spherical coordinates for frame/center conversions
* Altitude API docs now explicitly state that `ModifiedJulianDate` / `Period<MJD>` are interpreted on the TT axis
* Proper-motion API now requires explicit RA convention selection via `ProperMotion::from_mu_alpha(...)` (`µα`) or `ProperMotion::from_mu_alpha_star(...)` (`µα⋆ = µα cosδ`)

### Removed
* Legacy geodetic site type `coordinates::centers::ObserverSite` (and `ObserverSiteError`) in favor of `Geodetic<ECEF>`
* Legacy spherical geographic aliases (`Geographic`, `GeographicDir`, `GeographicPos`, `GeographicCartDir`, `GeographicCartPos`) in favor of `Geodetic<ECEF>` and `Ecef*` aliases

### Fixed
* Center-shift transforms now correctly convert AU shifts into the destination length unit
* Fixed `tests/test_astro_nights_roque_2026.rs` UTC-vs-TT regression by converting reference UTC-MJD endpoints (and test windows) to TT before comparison
* Fixed RA proper-motion propagation when using catalog `µα⋆` values by converting to true `µα` with epoch declination, with an explicit near-pole error guard (`cos(dec)≈0`)
* Corrected `astro::light_deflection::solar_deflection_magnitude` normalization to the physically consistent first-order model `Δθ = (2GM/c²R) * cot(θ/2)`:
  - at `90°` elongation and `1 AU`: `~0.00407″` (4.07 mas), not `~1.75″`
  - at the solar limb (`~959.63″`): `~1.75″`
* Updated light-deflection documentation and tests to match the corrected physics, including scalar-vs-vector consistency checks
* **Offline-by-default build**: VSOP87, ELP2000, and IERS EOP tables are now pre-generated and committed under `src/generated/`; `build.rs` no longer downloads anything during a normal build (fixes docs.rs build failure)
* Generated files are refreshed via `SIDERUST_REGEN=1 cargo build` or the helper script `scripts/update_generated_tables.sh`
* Added `.github/workflows/update-datasets.yml`, scheduled (weekly) workflow that regenerates committed tables and opens a PR when sources change

## [0.5.1] - 19/02/2026

### Changed
* Development-only git submodules (`affn`, `cheby`, `qtty`, `tempoch`) are now under `dev-deps/` and excluded from the published crate package.
* Added `dev-deps/README.md` clarifying that published `siderust` must always use crates.io dependencies (no committed `path` overrides).

## [0.5.0] - 12/02/2026

### Added

#### Ephemeris Backends
* `Ephemeris` trait abstraction (`calculus::ephemeris`) with pluggable backends and `DefaultEphemeris` auto-selection (DE441 > DE440 > VSOP87)
* `de440` Cargo feature with `De440Ephemeris` backend (JPL DE440, 1550–2650 CE)
* `de441` Cargo feature with `De441Ephemeris` backend (JPL DE441 part-2, extended coverage)
* Build-time DE440/DE441 pipelines under `scripts/jpl/` with NAIF `.bsp` support and Git LFS datasets
* `calculus::jpl` shared DE4xx infrastructure: `DeData` trait, `DeEphemeris<D>` generic backend, Chebyshev segment evaluation
* `cheby` workspace crate, Chebyshev polynomial toolkit (DCT fitting, Clenshaw evaluation, segment tables)
* Build system `SIDERUST_JPL_STUB` env var and `siderust_mock_de441` cfg flag for CI stub backends

#### Unified Altitude API
* `calculus::altitude` module: `AltitudePeriodsProvider` trait for finding time intervals when celestial bodies are within specific altitude ranges
* `AltitudePeriodsProvider` implementations for `Sun`, `Moon`, `Star<'_>`, and `direction::ICRS`
* Free functions: `crossings()`, `culminations()`, `altitude_ranges()`, `above_threshold()`, `below_threshold()`, `altitude_periods()`
* `AltitudeQuery`, `SearchOpts`, `CrossingEvent`, `CrossingDirection`, `CulminationEvent`, `CulminationKind` types
* Crate-root re-exports of the entire altitude API (`siderust::{above_threshold, crossings, ...}`)

#### Body-Specific Altitude Engines
* `calculus::stellar`, analytical sinusoidal model exploiting Earth's rotation for fixed-star altitude periods
* `calculus::lunar`, Moon altitude functions with topocentric parallax (`find_moon_above_horizon`, `find_moon_below_horizon`, `find_moon_altitude_range`)
* Moon Chebyshev cache (`moon_cache`) for optimized repeated ephemeris evaluation
* `calculus::horizontal`, shared equatorial→horizontal coordinate pipeline factored out of Sun/Moon engines

#### Numerical Engine
* `calculus::math_core` module: astronomy-agnostic numerical algorithms
  - `root_finding`: Brent's method with pre-computed endpoint values, plus bisection solver
  - `extrema`: golden-section minimiser/maximiser
  - `intervals`: interval assembly from roots
  - `bracketing`: seed/bracket generation policies

#### Celestial Bodies
* `bodies::Asteroid` type with `AsteroidBuilder`, `AsteroidClass` enum, and presets: `CERES_AST`, `BENNU`, `APOPHIS`
* `bodies::Comet` type with `CometBuilder`, `OrbitFrame` enum, `period_years()` helper, and presets: `HALLEY`, `ENCKE`, `HALE_BOPP`
* `calculus::pluto`, Meeus/Williams Pluto heliocentric ephemeris (42-43 periodic terms, ~0.5″ accuracy 1885–2099)

#### Coordinates & Frames
* `Galactic` reference frame (re-exported from `affn`)
* `coordinates::types` module with concise type aliases (`IcrsDir`, `EclipticDir`, `GeographicPos`, `HorizontalPos`, etc.) and prelude
* `coordinates::observation` module: `Astrometric<D>` / `Apparent<D>` wrapper types and `ObserverState` for explicit geometric/observed direction separation
* Horizontal coordinate convention helpers following IAU Alt-Az convention

#### Observatories
* La Silla Observatory (`observatories::LA_SILLA_OBSERVATORY`, ESO, Chile: −29.2584°, −70.7346°, 2400 m)

#### Examples
* New `jpl_precise_ephemeris`, unified DE440/DE441 backend comparison (replaces separate DE440/DE441 examples)
* New `altitude_periods_trait`, comprehensive `AltitudePeriodsProvider` trait demonstration
* New `compare_sun_moon_star`, generic body comparison via trait polymorphism
* New `night_quality_scoring`, practical observing planner scoring nights by darkness and Moon interference
* New `star_observability`, multi-star observing planner with visibility windows and peak altitudes
* New `find_night_periods_365day`, full-year astronomical night search with CLI support

#### Benchmarks
* New `ephemeris_comparison`, comparative benchmark: VSOP87 vs DE440 vs DE441 across all `Ephemeris` trait methods
* New `altitude_comparison`, comparative benchmark: Sun vs Moon vs Star for single-point eval and period searches (7/30/365-day)
* New `moon_altitude`, detailed lunar altitude benchmarks (single eval, above/below horizon, altitude ranges, algorithm comparison)
* New `star_altitude`, fixed-star altitude benchmarks (single eval, thresholds, crossings)
* New `elp2000`, ELP2000-82B evaluation benchmarks at multiple epochs
* New `de441`, DE441 ephemeris body-query benchmarks (feature-gated)

#### Tests
* Comprehensive ephemeris backend tests (`test_ephemeris.rs`, 879 lines) covering all backends and multiple epochs
* Unified altitude API tests (`test_altitude_api.rs`, `test_altitude_provider.rs`), crossings, culminations, thresholds, ranges
* Domain B invariant tests (`test_domain_b.rs`), aberration separation, Astrometric/Apparent states, topocentric parallax
* Stellar engine tests (`test_stellar.rs`), circumpolar, rise/set, never-visible edge cases
* Asteroid/Comet body tests (`test_asteroid.rs`)

### Changed
* `Satellite` type now uses `Cow<'a, str>` for the name field with `new_const()` and `new()` constructors
* `DefaultEphemeris` selection prefers DE441 when `de441` is enabled, then DE440, then VSOP87
* DE440/DE441 examples consolidated into single `jpl_precise_ephemeris.rs` with `#[cfg]` gates
* `solar_altitude_culminations` example removed (functionality covered by `find_night_periods_365day`)
* Examples README reorganized by theme: Getting Started, Observational Astronomy, Solar System, Ephemeris Backends, Serialization
* Benchmarks README updated to document comparative vs per-module benchmarks
* `time` module fully migrated to `tempoch` crate with typed time scales (`Time<S>`)
* Quantities throughout the codebase migrated from raw `f64` to `qtty` typed quantities

### Fixed
* README: stale import paths (`astro::JulianDate` → `time::JulianDate`), placeholder `#TBD` accuracy data removed, removed references to deleted `units` module, updated crate layout to reflect actual module structure
* `src/lib.rs` docs: removed false claims about `#![no_std]` support and `f128` quad precision
* CHANGELOG 0.3.0: merged duplicate "Changed" section, fixed typos ("Geographic", "Conversion")

## [0.4.0] - 19/01/2026

### Added
* External `qtty` crate dependency for dimensional analysis
* New `affn` geometry kernel integration for coordinate algebra (cartesian + spherical)
* Reference frames (`coordinates::frames`) and centers (`coordinates::centers`), plus astronomical type aliases in `coordinates::{cartesian,spherical}`
* Body-centric and topocentric transformations (including Geocentric ↔ Topocentric and Horizontal frame support)
* Transformation context/providers API (`coordinates::transform`) and observer-dependent helpers (`coordinates::observation`)
* Extensive coordinate documentation (`src/coordinates/README.md`) and new runnable examples under `examples/`
* Local CI helper script `ci-local.sh`
* Explicit equatorial frame split into `EquatorialMeanJ2000`, `EquatorialMeanOfDate`, and `EquatorialTrueOfDate`, including frame-bias (ICRS↔J2000), precession, and nutation rotations with helper matrices (e.g., `precession_rotation_from_j2000`)
* Velocity-aware aberration helpers (`apply/remove_aberration_*_with_velocity`) for supplying arbitrary observer velocities alongside the VSOP87E annual model
* Added the generic `calculus::events::altitude_periods` engine with Sun-specific wrappers (`calculus::solar::altitude_periods`), high-precision root-finding helpers, a new `examples/astronomical_night.rs`, and Roque de los Muchachos regression tests backed by JSON reference data for night/day/twilight windows
* Introduced the `time` module with the `TimeInstant` trait, generic `Period<T>` intervals, serde-serializable `ModifiedJulianDate`, and the `examples/time_periods.rs` showcase (plus the new `serde`/`serde_json` tooling for reference data)

### Changed
* Migrated from internal `units` module to external `qtty` crate (v0.2.0), updating APIs, docs, and examples
* Direction types are now frame-only (no reference center), and the coordinate module layout was consolidated around `Position`, `Displacement`, `Vector`, `Velocity`, and `Direction`
* Refactored coordinate transformations toward a unified `Transform`-based API with extension traits in `coordinates::prelude`
* Annual aberration now uses a Lorentz transform with exact AU/day light speed and the barycentric VSOP87E Earth velocity, operating in mean-J2000 equatorial space with machine-precision round-trips
* Frame rotations now include the ICRS↔J2000 frame bias, use mean obliquity of date for ecliptic transforms, and route equatorial transforms through the new mean-of-date/true-of-date frames

### Removed
* Internal `src/units/` module (angular.rs, frequency.rs, length.rs, mass.rs, power.rs, time.rs, unitless.rs, velocity.rs)
* Deprecated transformation APIs: `TransformToTopocentric` and the legacy `to_horizontal` helper
* Deprecated spherical coordinate extension modules replaced by the `coordinates::spherical` wrapper types

### Fixed
* Build/CI tooling improvements (Git LFS support, stub cleanup, improved local CI scripts)

## [0.3.1-3] - 04/09/2025

### Added
* LFS (vsop/elp) data

### Changed
* Build process offline

### Removed
* Unused imports

### Fixed
* Doc compilation

## [0.3.0] - 04/09/2025

### Added
- Coordinate macro checkers
- `spherical::Direction`, `cartesian::Direction`, `cartesian::Velocity`
- Build scripts to fetch VSOP87 coefficients and auto-generate static Rust arrays
- Conversion functions from/to `Position`/`Direction`
- `apply_aberration_to_direction` / `remove_aberration_from_direction`
- Compute VSOP87 velocity
- Implement `Unit` trait
- Expand bodies catalog (Lagrange Points and other major bodies)
- `Direction` implements `Display`

### Changed
- Geographic constructor order
- Refactor `radial_distance` → `distance`
- Refactor `distance_from_origin()` → `distance()`
- `SphericalCoord` → `spherical::Position`
- `CartesianCoord` → `cartesian::Position`
- Compute aberration using VSOP87 Earth velocity
- Spherical distance is no longer `Optional`
- `JulianDay` & `ModifiedJulianDay` → `JulianDate` & `ModifiedJulianDate`
- Coordinate constructors handle any angular/distance unit

### Removed
- `SphericalBuilder`
- Hardcoded VSOP coefficients
- Ron–Vondrák velocity series

### Fixed
- Mismatch addition of units in coordinate center transformation

## [0.2.0] - 31/05/2025

### Added
- Implement arithmetic operators `Add`, `Sub`, `Div`, `Mul` for `CartesianCoord`

### Changed
- Geocentric coordinates now account for aberration

## [0.1.0] - 15/05/2025

### Added
- Initial commit
- Target, Coordinates, Astro, Calculus, Observatories and Units modules
- AGPL-3 license
