# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
* Regression suite `tests/test_high_precision_earth_rotation_regression.rs` with ERFA/SOFA reference vectors for true-of-date horizontal and topocentric site-vector paths
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
* `cheby` workspace crate — Chebyshev polynomial toolkit (DCT fitting, Clenshaw evaluation, segment tables)
* Build system `SIDERUST_JPL_STUB` env var and `siderust_mock_de441` cfg flag for CI stub backends

#### Unified Altitude API
* `calculus::altitude` module: `AltitudePeriodsProvider` trait for finding time intervals when celestial bodies are within specific altitude ranges
* `AltitudePeriodsProvider` implementations for `Sun`, `Moon`, `Star<'_>`, and `direction::ICRS`
* Free functions: `crossings()`, `culminations()`, `altitude_ranges()`, `above_threshold()`, `below_threshold()`, `altitude_periods()`
* `AltitudeQuery`, `SearchOpts`, `CrossingEvent`, `CrossingDirection`, `CulminationEvent`, `CulminationKind` types
* Crate-root re-exports of the entire altitude API (`siderust::{above_threshold, crossings, ...}`)

#### Body-Specific Altitude Engines
* `calculus::stellar` — analytical sinusoidal model exploiting Earth's rotation for fixed-star altitude periods
* `calculus::lunar` — Moon altitude functions with topocentric parallax (`find_moon_above_horizon`, `find_moon_below_horizon`, `find_moon_altitude_range`)
* Moon Chebyshev cache (`moon_cache`) for optimized repeated ephemeris evaluation
* `calculus::horizontal` — shared equatorial→horizontal coordinate pipeline factored out of Sun/Moon engines

#### Numerical Engine
* `calculus::math_core` module: astronomy-agnostic numerical algorithms
  - `root_finding`: Brent's method with pre-computed endpoint values, plus bisection solver
  - `extrema`: golden-section minimiser/maximiser
  - `intervals`: interval assembly from roots
  - `bracketing`: seed/bracket generation policies

#### Celestial Bodies
* `bodies::Asteroid` type with `AsteroidBuilder`, `AsteroidClass` enum, and presets: `CERES_AST`, `BENNU`, `APOPHIS`
* `bodies::Comet` type with `CometBuilder`, `OrbitFrame` enum, `period_years()` helper, and presets: `HALLEY`, `ENCKE`, `HALE_BOPP`
* `calculus::pluto` — Meeus/Williams Pluto heliocentric ephemeris (42-43 periodic terms, ~0.5″ accuracy 1885–2099)

#### Coordinates & Frames
* `Galactic` reference frame (re-exported from `affn`)
* `coordinates::types` module with concise type aliases (`IcrsDir`, `EclipticDir`, `GeographicPos`, `HorizontalPos`, etc.) and prelude
* `coordinates::observation` module: `Astrometric<D>` / `Apparent<D>` wrapper types and `ObserverState` for explicit geometric/observed direction separation
* Horizontal coordinate convention helpers following IAU Alt-Az convention

#### Observatories
* La Silla Observatory (`observatories::LA_SILLA_OBSERVATORY` — ESO, Chile: −29.2584°, −70.7346°, 2400 m)

#### Examples
* New `jpl_precise_ephemeris` — unified DE440/DE441 backend comparison (replaces separate DE440/DE441 examples)
* New `altitude_periods_trait` — comprehensive `AltitudePeriodsProvider` trait demonstration
* New `compare_sun_moon_star` — generic body comparison via trait polymorphism
* New `night_quality_scoring` — practical observing planner scoring nights by darkness and Moon interference
* New `star_observability` — multi-star observing planner with visibility windows and peak altitudes
* New `find_night_periods_365day` — full-year astronomical night search with CLI support

#### Benchmarks
* New `ephemeris_comparison` — comparative benchmark: VSOP87 vs DE440 vs DE441 across all `Ephemeris` trait methods
* New `altitude_comparison` — comparative benchmark: Sun vs Moon vs Star for single-point eval and period searches (7/30/365-day)
* New `moon_altitude` — detailed lunar altitude benchmarks (single eval, above/below horizon, altitude ranges, algorithm comparison)
* New `star_altitude` — fixed-star altitude benchmarks (single eval, thresholds, crossings)
* New `elp2000` — ELP2000-82B evaluation benchmarks at multiple epochs
* New `de441` — DE441 ephemeris body-query benchmarks (feature-gated)

#### Tests
* Comprehensive ephemeris backend tests (`test_ephemeris.rs`, 879 lines) covering all backends and multiple epochs
* Unified altitude API tests (`test_altitude_api.rs`, `test_altitude_provider.rs`) — crossings, culminations, thresholds, ranges
* Domain B invariant tests (`test_domain_b.rs`) — aberration separation, Astrometric/Apparent states, topocentric parallax
* Stellar engine tests (`test_stellar.rs`) — circumpolar, rise/set, never-visible edge cases
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
