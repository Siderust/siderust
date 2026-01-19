# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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

### Deprecated

### Removed
* Internal `src/units/` module (angular.rs, frequency.rs, length.rs, mass.rs, power.rs, time.rs, unitless.rs, velocity.rs)
* Deprecated transformation APIs: `TransformToTopocentric` and the legacy `to_horizontal` helper
* Deprecated spherical coordinate extension modules replaced by the `coordinates::spherical` wrapper types

### Fixed
* Build/CI tooling improvements (Git LFS support, stub cleanup, improved local CI scripts)

### Security

## [0.3.1-3]

### Added
* LFS (vsop/elp) data

### Changed
* Build process offline

### Deprecated

### Removed
* unused imports

### Fixed
* Doc compilation

### Security

## [0.3.0]

### Added
- Coordinate macro checkers.
- spherical::Direction, cartesian::Direction, cartesian::Velocity.
- build scripts to fetch VSOP87 coeffs and auto-generate static rust arrays.
- Converstion functions from/to Position/Direction.
- apply_aberration_to_direction / remove_aberration_from_direction.
- Compute VSOP87 velocity.
- Implement Unit trait.
- Expand bodies catalog (Lagrange Points and other major bodies).
- Direction implements Display.

### Changed
- Geographic constructor order.
- spherical::Direction, cartesian::Direction.
- build scripts to fetch VSOP87 coeffs and auto-generate static rust arrays.
- Converstion functions from/to Position/Direction.
- apply_aberration_to_direction / remove_aberration_from_direction.
- Compute VSOP87 velocity

### Changed
- Geogreaphic constructor order.
- Refactor radial_distance -> distance.
- Refactor distance_from_origin() -> distance().
- SphericalCoord -> spherical::Position.
- CartesianCoord -> cartesian::Position.
- Compute Aberration using VSOP87 earth velocity.
- Spherical distance is no longer Optional.
- JulianDay & ModifiedJulianDay -> JuliaDate & ModifiedJulianDate.
- Coordinate constructors handle any angular/distance unit.

### Deprecated

### Removed
- SphericalBuilder.
- Hardcoded vsop coefficients.
- Ron–Vondrák velocity series.

### Fixed
- Missmatch addition of units in Coordinate center transformation.

### Security

## [0.2.0]

### Added
- Implement Arithmetic Operator Add, Sub, Div, Mult in CartesianCoord.

### Changed
- Geocentric coordinates now account for Aberration.

## [0.1.0]

### Added
- Initial commit.
- Target, Coordinates, Astro, Calculus, Observatories and Units modules.
- AGPL-3 license.
