# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

### Changed

### Deprecated

### Removed

### Fixed

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
