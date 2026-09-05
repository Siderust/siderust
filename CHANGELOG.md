# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Add `ObservatoryCatalog` with deterministic bundled catalog access and
  validated runtime loading from `[[observatory]]` TOML records.

### Changed

- Generate the existing named observatory constants from the canonical
  `data/observatories.toml` source while preserving their public names and
  scientific values. `Observatory::name` is now `Cow<'static, str>` so the
  same domain type can safely own names loaded at runtime. This is a
  source-breaking Rust API change for callers that construct `Observatory`
  with struct literals and must ship in a SemVer-compatible pre-1.0 minor
  release.

## [0.11.1] - 2026-09-03

### Added

- Add ergonomic `EquatorialTrueOfDate` ↔ `Horizontal` transformations for
  Cartesian and spherical directions, including precise UT1+TT variants,
  reverse transformations, and topocentric position aliases.
- Add `HealpixGrid::pixel_center_spherical` for retrieving HEALPix pixel
  centers directly as typed spherical coordinates.

### Changed

- Align `siderust` and `siderust-ffi` release metadata with version `0.11.1`
  and report FFI version `1101` from `siderust_ffi_version()`.
