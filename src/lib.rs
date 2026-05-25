// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Siderust
//!
//! **Precision positional astrometry, atmospheric optics, ephemerides,
//! sky-brightness and satellite mechanics in Rust.**
//!
//! `siderust` is a research-grade astronomy toolkit built on two explicit pillars:
//!
//! 1. **Typed quantities everywhere** — every physical quantity at a public API
//!    boundary is a [`qtty`] newtype (e.g. `Meters`, `Radians`, `JulianDate`,
//!    `Airmass`, `OpticalDepth`, `Albedo`, `IlluminationFraction`, `Refractivity`,
//!    `CipCoordinate`, …).  Bare `f64` is confined to internal math kernels.
//!
//! 2. **Compile-time model selection** — when a function can use one of several
//!    algorithms or conventions (nutation theory, airmass formula, extinction
//!    parameterisation, …) the choice is expressed as a zero-sized phantom-type
//!    parameter, not a runtime enum.  For example:
//!    ```rust,ignore
//!    // picks Iau2006A nutation at compile time — no runtime dispatch
//!    let equatorial = ecliptic.to_frame_as::<Equatorial, Iau2006A>(ctx);
//!    ```
//!
//! Together these two properties make type errors in physical calculations into
//! compile-time errors and eliminate whole classes of unit-confusion bugs.
//!
//! See [`doc/conventions.md`](../doc/conventions.md) for the mandatory
//! authoring rules: the academic-style module-doc template
//! (*Scientific scope / Technical scope / References*), typed-quantity
//! guidelines, and the phantom-type model-selection pattern.
//!
//! ## Scientific scope
//!
//! - Strongly-typed coordinates (center + frame + unit encoded in the type system)
//! - Ephemerides (VSOP87/ELP2000 always available; optional JPL DE4xx backends)
//! - Observation planning (altitude periods, crossings, culminations, azimuth windows)
//! - Atmospheric refraction and extinction (airmass, optical depth, scattering)
//! - Sky brightness and lunar photometry (phase geometry, albedo, spectral radiance)
//! - Time handling (via the `tempoch` crate, re-exported as [`time`])
//! - Keplerian / conic orbits and satellite mechanics
//! - Spectral utilities (optional `spectra` feature)
//!
//! ## API pillars
//!
//! - **Coordinates**: `cartesian::{Position, Direction, Velocity, ...}` and
//!   `spherical::{Position, Direction}` parameterized by `Center`, `Frame`, and `Unit`.
//! - **Transforms**: frame rotations + center shifts with compile-time guarantees.
//! - **Altitude API**: `AltitudePeriodsProvider` + free functions to compute crossings,
//!   culminations, altitude ranges, and above/below-threshold windows.
//! - **Ephemeris backends**: `Ephemeris` trait with VSOP87/ELP2000 and optional DE440/DE441.
//! - **Serde**: optional `serde` feature for public types.
//!
//! ## Crate Modules
//!
//! - [`coordinates`]             : Cartesian & Spherical coordinate types and transformations; includes [`SkyGrid`] sampling utility
//! - [`targets`]                 : `CoordinateWithPM<T>` + `Trackable` trait for observation targets
//! - [`time`]                    : Time types, scales, and typed `tempoch::Interval<T>` values
//! - [`astro`]                   : Aberration, nutation, precession, sidereal time, conic helpers, event support, orbits, orbital mechanics
//! - [`ephemeris`]               : Ephemeris traits and backends (VSOP87, ELP2000, DE4xx, Pluto)
//! - [`event`]                   : Altitude/azimuth/lunar/solar/stellar event-search APIs
//! - [`numeric`]                 : Reusable numerical kernels (root-finding, extrema, intervals, bracketing)
//! - [`bodies`]                  : Planets, stars, satellites, asteroids, comets, and built-in catalogs
//! - [`observatories`]           : Predefined observatory locations (Roque, Paranal, Mauna Kea, La Silla)
//! - [`qtty`]                    : Re-exports of typed quantity newtypes from the `qtty` crate (including `OpticalDepth`, `Airmass`, `Albedo`, `IlluminationFraction`, `Refractivity`, `CipCoordinate`)
//! - [`data`]                    : Scientific dataset catalog, provenance, checksums, compiled tables, and optional runtime download/cache manager
//! - [`formats`]                 : Low-level binary file-format parsers (e.g. SPICE DAF/SPK); no dataset-catalog knowledge
//! - `atmosphere` *(optional)* : Atmospheric refraction, extinction, airmass, and optical-depth models (`atmosphere` feature)
//! - `spectra` *(optional)*    : Spectral response and photometric bandpass utilities (`spectra` feature)
//! - `tables` *(optional)*     : Tabulated data loaders (`tables` feature)
//!
//! ## Error-handling conventions
//!
//! `siderust` deliberately uses three distinct error-reporting patterns;
//! they are not interchangeable.  The choice tells you something about the
//! contract of the function:
//!
//! | Pattern                        | Meaning                                                        |
//! |--------------------------------|----------------------------------------------------------------|
//! | `-> Result<T, ConcreteError>`  | A *recoverable, domain-specific* failure (bad input, no convergence, dataset out of range). The error type is module-local and documented. |
//! | `-> Option<T>`                 | A *lookup* that may legitimately have no answer (e.g. no event in the requested window, parameter outside the table's interpolation range). The `None` semantics are documented in the `# Returns` section of every public `Option`-returning function. |
//! | `panic!` / `unwrap` / `expect` | A **bug or contract violation**: a non-finite Julian Date, an inverse of a singular rotation, an inconsistent type-system state, or a corrupt build-time table. Public functions that may panic carry a `# Panics` section. |
//!
//! When extending the public API, prefer `Result` over `Option` whenever the
//! caller might want to know *why* there is no answer (e.g. "out of EOP
//! range" vs "no event").  Reserve panics for true invariant violations,
//! and document the trigger in the `# Panics` section.
//!
//! ## Minimal Example
//!
//! ```rust
//! use siderust::{bodies::Mars, time::{JulianDate, JD, TT, Time, UTC}};
//! use chrono::prelude::*;
//!
//! // 1. Select an epoch (UTC now to JD)
//! let jd: JulianDate = Time::<UTC>::from_chrono(Utc::now()).to::<TT>().to::<JD>();
//!
//! // 2. Compute barycentric ecliptic coordinates via VSOP87
//! let mars = Mars::vsop87e(jd);
//!
//! // 3. Print Mars's barycentric ecliptic position (AstronomicalUnits)
//! println!("{:?}", mars);
//! ```
//!
//! For a runnable tour of the library, see the `examples/` directory.

pub(crate) use ::qtty as ext_qtty;

pub mod astro;
#[cfg(feature = "atmosphere")]
pub mod atmosphere;
pub mod bodies;
pub mod catalogs;
pub mod coordinates;
pub mod data;
pub mod ephemeris;
pub mod event;
pub mod formats;
pub mod instruments;
#[cfg(any(feature = "spectra", feature = "atmosphere"))]
pub(crate) mod interp;
pub mod mission_context;
pub mod mission_geometry;
pub mod numeric;
pub mod observatories;
#[cfg(feature = "pod")]
pub mod pod;
pub mod qtty;
#[cfg(feature = "spectra")]
pub mod spectra;
pub mod targets;
pub mod time;

// Ergonomic re-exports of common time markers / epoch (`siderust::J2000` in rustdoc examples).
pub use time::{JulianDate, ModifiedJulianDate, J2000, JD, MJD};

// Convenience re-export: sky sampling utilities.
pub use coordinates::{SkyGrid, SkyGridCell};

#[doc(hidden)]
pub use data::checksum;
pub(crate) mod macros;

// ---------------------------------------------------------------------------
// Convenience re‑exports: interval utilities
// ---------------------------------------------------------------------------
pub use numeric::intervals::intersect as intersect_periods;

// ---------------------------------------------------------------------------
// Convenience re‑exports: unified azimuth API
// ---------------------------------------------------------------------------
pub use event::azimuth::{
    azimuth_crossings, azimuth_extrema, azimuth_periods as compute_azimuth_periods, azimuth_ranges,
    in_azimuth_range, outside_azimuth_range, AzimuthCrossingDirection, AzimuthCrossingEvent,
    AzimuthExtremum, AzimuthExtremumKind, AzimuthProvider, AzimuthQuery,
};
pub use event::solar::{twilight, Twilight};
pub use event::solar::{twilight_classification, TwilightPhase};

// ---------------------------------------------------------------------------
// Convenience re‑exports: unified altitude API
// ---------------------------------------------------------------------------
pub use affn::conic::ConicKind;
pub use astro::conic::{ConicError, ConicOrbit, MeanMotionOrbit};
pub use astro::orbit::{KeplerianOrbit, PreparedOrbit};
pub use event::altitude::{
    above_threshold, altitude_periods as compute_altitude_periods, altitude_ranges,
    below_threshold, crossings, culminations, AltitudePeriodsProvider, AltitudeQuery,
    CrossingDirection, CrossingEvent, CulminationEvent, CulminationKind, SearchOpts,
};

// ---------------------------------------------------------------------------
// Convenience re‑exports: lunar phase API
// ---------------------------------------------------------------------------
pub use event::lunar::phase::{
    find_phase_events, illumination_above, illumination_below, illumination_range,
    moon_phase_geocentric, moon_phase_topocentric, MoonPhaseGeometry, MoonPhaseLabel,
    MoonPhaseSeries, PhaseEvent, PhaseKind, PhaseSearchOpts, PhaseThresholds,
};
pub use event::lunar::photometry::{
    lunar_albedo_jones2013, lunar_full_moon_albedo_jones2013, lunar_phase_attenuation_jones2013,
    reflected_lunar_spectral_radiance_jones2013,
};

// ---------------------------------------------------------------------------
// Convenience re‑exports: target abstractions
// ---------------------------------------------------------------------------
pub use targets::{CoordinateWithPM, Trackable};
