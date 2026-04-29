// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Siderust
//!
//! **Precision astronomy & satellite mechanics in Rust.**
//!
//! `siderust` is a research-grade astronomy toolkit focused on:
//! - Strongly-typed coordinates (center + frame + unit encoded in the type system)
//! - Ephemerides (VSOP87/ELP2000 always available; optional JPL DE4xx backends)
//! - Observation planning utilities (altitude periods, crossings, culminations)
//! - Time handling (via the `tempoch` crate, re-exported as [`time`])
//!
//! ## Features
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
//! - [`coordinates`]   : Cartesian & Spherical coordinate types and transformations between reference centers & frames
//! - [`targets`]       : `CoordinateWithPM<T>` + `Trackable` trait for targets
//! - [`time`]          : Time types and scale-based `Period<S>` / generic `Interval<T>`
//! - [`astro`]         : Utilities for aberration, nutation, precession, sidereal time, and event searches
//! - [`calculus`]      : Numerical kernels (VSOP87, ELP2000, Pluto, DE4xx, altitude API, root-finding)
//! - [`bodies`]        : Planets, stars, satellites, asteroids, comets, and built-in catalogs
//! - [`observatories`] : Predefined observatory locations (Roque, Paranal, Mauna Kea, La Silla)
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
//! use siderust::{
//!     bodies::Mars,
//!     time::JulianDate,
//! };
//! use chrono::prelude::*;
//!
//! // 1. Select an epoch (UTC now to JD)
//! let jd = JulianDate::from_utc(Utc::now());
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
pub mod bodies;
pub mod calculus;
pub mod coordinates;
pub mod data;
pub mod geometry;
pub mod interp;
pub mod observatories;
pub mod provenance;
pub mod qtty;
#[cfg(feature = "atmosphere")]
pub mod atmosphere;
#[cfg(feature = "spectra")]
pub mod spectra;
#[cfg(feature = "tables")]
pub mod tables;
pub mod targets;
pub mod time;

pub(crate) mod macros;

// ---------------------------------------------------------------------------
// Convenience re‑exports: unified azimuth API
// ---------------------------------------------------------------------------
pub use calculus::azimuth::{
    azimuth_crossings, azimuth_extrema, azimuth_periods as compute_azimuth_periods, azimuth_ranges,
    in_azimuth_range, outside_azimuth_range, AzimuthCrossingDirection, AzimuthCrossingEvent,
    AzimuthExtremum, AzimuthExtremumKind, AzimuthProvider, AzimuthQuery,
};
pub use calculus::solar::{twilight, Twilight};
pub use calculus::solar::{twilight_classification, TwilightPhase};

// ---------------------------------------------------------------------------
// Convenience re‑exports: unified altitude API
// ---------------------------------------------------------------------------
pub use affn::conic::ConicKind;
pub use astro::conic::{ConicError, ConicOrbit, MeanMotionOrbit};
pub use astro::orbit::{KeplerianOrbit, PreparedOrbit};
pub use calculus::altitude::{
    above_threshold, altitude_periods as compute_altitude_periods, altitude_ranges,
    below_threshold, crossings, culminations, AltitudePeriodsProvider, AltitudeQuery,
    CrossingDirection, CrossingEvent, CulminationEvent, CulminationKind, SearchOpts,
};

// ---------------------------------------------------------------------------
// Convenience re‑exports: lunar phase API
// ---------------------------------------------------------------------------
pub use calculus::lunar::phase::{
    find_phase_events, illumination_above, illumination_below, illumination_range,
    moon_phase_geocentric, moon_phase_topocentric, MoonPhaseGeometry, MoonPhaseLabel,
    MoonPhaseSeries, PhaseEvent, PhaseKind, PhaseSearchOpts, PhaseThresholds,
};

// ---------------------------------------------------------------------------
// Convenience re‑exports: target abstractions
// ---------------------------------------------------------------------------
pub use targets::{CoordinateWithPM, Trackable};
