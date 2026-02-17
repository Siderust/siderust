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
//! - [`targets`]       : `Target<T>` tracking with time and proper motion
//! - [`time`]          : Time types and scale-based `Period<S>` / generic `Interval<T>`
//! - [`astro`]         : Utilities for aberration, nutation, precession, sidereal time, and event searches
//! - [`calculus`]      : Numerical kernels (VSOP87, ELP2000, Pluto, DE4xx, altitude API, root-finding)
//! - [`bodies`]        : Planets, stars, satellites, asteroids, comets, and built-in catalogs
//! - [`observatories`] : Predefined observatory locations (Roque, Paranal, Mauna Kea, La Silla)
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
//! println!("{}", mars.position);
//! ```
//!
//! For a runnable tour of the library, see the `examples/` directory.

pub mod astro;
pub mod bodies;
pub mod calculus;
pub mod coordinates;
pub mod observatories;
pub mod targets;
pub use tempoch as time;

pub(crate) mod macros;

// ---------------------------------------------------------------------------
// Convenience re‑exports: unified altitude API
// ---------------------------------------------------------------------------
pub use calculus::altitude::{
    above_threshold, altitude_periods as compute_altitude_periods, altitude_ranges,
    below_threshold, crossings, culminations, AltitudePeriodsProvider, AltitudeQuery,
    CrossingDirection, CrossingEvent, CulminationEvent, CulminationKind, SearchOpts,
};
