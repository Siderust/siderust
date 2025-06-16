//! # Siderust
//!
//! **Precision astronomy & satellite mechanics in Rust.**
//!
//! The `siderust` crate provides high-accuracy celestial and orbital
//! calculations with zero-allocation performance and strong type safety.
//!
//! ## Features
//!
//! - **Coordinate Systems**: Strongly-typed `CartesianCoord` and `SphericalCoord` types,
//!   parameterized by phantom `ReferenceCenter` (e.g., Sun, Earth) and
//!   `ReferenceFrame` (e.g., ICRS, Ecliptic, Equatorial, Horizontal), and
//!   `Kind` (Position, Direction), enabling compile-time protection
//!   against mismatched frames and origins. Seamless conversions between
//!   coordinate systems are supported via `From`/`Into` and the `Transform` trait.
//! - **Target Tracking**: `Target<T>` couples any coordinate type with
//!   an observation `JulianDay` and optional `ProperMotion` for
//!   extrapolation and movement analysis.
//! - **Units & Quantities**: Strongly-typed `Mass`, `Length`, `Angle`,
//!   `Velocity`, `Duration`, etc., with operator overloading for
//!   natural arithmetic while preventing unit mistakes.
//! - **Astronomical Utilities**: Aberration, nutation, precession,
//!   sidereal time, apparent Sun/Moon positions, and rise/culmination
//!   search routines for celestial bodies and satellites.
//! - **Celestial Bodies & Catalogs**: Built-in data for Sun through
//!   Neptune, major moons, a starter star catalog, and helpers to load
//!   Gaia/Hipparcos datasets or custom star catalogs.
//! - **Numerical Kernels**: Kepler equation solvers, VSOP87 & ELP2000
//!   theories for planetary & lunar coordinates, and light-time
//!   corrections, validated against JPL Horizons & IMCCE to <1 mas.
//! - **Performance**: Zero heap allocations in core routines, SIMD
//!   optimizations, optional `f128` quad precision, and `#![no_std]`
//!   support with `libm` fallback.
//!
//! ## Crate Modules
//!
//! - [`units`]         : Strongly-typed physical quantities and angles
//! - [`coordinates`]   : Cartesian & Spherical coordinate types and
//!                       transformations between reference centers & frames
//! - [`targets`]       : `Target<T>` tracking with time and proper motion
//! - [`astro`]         : Utilities for aberration, nutation, precession,
//!                       sidereal time, and event searches
//! - [`calculus`]      : Numerical kernels (Kepler, VSOP87, ELP2000, etc.)
//! - [`bodies`]        : Data structures for planets, comets, stars,
//!                       satellites, and built-in catalogs
//! - [`observatories`] : Predefined observatory locations and helpers
//!
//! ## Minimal Example
//!
//! ```rust
//! use siderust::{
//!     bodies::Mars,
//!     units::JulianDay,
//! };
//! use chrono::prelude::*;
//!
//! // 1. Select an epoch (UTC now to JD)
//! let jd = JulianDay::from_utc(Utc::now());
//!
//! // 2. Compute heliocentric ecliptic coordinates via VSOP87
//! let mars = Mars::vsop87e(jd);
//!
//! // 3. Print Mars's heliocentric ecliptic position (AU)
//! println!("{}", mars.position);
//! ```
//!
//! ---
//!
//! **Note:** This documentation is generated and reviewed under @vpramon supervision.
//! For detailed usage and API, see the module-level docs and examples.

pub mod units;
pub mod coordinates;
pub mod targets;
pub mod astro;
pub mod calculus;
pub mod bodies;
pub mod observatories;

pub(crate) mod macros;
