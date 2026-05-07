// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Astro Module
//!
//! Top-level entry point for siderust's astronomical calculations, gathering
//! the IAU-aligned models for Earth rotation, sidereal time, precession,
//! nutation, aberration, light deflection, proper motion, and orbital
//! mechanics under a single namespace.
//!
//! ## Scientific scope
//!
//! The submodules implement the standard chain of effects required to relate
//! catalogue (BCRS/GCRS) directions to the apparent positions seen by a
//! ground- or space-based observer at a given epoch: precession and nutation
//! of the Earth's equator, Earth Orientation Parameters, sidereal/Earth-
//! rotation angles, polar motion, stellar aberration, gravitational light
//! deflection, proper motion of stars, and Keplerian/conic orbital
//! propagation. Together they provide an IAU 2006/2000A-compliant
//! astrometric reduction pipeline.
//!
//! ## Technical scope
//!
//! Each submodule is a thin, typed wrapper around the relevant SOFA/IERS
//! algorithms expressed in terms of `qtty` quantities and `tempoch` time
//! scales rather than raw `f64`. Models are selected at compile time via
//! marker types (e.g. the [`nutation`] family) so that transform paths
//! remain zero-cost. The [`orientation`] re-export exposes the IAU pole and
//! prime-meridian rotation parameters used for planetary body frames.
//!
//! ## Submodules
//!
//! - [`aberration`]: stellar aberration via the full Lorentz transform.
//! - [`cio`]: CIP `(X, Y)` coordinates and CIO locator `s`.
//! - [`conic`], [`orbit`]: conic-section and Keplerian orbital mechanics.
//! - [`earth_rotation`], [`earth_rotation_provider`], [`era`]: TT↔UT1, ERA,
//!   and the composite ITRS→equatorial rotation.
//! - [`eop`]: Earth Orientation Parameters backed by `tempoch` EOP data.
//! - [`iers_data`]: IERS EOP data re-exported from `tempoch`.
//! - [`light_deflection`]: GR deflection by the Sun and planets.
//! - [`nutation`]: IAU 2000A/2000B/2006A nutation models.
//! - [`orientation`]: IAU body pole / prime-meridian parameters.
//! - [`polar_motion`]: polar-motion matrix `W`.
//! - [`precession`]: IAU 2006 (Fukushima-Williams) precession.
//! - [`proper_motion`]: stellar proper-motion propagation.
//! - [`sidereal`]: GMST/GAST/LST.
//! - [`units`]: astronomical units (e.g. Gaussian year) not covered by `qtty`.
//!
//! ## References
//!
//! * IERS Conventions (2010)
//! * IAU 2000 Resolution B1, IAU 2006 Resolution B1
//! * SOFA software collection

pub mod aberration;
pub mod cio;
pub mod conic;
pub mod earth_rotation;
pub mod earth_rotation_provider;
pub mod eop;
pub mod era;
pub mod iers_data;
pub mod light_deflection;
pub mod nutation;
pub mod orbit;
pub mod orientation;
pub mod polar_motion;
pub mod precession;
pub mod proper_motion;
pub mod sidereal;
pub mod units;

pub use orientation::{HasIauRotation, IauRotationParams};
