// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Reference Frames Module
//!
//! This module re-exports astronomical reference frames defined in `affn`
//! and adds siderust-specific extensions.
//!
//! ## Predefined Frames
//!
//! The following reference frames are available (from `affn`):
//!
//! - [`ICRS`]: International Celestial Reference System (quasi-inertial, fundamental reference).
//! - [`Horizontal`]: Local horizon system (altitude-azimuth).
//! - [`EquatorialMeanJ2000`]: Mean equator/equinox of J2000 (FK5/J2000 mean).
//! - [`EquatorialMeanOfDate`]: Mean equator/equinox of date (precessed, no nutation).
//! - [`EquatorialTrueOfDate`]: True equator/equinox of date (precession + nutation).
//! - [`EclipticMeanJ2000`]: Mean ecliptic coordinate system at J2000.
//! - [`EclipticOfDate`]: Mean ecliptic coordinate system of date.
//! - [`EclipticMeanOfDate`]: Alias for [`EclipticOfDate`] (naming parity).
//! - [`EclipticTrueOfDate`]: True ecliptic coordinate system of date.
//! - [`ITRF`]: International Terrestrial Reference Frame (Earth-fixed).
//! - [`ECEF`]: Earth-Centered, Earth-Fixed (geocentric, rotating with the Earth).
//! - [`Galactic`]: Galactic coordinate system.
//!
//! Each frame type provides inherent named constructors and getters on
//! `Direction<F>` and `Position<C, F, U>`. For example:
//!
//! ```rust
//! use siderust::coordinates::frames::ICRS;
//! use affn::spherical::Direction;
//! use qtty::*;
//!
//! let d = Direction::<ICRS>::new(120.0 * DEG, 45.0 * DEG);
//! assert_eq!(d.ra(), 120.0 * DEG);
//! assert_eq!(d.dec(), 45.0 * DEG);
//! ```
//!
//! ## Extending
//!
//! To define a new reference frame, use the derive macro:
//!
//! ```rust
//! use affn::prelude::*;
//!
//! #[derive(Debug, Copy, Clone, ReferenceFrame)]
//! struct MyCustomFrame;
//! assert_eq!(MyCustomFrame::frame_name(), "MyCustomFrame");
//! ```

// Re-export trait infrastructure from affn
pub use affn::frames::{ReferenceFrame, SphericalNaming};

// Re-export all astronomical frame types from affn
pub use affn::frames::{
    EclipticMeanJ2000, EclipticMeanOfDate, EclipticOfDate, EclipticTrueOfDate, EquatorialMeanJ2000,
    EquatorialMeanOfDate, EquatorialTrueOfDate, Galactic, Horizontal, CIRS, ECEF, GCRS, ICRF, ICRS,
    ITRF, TIRS,
};

// NOTE: The `Horizontal` frame type uses the **North-clockwise** azimuth convention
// (0° = North, increasing through East). If you are importing data that uses a
// different convention (e.g. South-origin or counter-clockwise), use the helpers in
// [`crate::coordinates::horizontal`] to convert before constructing `Horizontal`
// coordinates.

// =============================================================================
// MutableFrame: Marker for Transformable Frames
// =============================================================================

/// Marker trait for frames that support time-independent mutual transformations.
///
/// This trait is implemented for frames between which coordinate
/// transformations are currently supported (ICRS, EclipticMeanJ2000, EquatorialMeanJ2000).
pub trait MutableFrame: ReferenceFrame {}

impl MutableFrame for ICRS {}
impl MutableFrame for EclipticMeanJ2000 {}
impl MutableFrame for EquatorialMeanJ2000 {}
