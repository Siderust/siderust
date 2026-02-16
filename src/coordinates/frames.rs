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
//! - [`Ecliptic`]: Ecliptic coordinate system (based on the plane of Earth's orbit at J2000).
//! - [`EclipticOfDate`]: Mean ecliptic coordinate system of date (siderust-specific).
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
    Ecliptic, EquatorialMeanJ2000, EquatorialMeanOfDate, EquatorialTrueOfDate, Galactic,
    Horizontal, CIRS, ECEF, GCRS, ICRF, ICRS, ITRF, TIRS,
};

// NOTE: The `Horizontal` frame type uses the **North-clockwise** azimuth convention
// (0° = North, increasing through East). If you are importing data that uses a
// different convention (e.g. South-origin or counter-clockwise), use the helpers in
// [`crate::coordinates::horizontal`] to convert before constructing `Horizontal`
// coordinates.

// =============================================================================
// Siderust-Specific Frame Definitions
// =============================================================================

use affn::prelude::ReferenceFrame as DeriveReferenceFrame;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Mean ecliptic coordinate system of date.
///
/// This frame uses the mean ecliptic plane (obliquity of date) without nutation.
/// Transformations to this frame are time-dependent and require IAU 2006 precession.
///
/// Use the traits in [`crate::coordinates::transform::ecliptic_of_date`] to convert
/// between this frame and equatorial frames.
#[derive(Debug, Copy, Clone, DeriveReferenceFrame)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct EclipticOfDate;

// =============================================================================
// MutableFrame: Marker for Transformable Frames
// =============================================================================

/// Marker trait for frames that support time-independent mutual transformations.
///
/// This trait is implemented for frames between which coordinate
/// transformations are currently supported (ICRS, Ecliptic, EquatorialMeanJ2000).
pub trait MutableFrame: ReferenceFrame {}

impl MutableFrame for ICRS {}
impl MutableFrame for Ecliptic {}
impl MutableFrame for EquatorialMeanJ2000 {}
