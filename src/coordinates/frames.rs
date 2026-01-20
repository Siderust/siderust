// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Reference Frames Module
//!
//! This module defines astronomical reference frames for coordinate systems.
//! A reference frame specifies the orientation of the axes used to describe positions in space.
//!
//! ## Architecture
//!
//! All frame types implement the [`ReferenceFrame`] trait from `affn`, which provides
//! a common interface. The trait itself is re-exported from `affn` for convenience.
//!
//! ## Predefined Frames
//!
//! The following reference frames are provided:
//!
//! - [`ICRS`]: International Celestial Reference System (quasi-inertial, fundamental reference).
//! - [`Horizontal`]: Local horizon system (altitude-azimuth).
//! - [`EquatorialMeanJ2000`]: Mean equator/equinox of J2000 (FK5/J2000 mean).
//! - [`EquatorialMeanOfDate`]: Mean equator/equinox of date (precessed, no nutation).
//! - [`EquatorialTrueOfDate`]: True equator/equinox of date (precession + nutation).
//! - [`Ecliptic`]: Ecliptic coordinate system (based on the plane of Earth's orbit).
//! - [`ITRF`]: International Terrestrial Reference Frame (Earth-fixed).
//! - [`ECEF`]: Earth-Centered, Earth-Fixed (geocentric, rotating with the Earth).
//! - [`Galactic`]: Galactic coordinate system.
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
//!
//! ## Example
//!
//! ```rust
//! use siderust::coordinates::frames::{ReferenceFrame, ICRS};
//!
//! let name = ICRS::frame_name();
//! assert_eq!(name, "ICRS");
//! ```

// Re-export the core trait from affn
pub use affn::frames::ReferenceFrame;
// Import derive from prelude for use in this module
use affn::prelude::ReferenceFrame as DeriveReferenceFrame;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

// =============================================================================
// Astronomical Reference Frames
// =============================================================================

/// International Celestial Reference System.
///
/// The fundamental celestial reference frame used in modern astronomy.
/// It is quasi-inertial and centered at the solar system barycenter.
/// The axes are defined by the positions of distant quasars.
#[derive(Debug, Copy, Clone, PartialEq, Eq, Default, DeriveReferenceFrame)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct ICRS;

/// Local horizon coordinate system.
///
/// A topocentric frame based on the observer's local horizon.
/// Uses altitude (elevation above horizon) and azimuth (bearing from north).
#[derive(Debug, Copy, Clone, PartialEq, Eq, Default, DeriveReferenceFrame)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Horizontal;

/// Mean equator and equinox of J2000.0 (FK5/J2000 mean).
///
/// Earth-based mean equator/equinox at epoch J2000.0, with nutation removed.
/// This is the classic "J2000 equatorial" frame used by many catalogs.
#[derive(Debug, Copy, Clone, PartialEq, Eq, Default, DeriveReferenceFrame)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct EquatorialMeanJ2000;

/// Mean equator and equinox of date.
///
/// Earth-based mean equator/equinox at a given epoch (precession applied,
/// nutation removed). Requires a TT epoch for transformations.
#[derive(Debug, Copy, Clone, PartialEq, Eq, Default, DeriveReferenceFrame)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct EquatorialMeanOfDate;

/// True equator and equinox of date.
///
/// Earth-based true equator/equinox at a given epoch (precession + nutation).
/// Requires a TT epoch for transformations.
#[derive(Debug, Copy, Clone, PartialEq, Eq, Default, DeriveReferenceFrame)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct EquatorialTrueOfDate;

/// Ecliptic coordinate system.
///
/// Based on the plane of Earth's orbit around the Sun.
/// Uses ecliptic longitude and latitude.
#[derive(Debug, Copy, Clone, PartialEq, Eq, Default, DeriveReferenceFrame)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Ecliptic;

/// International Terrestrial Reference Frame.
///
/// A geocentric Earth-fixed frame that co-rotates with the Earth.
/// Used for geodetic and geophysical applications.
#[derive(Debug, Copy, Clone, PartialEq, Eq, Default, DeriveReferenceFrame)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct ITRF;

/// Earth-Centered, Earth-Fixed coordinate system.
///
/// A geocentric Cartesian coordinate system that rotates with the Earth.
/// The X-axis points to the intersection of the prime meridian and equator.
#[derive(Debug, Copy, Clone, PartialEq, Eq, Default, DeriveReferenceFrame)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct ECEF;

/// Galactic coordinate system.
///
/// Based on the plane of the Milky Way galaxy.
/// Uses galactic longitude and latitude, with the center
/// of the galaxy defining the origin of galactic longitude.
#[derive(Debug, Copy, Clone, PartialEq, Eq, Default, DeriveReferenceFrame)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Galactic;

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
