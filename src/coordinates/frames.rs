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
//! - [`Equatorial`]: Equatorial coordinate system (right ascension and declination).
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

// =============================================================================
// Astronomical Reference Frames
// =============================================================================

/// International Celestial Reference System.
///
/// The fundamental celestial reference frame used in modern astronomy.
/// It is quasi-inertial and centered at the solar system barycenter.
/// The axes are defined by the positions of distant quasars.
#[derive(Debug, Copy, Clone, PartialEq, Eq, Default, DeriveReferenceFrame)]
pub struct ICRS;

/// Local horizon coordinate system.
///
/// A topocentric frame based on the observer's local horizon.
/// Uses altitude (elevation above horizon) and azimuth (bearing from north).
#[derive(Debug, Copy, Clone, PartialEq, Eq, Default, DeriveReferenceFrame)]
pub struct Horizontal;

/// Equatorial coordinate system.
///
/// Based on Earth's equator and the vernal equinox.
/// Uses right ascension (RA) and declination (Dec).
/// May be fixed to a specific epoch (e.g., J2000) or rotating with precession.
#[derive(Debug, Copy, Clone, PartialEq, Eq, Default, DeriveReferenceFrame)]
pub struct Equatorial;

/// Ecliptic coordinate system.
///
/// Based on the plane of Earth's orbit around the Sun.
/// Uses ecliptic longitude and latitude.
#[derive(Debug, Copy, Clone, PartialEq, Eq, Default, DeriveReferenceFrame)]
pub struct Ecliptic;

/// International Terrestrial Reference Frame.
///
/// A geocentric Earth-fixed frame that co-rotates with the Earth.
/// Used for geodetic and geophysical applications.
#[derive(Debug, Copy, Clone, PartialEq, Eq, Default, DeriveReferenceFrame)]
pub struct ITRF;

/// Earth-Centered, Earth-Fixed coordinate system.
///
/// A geocentric Cartesian coordinate system that rotates with the Earth.
/// The X-axis points to the intersection of the prime meridian and equator.
#[derive(Debug, Copy, Clone, PartialEq, Eq, Default, DeriveReferenceFrame)]
pub struct ECEF;

/// Galactic coordinate system.
///
/// Based on the plane of the Milky Way galaxy.
/// Uses galactic longitude and latitude, with the center
/// of the galaxy defining the origin of galactic longitude.
#[derive(Debug, Copy, Clone, PartialEq, Eq, Default, DeriveReferenceFrame)]
pub struct Galactic;

// =============================================================================
// MutableFrame: Marker for Transformable Frames
// =============================================================================

/// Marker trait for frames that support mutual transformations.
///
/// This trait is implemented for frames between which coordinate
/// transformations are currently supported (ICRS, Ecliptic, Equatorial).
pub trait MutableFrame: ReferenceFrame {}

impl MutableFrame for ICRS {}
impl MutableFrame for Ecliptic {}
impl MutableFrame for Equatorial {}
