// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Reference Frames Module
//!
//! This module re-exports astronomical reference frames defined in `affn`
//! and adds siderust-specific extensions including planetary body-fixed
//! frames and historical catalog frames.
//!
//! ## Predefined Frames (from `affn`)
//!
//! The following reference frames are available (from `affn`):
//!
//! - [`ICRS`]: International Celestial Reference System (quasi-inertial, fundamental reference).
//! - [`ICRF`]: International Celestial Reference Frame (ICRS realization via VLBI).
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
//! - [`Galactic`]: Galactic coordinate system (IAU 1958).
//! - [`GCRS`]: Geocentric Celestial Reference System.
//! - [`CIRS`]: Celestial Intermediate Reference System.
//! - [`TIRS`]: Terrestrial Intermediate Reference System.
//!
//! ## Siderust-Specific Frames
//!
//! All frames are now defined in `affn` and re-exported here.
//!
//! ## Planetary Body-Fixed Frames (from `affn`)
//!
//! The following body-fixed frame marker types are defined in `affn`
//! (behind `feature = "astro"`) and re-exported here. The generic
//! [`IauRotationParams`] type comes from [`crate::astro`], while body-specific
//! constants (e.g. `MARS_ROTATION`) are defined in [`crate::bodies::solar_system`]
//! and re-exported by [`planetary`]. `FrameRotationProvider`
//! implementations live in `transform::providers`.
//!
//! - [`MercuryFixed`]: Mercury IAU 2015 body-fixed rotation model.
//! - [`VenusFixed`]: Venus IAU 2015 body-fixed rotation model.
//! - [`MarsFixed`]: Mars IAU 2015 body-fixed rotation model.
//! - [`MoonPrincipalAxes`]: Moon principal axes (selenocentric) frame.
//! - [`JupiterSystemIII`]: Jupiter System III magnetic-field-based rotation.
//! - [`SaturnFixed`]: Saturn IAU body-fixed rotation model.
//! - [`UranusFixed`]: Uranus IAU body-fixed rotation model.
//! - [`NeptuneFixed`]: Neptune IAU body-fixed rotation model.
//! - [`PlutoFixed`]: Pluto IAU body-fixed rotation model.
//!
//! Each frame type provides its canonical name via [`ReferenceFrame::frame_name()`].
//! Planetocentric body-fixed frames use latitude/longitude/radius naming
//! via [`SphericalNaming`].
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

/// Planetary body-fixed frame markers and their IAU rotation-parameter constants.
pub mod planetary {
    // Re-export frame types from affn so existing `use …::planetary::MarsFixed`
    // paths continue to work.
    pub use affn::frames::{
        JupiterSystemIII, MarsFixed, MercuryFixed, MoonPrincipalAxes, NeptuneFixed, PlutoFixed,
        SaturnFixed, UranusFixed, VenusFixed,
    };

    // Re-export planetary IAU rotation parameters from solar-system body metadata.
    pub use crate::astro::{HasIauRotation, IauRotationParams};
    pub use crate::bodies::solar_system::{
        JUPITER_ROTATION, MARS_ROTATION, MERCURY_ROTATION, MOON_ROTATION, NEPTUNE_ROTATION,
        PLUTO_ROTATION, SATURN_ROTATION, URANUS_ROTATION, VENUS_ROTATION,
    };
}

// Re-export trait infrastructure from affn
pub use affn::frames::{ReferenceFrame, SphericalNaming};

// Re-export ellipsoid traits and predefined ellipsoids from affn
pub use affn::ellipsoid::{Ellipsoid, Grs80, HasEllipsoid, Wgs84};

// Re-export the generic ellipsoidal Position type for direct use
pub use affn::ellipsoidal::Position as EllipsoidalPosition;

// Re-export all astronomical frame types from affn
pub use affn::frames::{
    EclipticMeanJ2000, EclipticMeanOfDate, EclipticOfDate, EclipticTrueOfDate, EquatorialMeanJ2000,
    EquatorialMeanOfDate, EquatorialTrueOfDate, Galactic, Horizontal, CIRS, ECEF, FK4B1950, GCRS,
    ICRF, ICRS, ITRF, TEME, TIRS,
};

// Re-export siderust-specific frame types
pub use planetary::{
    JupiterSystemIII, MarsFixed, MercuryFixed, MoonPrincipalAxes, NeptuneFixed, PlutoFixed,
    SaturnFixed, UranusFixed, VenusFixed,
};

// NOTE: The `Horizontal` frame type uses the **North-clockwise** azimuth convention
// (0deg = North, increasing through East). If you are importing data that uses a
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
impl MutableFrame for FK4B1950 {}
impl MutableFrame for TEME {}
impl MutableFrame for Galactic {}
impl MutableFrame for GCRS {}
