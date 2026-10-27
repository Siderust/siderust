// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Generic celestial-body surface observers
//!
//! ## Scientific scope
//!
//! A surface observer is a point fixed relative to a celestial body's
//! reference surface and body-fixed frame. Geometric altitude and azimuth
//! are measured in the local east–north–up (ENU) system at that site:
//! altitude above the local geometric horizon (zenith = +90°) and azimuth
//! from north through east (0° = north, 90° = east).
//!
//! This module provides body-agnostic surface geometry for any body that
//! publishes a body-fixed orientation model and an explicit reference
//! surface. Earth, the Moon, and Mars are first-class; additional rocky
//! bodies reuse the same abstractions. Gas giants require an explicitly
//! documented reference surface (e.g. a 1-bar level) and are rejected by
//! default.
//!
//! Geometric Alt/Az computed here does **not** include atmospheric
//! refraction, aberration, light deflection, light-time, or terrain
//! occultation. Those belong in the observation / event pipelines.
//!
//! ## Technical scope
//!
//! - [`ReferenceSurface`] — spherical or ellipsoidal reference surface
//! - [`LatitudeConvention`] — planetocentric vs planetographic latitude
//! - [`BodySurfaceModel`] — compile-time link of body → center, fixed frame,
//!   and default surface
//! - [`SurfaceSite`] — typed observer site on body `B`
//! - [`SurfaceObserver`] — parameterized reference center carrying site params
//! - Local ENU basis and body-fixed Cartesian conversion
//!
//! Horizontal transforms live in
//! [`crate::coordinates::transform::surface_horizontal`].
//!
//! ## References
//!
//! - IAU Working Group on Cartographic Coordinates and Rotational Elements
//!   (2015). *Celestial Mechanics and Dynamical Astronomy* 130, 22.
//! - Archinal et al. (2018). Report of the IAU WGCCRE. *Celestial Mechanics
//!   and Dynamical Astronomy* 130, 22.
//! - See also [`crate::coordinates::transform::horizontal`] for the
//!   Earth-specific GAST path retained for backward compatibility.

pub mod body;
pub mod center;
pub mod enu;
pub mod error;
pub mod reference_surface;
pub mod site;

pub use body::{
    BodySurfaceModel, EarthSurface, LatitudeConvention, MarsSurface, MercurySurface, MoonSurface,
    PlutoSurface, VenusSurface, EARTH_WGS84_SURFACE, MARS_IAU_ELLIPSOID, MERCURY_IAU_SPHERE,
    MOON_IAU_SPHERE, PLUTO_IAU_SPHERE, VENUS_IAU_SPHERE,
};
pub use center::{SurfaceObserver, SurfaceObserverParams};
pub use enu::{body_fixed_position, local_enu_basis, EnuBasis};
pub use error::SurfaceError;
pub use reference_surface::ReferenceSurface;
pub use site::SurfaceSite;
