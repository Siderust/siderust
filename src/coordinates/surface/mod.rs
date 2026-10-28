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
//! implements [`BodySurfaceModel`] with a documented reference surface and
//! body-fixed orientation. Earth, the Moon, Mars, Mercury, Venus, and Pluto
//! are provided. Gas giants do not implement the trait unless a mission
//! supplies an explicit reference surface (e.g. a 1-bar level).
//!
//! **Longitude** is always stored east-positive. This matches planetocentric /
//! body-fixed Cartesian conventions. External **planetographic** longitudes
//! may use a different sense (west-positive for Mars in NAIF SPICE); convert
//! explicitly with [`east_positive_from_planetographic`] or
//! [`BodySurfaceModel::east_positive_longitude`]. No hidden sign flips occur
//! inside constructors. Latitude uses [`LatitudeConvention`] independently.
//!
//! Negative heights are permitted until the inward reference-surface normal
//! reaches a curvature/centre singularity; invalid interior sites are rejected
//! during construction or deserialization. Cartesian Horizontal axes follow
//! affn's spherical convention: +X north, +Y east and +Z local up. At exact
//! zenith/nadir, azimuth is set to 0 by convention because it is physically
//! undefined.
//!
//! Lunar fixed orientation uses the NAIF `IAU_MOON` analytical Mean-Earth
//! model with 13 periodic arguments ([`crate::coordinates::frames::MoonIauMeanEarth`]),
//! **not** a DE-specific `MOON_PA` / binary-PCK `MOON_ME` realisation. The
//! legacy name `MoonPrincipalAxes` is a compatibility alias for the same
//! marker. Public epochs are TT; tempoch converts them to TDB for lunar
//! orientation.
//!
//! Geometric Alt/Az computed here does **not** include atmospheric
//! refraction, aberration, light deflection, light-time, or terrain
//! occultation. Those belong in the observation / event pipelines.
//!
//! ## Technical scope
//!
//! - [`ReferenceSurface`] — spherical or ellipsoidal reference surface
//! - [`LatitudeConvention`] — planetocentric vs planetographic latitude
//! - [`PlanetographicLongitudeSense`] — external planetographic longitude sense
//! - [`BodySurfaceModel`] — compile-time link of body → center, fixed frame,
//!   and default surface
//! - [`SurfaceSite`] — typed observer site on body `B` (e.g. `SurfaceSite::<Moon>`)
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
//! - NAIF, `pck00011.tpc`, IAU_MOON analytical Mean-Earth orientation;
//!   planetographic longitude sense via SPICE `pgrrec` / `plnsns`.
//! - See also [`crate::coordinates::transform::horizontal`] for the
//!   Earth-specific GAST path retained for backward compatibility.

pub mod body;
pub mod center;
pub mod enu;
pub mod error;
pub mod reference_surface;
pub mod site;

pub use body::{
    east_positive_from_planetographic, planetographic_from_east_positive, BodySurfaceModel,
    LatitudeConvention, PlanetographicLongitudeSense, EARTH_WGS84_SURFACE, MARS_IAU_ELLIPSOID,
    MERCURY_IAU_SPHERE, MOON_IAU_SPHERE, PLUTO_IAU_SPHERE, VENUS_IAU_SPHERE,
};
pub use center::{SurfaceObserver, SurfaceObserverParams};
pub use enu::{body_fixed_position, local_enu_basis, EnuBasis};
pub use error::SurfaceError;
pub use reference_surface::ReferenceSurface;
pub use site::SurfaceSite;
