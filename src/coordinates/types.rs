// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Concise re-exports of common coordinate types.
//!
//! This module provides short, unambiguous aliases for the most frequently used
//! coordinate types so that downstream code can avoid deep
//! `representation::semantic` import paths.
//!
//! # Naming convention
//!
//! | Suffix  | Meaning                                                    |
//! |---------|------------------------------------------------------------|
//! | `Dir`   | Spherical **direction** (frame-only, no center, unitless)  |
//! | `Pos`   | Spherical **position** (center + frame + distance unit)    |
//! | `CartDir`| Cartesian **direction** (unit vector, frame-only)         |
//! | `CartPos`| Cartesian **position** (affine point, center + frame + unit)|
//!
//! Spherical types are treated as the primary representation because they are
//! the most common in observational astronomy.
//!
//! # Example
//!
//! Before:
//! ```rust,ignore
//! use siderust::coordinates::spherical::direction::ICRS;
//! use siderust::coordinates::spherical::position::Geographic;
//! use siderust::coordinates::cartesian::position::Ecliptic;
//! ```
//!
//! After:
//! ```rust,ignore
//! use siderust::coordinates::types::{IcrsDir, GeographicPos, EclipticCartPos};
//! ```
//!
//! Or import everything at once:
//! ```rust,ignore
//! use siderust::coordinates::prelude::*;
//! ```

// =============================================================================
// Spherical direction aliases (frame-only, no center, dimensionless)
// =============================================================================

/// **Ecliptic** spherical direction (longitude, latitude).
pub use super::spherical::direction::Ecliptic as EclipticDir;
/// **Equatorial mean J2000** spherical direction (right-ascension, declination).
pub use super::spherical::direction::EquatorialMeanJ2000 as EquatorialJ2000Dir;
/// **Equatorial mean of date** spherical direction (right-ascension, declination).
pub use super::spherical::direction::EquatorialMeanOfDate as EquatorialMeanOfDateDir;
/// **Equatorial true of date** spherical direction (right-ascension, declination).
pub use super::spherical::direction::EquatorialTrueOfDate as EquatorialTrueOfDateDir;
/// **Galactic** spherical direction (l, b).
pub use super::spherical::direction::Galactic as GalacticDir;
/// **Geographic (ECEF)** spherical direction (longitude, latitude).
pub use super::spherical::direction::Geographic as GeographicDir;
/// **Horizontal** spherical direction (altitude, azimuth).
pub use super::spherical::direction::Horizontal as HorizontalDir;
/// **ICRS** spherical direction (right-ascension, declination).
pub use super::spherical::direction::ICRS as IcrsDir;

// =============================================================================
// Spherical position aliases (center + frame + distance unit)
// =============================================================================

/// **Heliocentric Ecliptic** spherical position (λ, β, R).
pub use super::spherical::position::Ecliptic as EclipticPos;
/// **Geocentric Equatorial mean J2000** spherical position (α, δ, d).
pub use super::spherical::position::EquatorialMeanJ2000 as EquatorialJ2000Pos;
/// **Geocentric Equatorial mean of date** spherical position.
pub use super::spherical::position::EquatorialMeanOfDate as EquatorialMeanOfDatePos;
/// **Geocentric Equatorial true of date** spherical position.
pub use super::spherical::position::EquatorialTrueOfDate as EquatorialTrueOfDatePos;
/// **Geographic (ECEF)** spherical position (longitude, latitude, altitude in km).
pub use super::spherical::position::Geographic as GeographicPos;
/// **Topocentric Horizontal** spherical position (Alt, Az, d).
pub use super::spherical::position::Horizontal as HorizontalPos;
/// **Geocentric ICRS** spherical position.
pub use super::spherical::position::GCRS as GcrsPos;
/// **Heliocentric ICRS** spherical position.
pub use super::spherical::position::HCRS as HcrsPos;
/// **Barycentric ICRS** spherical position.
pub use super::spherical::position::ICRS as IcrsPos;

// =============================================================================
// Cartesian direction aliases (unit vector, frame-only)
// =============================================================================

/// **Ecliptic** cartesian direction (unit vector).
pub use super::cartesian::direction::Ecliptic as EclipticCartDir;
/// **Equatorial mean J2000** cartesian direction (unit vector).
pub use super::cartesian::direction::EquatorialMeanJ2000 as EquatorialJ2000CartDir;
/// **Equatorial mean of date** cartesian direction (unit vector).
pub use super::cartesian::direction::EquatorialMeanOfDate as EquatorialMeanOfDateCartDir;
/// **Equatorial true of date** cartesian direction (unit vector).
pub use super::cartesian::direction::EquatorialTrueOfDate as EquatorialTrueOfDateCartDir;
/// **Geographic (ECEF)** cartesian direction (unit vector).
pub use super::cartesian::direction::Geographic as GeographicCartDir;
/// **Horizontal** cartesian direction (unit vector).
pub use super::cartesian::direction::Horizontal as HorizontalCartDir;
/// **ICRS** cartesian direction (unit vector).
pub use super::cartesian::direction::ICRS as IcrsCartDir;

// =============================================================================
// Cartesian position aliases (affine point, center + frame + unit)
// =============================================================================

/// **Heliocentric Ecliptic** cartesian position.
pub use super::cartesian::position::Ecliptic as EclipticCartPos;
/// **Geocentric Equatorial mean J2000** cartesian position.
pub use super::cartesian::position::EquatorialMeanJ2000 as EquatorialJ2000CartPos;
/// **Geocentric Equatorial mean of date** cartesian position.
pub use super::cartesian::position::EquatorialMeanOfDate as EquatorialMeanOfDateCartPos;
/// **Geocentric Equatorial true of date** cartesian position.
pub use super::cartesian::position::EquatorialTrueOfDate as EquatorialTrueOfDateCartPos;
/// **Geocentric Geographic (ECEF)** cartesian position.
pub use super::cartesian::position::Geographic as GeographicCartPos;
/// **Topocentric Horizontal** cartesian position.
pub use super::cartesian::position::Horizontal as HorizontalCartPos;
/// **Geocentric ICRS** cartesian position.
pub use super::cartesian::position::GCRS as GcrsCartPos;
/// **Heliocentric ICRS** cartesian position.
pub use super::cartesian::position::HCRS as HcrsCartPos;
/// **Barycentric ICRS** cartesian position.
pub use super::cartesian::position::ICRS as IcrsCartPos;
/// **Topocentric ICRS** cartesian position.
pub use super::cartesian::position::TCRS as TcrsCartPos;

// =============================================================================
// Backward-compatible legacy-style namespaces
// =============================================================================

/// Legacy-style spherical direction aliases using the original names.
///
/// This preserves imports like:
/// `use siderust::coordinates::types::direction::ICRS;`
pub mod direction {
    pub use super::EclipticDir as Ecliptic;
    pub use super::EquatorialJ2000Dir as EquatorialMeanJ2000;
    pub use super::EquatorialMeanOfDateDir as EquatorialMeanOfDate;
    pub use super::EquatorialTrueOfDateDir as EquatorialTrueOfDate;
    pub use super::GalacticDir as Galactic;
    pub use super::GeographicDir as Geographic;
    pub use super::HorizontalDir as Horizontal;
    pub use super::IcrsDir as ICRS;
}

/// Legacy-style spherical position aliases using the original names.
///
/// This preserves imports like:
/// `use siderust::coordinates::types::position::ICRS;`
pub mod position {
    pub use super::EclipticPos as Ecliptic;
    pub use super::EquatorialJ2000Pos as EquatorialMeanJ2000;
    pub use super::EquatorialMeanOfDatePos as EquatorialMeanOfDate;
    pub use super::EquatorialTrueOfDatePos as EquatorialTrueOfDate;
    pub use super::GcrsPos as GCRS;
    pub use super::GeographicPos as Geographic;
    pub use super::HcrsPos as HCRS;
    pub use super::HorizontalPos as Horizontal;
    pub use super::IcrsPos as ICRS;
}
