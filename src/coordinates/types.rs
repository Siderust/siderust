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
//! use siderust::coordinates::cartesian::position::EclipticMeanJ2000;
//! ```
//!
//! After:
//! ```rust,ignore
//! use siderust::coordinates::types::{IcrsDir, EclipticCartPos};
//! ```
//!
//! Or import everything at once:
//! ```rust,ignore
//! use siderust::coordinates::prelude::*;
//! ```

// =============================================================================
// Spherical direction aliases (frame-only, no center, dimensionless)
// =============================================================================

/// **ECEF** spherical direction (unit vector in Earth-fixed frame).
pub use super::spherical::direction::EcefDir as EcefSphericalDir;
/// **EclipticMeanJ2000** spherical direction (longitude, latitude).
pub use super::spherical::direction::EclipticMeanJ2000 as EclipticDir;
/// **Equatorial mean J2000** spherical direction (right-ascension, declination).
pub use super::spherical::direction::EquatorialMeanJ2000 as EquatorialJ2000Dir;
/// **Equatorial mean of date** spherical direction (right-ascension, declination).
pub use super::spherical::direction::EquatorialMeanOfDate as EquatorialMeanOfDateDir;
/// **Equatorial true of date** spherical direction (right-ascension, declination).
pub use super::spherical::direction::EquatorialTrueOfDate as EquatorialTrueOfDateDir;
/// **Galactic** spherical direction (l, b).
pub use super::spherical::direction::Galactic as GalacticDir;
/// **Horizontal** spherical direction (altitude, azimuth).
pub use super::spherical::direction::Horizontal as HorizontalDir;
/// **FK4 B1950** spherical direction (right-ascension, declination in FK4).
pub use super::spherical::direction::FK4B1950 as Fk4B1950Dir;
/// **GCRS** spherical direction (right-ascension, declination in GCRS).
pub use super::spherical::direction::GCRS as GcrsFrameDir;
/// **ICRS** spherical direction (right-ascension, declination).
pub use super::spherical::direction::ICRS as IcrsDir;
/// **TEME** spherical direction (right-ascension, declination in TEME).
pub use super::spherical::direction::TEME as TemeDir;

// =============================================================================
// Spherical position aliases (center + frame + distance unit)
// =============================================================================

/// **Heliocentric EclipticMeanJ2000** spherical position (λ, β, R).
pub use super::spherical::position::EclipticMeanJ2000 as EclipticPos;
/// **Geocentric Equatorial mean J2000** spherical position (α, δ, d).
pub use super::spherical::position::EquatorialMeanJ2000 as EquatorialJ2000Pos;
/// **Geocentric Equatorial mean of date** spherical position.
pub use super::spherical::position::EquatorialMeanOfDate as EquatorialMeanOfDatePos;
/// **Geocentric Equatorial true of date** spherical position.
pub use super::spherical::position::EquatorialTrueOfDate as EquatorialTrueOfDatePos;
// NOTE: `GeographicPos` (spherical position alias) has been removed.
// The old `Geographic = Position<Geocentric, ECEF, Kilometer>` was a correctness
// footgun: spherical `distance` is radial distance, not ellipsoidal height.
// Use `Geodetic<ECEF>` (= `ellipsoidal::Position<Geocentric, ECEF, U>`) for geodetic (lon/lat/h) positions.
/// **Topocentric Horizontal** spherical position (Alt, Az, d).
pub use super::spherical::position::Horizontal as HorizontalPos;
/// **Geocentric ICRS** spherical position.
pub use super::spherical::position::GCRS as GcrsPos;
/// **Heliocentric ICRS** spherical position.
pub use super::spherical::position::HCRS as HcrsPos;
/// **Barycentric ICRS** spherical position.
pub use super::spherical::position::ICRS as IcrsPos;

// =============================================================================
// New system-level spherical position aliases
// =============================================================================

/// **Galactic** barycentric spherical position (l, b, d).
pub use super::spherical::position::Galactic as GalacticPos;
/// **Geocentric J2000** equatorial spherical position (FK5-compatible ECI).
pub use super::spherical::position::GeocentricJ2000 as GeocentricJ2000Pos;
/// **Heliocentric Ecliptic J2000** spherical position (λ, β, R).
pub use super::spherical::position::HeliocentricEclipticJ2000 as HeliocentricEclipticJ2000Pos;
/// **Heliocentric J2000** equatorial spherical position.
pub use super::spherical::position::HeliocentricJ2000 as HeliocentricJ2000Pos;
/// **Barycentric ICRF** spherical position (BCRS).
pub use super::spherical::position::BCRS as BcrsPos;
/// **FK4 B1950** geocentric spherical position.
pub use super::spherical::position::FK4B1950 as Fk4B1950Pos;
/// **Geocentric TEME** spherical position (SGP4/TLE).
pub use super::spherical::position::TEME as TemePos;

// --- Planetocentric body-fixed spherical positions ---
/// **Jovicentric** System III spherical position.
pub use super::spherical::position::JupiterSystemIII as JupiterSystemIIIPos;
/// **Marscentric** body-fixed spherical position.
pub use super::spherical::position::MarsFixed as MarsFixedPos;
/// **Mercurycentric** body-fixed spherical position.
pub use super::spherical::position::MercuryFixed as MercuryFixedPos;
/// **Selenocentric** Moon principal axes spherical position.
pub use super::spherical::position::MoonPrincipalAxes as MoonPrincipalAxesPos;
/// **Neptunocentric** body-fixed spherical position.
pub use super::spherical::position::NeptuneFixed as NeptuneFixedPos;
/// **Plutocentric** body-fixed spherical position.
pub use super::spherical::position::PlutoFixed as PlutoFixedPos;
/// **Saturnocentric** body-fixed spherical position.
pub use super::spherical::position::SaturnFixed as SaturnFixedPos;
/// **Uranocentric** body-fixed spherical position.
pub use super::spherical::position::UranusFixed as UranusFixedPos;
/// **Venuscentric** body-fixed spherical position.
pub use super::spherical::position::VenusFixed as VenusFixedPos;

// =============================================================================
// Cartesian direction aliases (unit vector, frame-only)
// =============================================================================

/// **ECEF** cartesian direction (unit vector in Earth-fixed frame).
pub use super::cartesian::direction::EcefCartDir;
/// **EclipticMeanJ2000** cartesian direction (unit vector).
pub use super::cartesian::direction::EclipticMeanJ2000 as EclipticCartDir;
/// **Equatorial mean J2000** cartesian direction (unit vector).
pub use super::cartesian::direction::EquatorialMeanJ2000 as EquatorialJ2000CartDir;
/// **Equatorial mean of date** cartesian direction (unit vector).
pub use super::cartesian::direction::EquatorialMeanOfDate as EquatorialMeanOfDateCartDir;
/// **Equatorial true of date** cartesian direction (unit vector).
pub use super::cartesian::direction::EquatorialTrueOfDate as EquatorialTrueOfDateCartDir;
/// **Horizontal** cartesian direction (unit vector).
pub use super::cartesian::direction::Horizontal as HorizontalCartDir;
/// **ICRS** cartesian direction (unit vector).
pub use super::cartesian::direction::ICRS as IcrsCartDir;

// =============================================================================
// Cartesian position aliases (affine point, center + frame + unit)
// =============================================================================

/// **Geocentric ECEF** cartesian position (Earth-Centered Earth-Fixed).
pub use super::cartesian::position::EcefPos;
/// **Heliocentric EclipticMeanJ2000** cartesian position.
pub use super::cartesian::position::EclipticMeanJ2000 as EclipticCartPos;
/// **Geocentric Equatorial mean J2000** cartesian position.
pub use super::cartesian::position::EquatorialMeanJ2000 as EquatorialJ2000CartPos;
/// **Geocentric Equatorial mean of date** cartesian position.
pub use super::cartesian::position::EquatorialMeanOfDate as EquatorialMeanOfDateCartPos;
/// **Geocentric Equatorial true of date** cartesian position.
pub use super::cartesian::position::EquatorialTrueOfDate as EquatorialTrueOfDateCartPos;
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
// New system-level cartesian position aliases
// =============================================================================

/// **Galactic** barycentric cartesian position.
pub use super::cartesian::position::Galactic as GalacticCartPos;
/// **Heliocentric J2000** equatorial cartesian position.
pub use super::cartesian::position::HeliocentricJ2000 as HeliocentricJ2000CartPos;
/// **Barycentric ICRF** cartesian position (BCRS).
pub use super::cartesian::position::BCRS as BcrsCartPos;
/// **FK4 B1950** geocentric cartesian position.
pub use super::cartesian::position::FK4B1950 as Fk4B1950CartPos;
/// **Geocentric ITRF** cartesian position.
pub use super::cartesian::position::ITRF as ItrfCartPos;
/// **Geocentric TEME** cartesian position (SGP4/TLE).
pub use super::cartesian::position::TEME as TemeCartPos;

// --- Planetocentric body-fixed cartesian positions ---
/// **Jovicentric** System III cartesian position.
pub use super::cartesian::position::JupiterSystemIII as JupiterSystemIIICartPos;
/// **Marscentric** body-fixed cartesian position.
pub use super::cartesian::position::MarsFixed as MarsFixedCartPos;
/// **Mercurycentric** body-fixed cartesian position.
pub use super::cartesian::position::MercuryFixed as MercuryFixedCartPos;
/// **Selenocentric** cartesian position in Moon principal axes frame.
pub use super::cartesian::position::MoonPrincipalAxes as MoonPrincipalAxesCartPos;
/// **Neptunocentric** body-fixed cartesian position.
pub use super::cartesian::position::NeptuneFixed as NeptuneFixedCartPos;
/// **Plutocentric** body-fixed cartesian position.
pub use super::cartesian::position::PlutoFixed as PlutoFixedCartPos;
/// **Saturnocentric** body-fixed cartesian position.
pub use super::cartesian::position::SaturnFixed as SaturnFixedCartPos;
/// **Uranocentric** body-fixed cartesian position.
pub use super::cartesian::position::UranusFixed as UranusFixedCartPos;
/// **Venuscentric** body-fixed cartesian position.
pub use super::cartesian::position::VenusFixed as VenusFixedCartPos;

// =============================================================================
// Backward-compatible legacy-style namespaces
// =============================================================================

/// Legacy-style spherical direction aliases using the original names.
///
/// This preserves imports like:
/// `use siderust::coordinates::types::direction::ICRS;`
pub mod direction {
    pub use super::EclipticDir as EclipticMeanJ2000;
    pub use super::EquatorialJ2000Dir as EquatorialMeanJ2000;
    pub use super::EquatorialMeanOfDateDir as EquatorialMeanOfDate;
    pub use super::EquatorialTrueOfDateDir as EquatorialTrueOfDate;
    pub use super::Fk4B1950Dir as FK4B1950;
    pub use super::GalacticDir as Galactic;
    pub use super::GcrsFrameDir as GCRS;
    pub use super::HorizontalDir as Horizontal;
    pub use super::IcrsDir as ICRS;
    pub use super::TemeDir as TEME;
}

/// Legacy-style spherical position aliases using the original names.
///
/// This preserves imports like:
/// `use siderust::coordinates::types::position::ICRS;`
pub mod position {
    pub use super::BcrsPos as BCRS;
    pub use super::EclipticPos as EclipticMeanJ2000;
    pub use super::EquatorialJ2000Pos as EquatorialMeanJ2000;
    pub use super::EquatorialMeanOfDatePos as EquatorialMeanOfDate;
    pub use super::EquatorialTrueOfDatePos as EquatorialTrueOfDate;
    pub use super::Fk4B1950Pos as FK4B1950;
    pub use super::GalacticPos as Galactic;
    pub use super::GcrsPos as GCRS;
    pub use super::HcrsPos as HCRS;
    pub use super::HorizontalPos as Horizontal;
    pub use super::IcrsPos as ICRS;
    pub use super::TemePos as TEME;
}
