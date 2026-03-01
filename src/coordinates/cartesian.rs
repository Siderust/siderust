// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Cartesian coordinate type aliases for common astronomical systems.
//!
//! This module provides convenient type aliases that combine algebraic cartesian
//! types with standard astronomical reference frames and centers.
//!
//! ## Semantic Types
//!
//! - [`Position`]: Affine points with center, frame, and unit
//! - [`Displacement`]: Displacement vectors with frame and unit (center-independent)
//! - [`Direction`]: Unit vectors with frame only (dimensionless)
//! - [`Velocity`]: Rate-of-change vectors with frame and unit
//! - [`Vector`]: Generic free vector (base type for Displacement and Velocity)

use crate::coordinates::{centers, frames};

/// Re-export the algebraic types
pub use affn::cartesian::*;

// =============================================================================
// Direction type aliases (frame-only, no center, dimensionless)
// =============================================================================

pub mod direction {
    use super::frames;
    pub use super::Direction;

    /// **EclipticMeanJ2000** cartesian direction (unit vector).
    pub type EclipticMeanJ2000 = Direction<frames::EclipticMeanJ2000>;
    /// **Equatorial mean J2000** cartesian direction (unit vector).
    pub type EquatorialMeanJ2000 = Direction<frames::EquatorialMeanJ2000>;
    /// **Equatorial mean of date** cartesian direction (unit vector).
    pub type EquatorialMeanOfDate = Direction<frames::EquatorialMeanOfDate>;
    /// **Equatorial true of date** cartesian direction (unit vector).
    pub type EquatorialTrueOfDate = Direction<frames::EquatorialTrueOfDate>;
    /// **Horizontal** cartesian direction (unit vector).
    pub type Horizontal = Direction<frames::Horizontal>;
    /// **Geocentric Earth-fixed (ECEF)** cartesian direction (unit vector).
    ///
    /// For geodetic (lon/lat/h) positions, use [`Geodetic<ECEF>`](crate::coordinates::centers::Geodetic)
    /// instead; this type is for unit vectors in the Earth-fixed frame.
    /// Cartesian unit vector in the **Earth-Centred Earth-Fixed** frame.
    pub type EcefCartDir = Direction<frames::ECEF>;
    /// **ICRS** cartesian direction (unit vector).
    pub type ICRS = Direction<frames::ICRS>;
    /// **ICRF** cartesian direction (unit vector).
    pub type ICRF = Direction<frames::ICRF>;
    /// **FK4 B1950** cartesian direction (unit vector).
    pub type FK4B1950 = Direction<frames::FK4B1950>;
    /// **TEME** cartesian direction (unit vector).
    pub type TEME = Direction<frames::TEME>;
    /// **GCRS** cartesian direction (unit vector).
    pub type GCRSDir = Direction<frames::GCRS>;
    /// **Galactic** cartesian direction (unit vector).
    pub type Galactic = Direction<frames::Galactic>;
}

// =============================================================================
// Displacement type aliases (frame + unit, no center)
// =============================================================================

pub mod displacement {
    use super::frames;
    pub use super::Displacement;

    /// **EclipticMeanJ2000** displacement vector.
    pub type EclipticMeanJ2000<U> = Displacement<frames::EclipticMeanJ2000, U>;
    /// **Equatorial mean J2000** displacement vector.
    pub type EquatorialMeanJ2000<U> = Displacement<frames::EquatorialMeanJ2000, U>;
    /// **Equatorial mean of date** displacement vector.
    pub type EquatorialMeanOfDate<U> = Displacement<frames::EquatorialMeanOfDate, U>;
    /// **Equatorial true of date** displacement vector.
    pub type EquatorialTrueOfDate<U> = Displacement<frames::EquatorialTrueOfDate, U>;
    /// **Horizontal** displacement vector.
    pub type Horizontal<U> = Displacement<frames::Horizontal, U>;
    /// **ICRS** displacement vector.
    pub type ICRS<U> = Displacement<frames::ICRS, U>;
    /// **ICRF** displacement vector.
    pub type ICRF<U> = Displacement<frames::ICRF, U>;
    /// **FK4 B1950** displacement vector.
    pub type FK4B1950<U> = Displacement<frames::FK4B1950, U>;
    /// **TEME** displacement vector.
    pub type TEME<U> = Displacement<frames::TEME, U>;
}

// =============================================================================
// Position type aliases (center + frame + unit)
// =============================================================================

pub mod position {
    pub use super::Position;
    use super::{centers, frames};

    /// **Heliocentric EclipticMeanJ2000** cartesian position.
    pub type EclipticMeanJ2000<U, C = centers::Heliocentric> =
        Position<C, frames::EclipticMeanJ2000, U>;
    /// **Geocentric Equatorial mean J2000** cartesian position.
    pub type EquatorialMeanJ2000<U, C = centers::Geocentric> =
        Position<C, frames::EquatorialMeanJ2000, U>;
    /// **Geocentric Equatorial mean of date** cartesian position.
    pub type EquatorialMeanOfDate<U, C = centers::Geocentric> =
        Position<C, frames::EquatorialMeanOfDate, U>;
    /// **Geocentric Equatorial true of date** cartesian position.
    pub type EquatorialTrueOfDate<U, C = centers::Geocentric> =
        Position<C, frames::EquatorialTrueOfDate, U>;
    /// **Topocentric Horizontal** cartesian position.
    pub type Horizontal<U, C = centers::Topocentric> = Position<C, frames::Horizontal, U>;
    /// **Geocentric Earth-Centered Earth-Fixed (ECEF)** cartesian position.
    ///
    /// For geodetic (lon/lat/h) positions, use [`Geodetic<ECEF>`](crate::coordinates::centers::Geodetic)
    /// instead; this type is for Cartesian XYZ in the Earth-fixed frame.
    /// The ellipsoid-correct conversion is [`to_cartesian`](affn::ellipsoidal::Position::to_cartesian).
    pub type EcefPos<U, C = centers::Geocentric> = Position<C, frames::ECEF, U>;
    /// **Barycentric ICRS** cartesian position.
    pub type ICRS<U, C = centers::Barycentric> = Position<C, frames::ICRS, U>;
    /// **Heliocentric ICRS** cartesian position.
    pub type HCRS<U> = Position<centers::Heliocentric, frames::ICRS, U>;
    /// **Geocentric ICRS** cartesian position.
    ///
    /// # Approximation
    ///
    /// This alias uses [`frames::ICRS`] as a first-order approximation for
    /// the Geocentric Celestial Reference System ([`frames::GCRS`]). The
    /// difference is < 1 mas for typical astronomical directions (neglected:
    /// geocentre offset, relativistic terms). For strictly IAU-correct GCRS,
    /// use `Position<Geocentric, frames::GCRS, U>` directly.
    pub type GCRS<U> = Position<centers::Geocentric, frames::ICRS, U>;
    /// **Topocentric ICRS** cartesian position.
    pub type TCRS<U> = Position<centers::Topocentric, frames::ICRS, U>;

    /// **Barycentric ICRF** cartesian position (BCRS).
    pub type BCRS<U> = Position<centers::Barycentric, frames::ICRF, U>;

    /// **Heliocentric J2000** equatorial cartesian position.
    pub type HeliocentricJ2000<U> = Position<centers::Heliocentric, frames::EquatorialMeanJ2000, U>;

    /// **Geocentric J2000** equatorial cartesian position (FK5-compatible ECI).
    pub type GeocentricJ2000<U> = Position<centers::Geocentric, frames::EquatorialMeanJ2000, U>;

    /// **Geocentric FK4 B1950** cartesian position.
    pub type FK4B1950<U> = Position<centers::Geocentric, frames::FK4B1950, U>;

    /// **Geocentric TEME** cartesian position (SGP4/TLE).
    pub type TEME<U> = Position<centers::Geocentric, frames::TEME, U>;

    /// **Geocentric ITRF** cartesian position.
    pub type ITRF<U> = Position<centers::Geocentric, frames::ITRF, U>;

    /// **Galactic** barycentric cartesian position.
    pub type Galactic<U> = Position<centers::Barycentric, frames::Galactic, U>;

    // --- Planetocentric body-fixed cartesian positions ---

    /// **Mercurycentric** body-fixed cartesian position.
    pub type MercuryFixed<U> = Position<centers::Mercurycentric, frames::MercuryFixed, U>;
    /// **Venuscentric** body-fixed cartesian position.
    pub type VenusFixed<U> = Position<centers::Venuscentric, frames::VenusFixed, U>;
    /// **Marscentric** body-fixed cartesian position.
    pub type MarsFixed<U> = Position<centers::Marscentric, frames::MarsFixed, U>;
    /// **Selenocentric** cartesian position in Moon principal axes frame.
    pub type MoonPrincipalAxes<U> = Position<centers::Selenocentric, frames::MoonPrincipalAxes, U>;
    /// **Jovicentric** System III cartesian position.
    pub type JupiterSystemIII<U> = Position<centers::Jovicentric, frames::JupiterSystemIII, U>;
    /// **Saturnocentric** body-fixed cartesian position.
    pub type SaturnFixed<U> = Position<centers::Saturnocentric, frames::SaturnFixed, U>;
    /// **Uranocentric** body-fixed cartesian position.
    pub type UranusFixed<U> = Position<centers::Uranocentric, frames::UranusFixed, U>;
    /// **Neptunocentric** body-fixed cartesian position.
    pub type NeptuneFixed<U> = Position<centers::Neptunocentric, frames::NeptuneFixed, U>;
    /// **Plutocentric** body-fixed cartesian position.
    pub type PlutoFixed<U> = Position<centers::Plutocentric, frames::PlutoFixed, U>;
}

// =============================================================================
// Velocity type aliases (frame + unit, no center)
// =============================================================================

pub mod velocity {
    use super::frames;
    pub use super::Velocity;

    /// **EclipticMeanJ2000** cartesian velocity vector.
    pub type EclipticMeanJ2000<U> = Velocity<frames::EclipticMeanJ2000, U>;
    /// **Equatorial mean J2000** cartesian velocity vector.
    pub type EquatorialMeanJ2000<U> = Velocity<frames::EquatorialMeanJ2000, U>;
    /// **Equatorial mean of date** cartesian velocity vector.
    pub type EquatorialMeanOfDate<U> = Velocity<frames::EquatorialMeanOfDate, U>;
    /// **Equatorial true of date** cartesian velocity vector.
    pub type EquatorialTrueOfDate<U> = Velocity<frames::EquatorialTrueOfDate, U>;
    /// **Horizontal** cartesian velocity vector.
    pub type Horizontal<U> = Velocity<frames::Horizontal, U>;
    /// **ICRS** cartesian velocity vector.
    pub type ICRS<U> = Velocity<frames::ICRS, U>;
    /// **ICRF** cartesian velocity vector.
    pub type ICRF<U> = Velocity<frames::ICRF, U>;
    /// **TEME** cartesian velocity vector.
    pub type TEME<U> = Velocity<frames::TEME, U>;
    /// **FK4 B1950** cartesian velocity vector.
    pub type FK4B1950<U> = Velocity<frames::FK4B1950, U>;
}
