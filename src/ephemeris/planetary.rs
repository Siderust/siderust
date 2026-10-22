// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Planetary heliocentric ephemeris extension
//!
//! ## Scientific scope
//!
//! The core [`Ephemeris`](super::Ephemeris) trait covers Sun, Earth, and Moon —
//! the bodies required by the coordinate-transform and aberration pipelines.
//! Apparent physical orientation of Mars and Saturn additionally needs each
//! planet's heliocentric ecliptic state so geocentric directions can be formed
//! with the **same** compile-time backend that supplies Earth.
//!
//! ## Technical scope
//!
//! [`PlanetaryEphemeris`] extends [`Ephemeris`](super::Ephemeris) with Mars and
//! Saturn heliocentric positions in [`EclipticMeanJ2000`](crate::coordinates::frames::EclipticMeanJ2000)
//! and astronomical units. Backends that lack planetary data simply omit this
//! trait; orientation entry points require `E: PlanetaryEphemeris` and therefore
//! cannot silently fall back to another ephemeris. [`Vsop87Ephemeris`](super::Vsop87Ephemeris)
//! implements this trait. [`RuntimeEphemeris`](super::RuntimeEphemeris) currently
//! indexes only Sun/Earth/Moon segments and does **not** implement it.
//!
//! ## References
//!
//! - Bretagnon, P. & Francou, G. (1988). "Planetary theories in rectangular and
//!   spherical variables: VSOP87 solutions". *A&A*, 202, 309–315.
//! - Archinal, B. A. et al., IAU Working Group on Cartographic Coordinates and
//!   Rotational Elements reports.

use super::Ephemeris;
use crate::coordinates::{cartesian::Position, centers::Heliocentric, frames::EclipticMeanJ2000};
use crate::qtty::AstronomicalUnit;
use crate::time::JulianDate;

/// Compile-time ephemeris backend that can supply Mars and Saturn heliocentric
/// states for orientation geometry.
///
/// Implementors must evaluate both planets with the **same** theory/data set
/// used by their [`Ephemeris`](super::Ephemeris) Earth methods. Do not mix
/// VSOP87 Earth with an unrelated Mars/Saturn source behind this trait.
pub trait PlanetaryEphemeris: Ephemeris {
    /// Fallible Mars heliocentric ecliptic position (AU).
    #[inline]
    fn try_mars_heliocentric(
        jd: JulianDate,
    ) -> Result<Position<Heliocentric, EclipticMeanJ2000, AstronomicalUnit>, super::EphemerisError>
    {
        Ok(Self::mars_heliocentric(jd))
    }

    /// Mars heliocentric ecliptic position (AU).
    fn mars_heliocentric(
        jd: JulianDate,
    ) -> Position<Heliocentric, EclipticMeanJ2000, AstronomicalUnit>;

    /// Fallible Saturn heliocentric ecliptic position (AU).
    #[inline]
    fn try_saturn_heliocentric(
        jd: JulianDate,
    ) -> Result<Position<Heliocentric, EclipticMeanJ2000, AstronomicalUnit>, super::EphemerisError>
    {
        Ok(Self::saturn_heliocentric(jd))
    }

    /// Saturn heliocentric ecliptic position (AU).
    fn saturn_heliocentric(
        jd: JulianDate,
    ) -> Position<Heliocentric, EclipticMeanJ2000, AstronomicalUnit>;
}
