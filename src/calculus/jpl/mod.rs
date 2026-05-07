// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # JPL Planetary and Lunar Ephemerides (DE4xx)
//!
//! ## Scientific scope
//!
//! The JPL DE (*Development Ephemeris*) series are numerical integrations of
//! the equations of motion for the Solar System bodies produced by the Jet
//! Propulsion Laboratory.  Each release (DE440, DE441, …) provides highly
//! accurate Chebyshev polynomial representations of barycentric and
//! heliocentric positions and velocities.
//!
//! The two variants supported here differ in time span and perturbation model:
//!
//! - **DE440** — modern fit (1550–2650), includes lunar laser ranging (LLR) data
//!   and relativistic corrections.  Recommended for near-term applications.
//! - **DE441** — extended fit (−13 200 to +17 191), sacrificing minor accuracy
//!   improvements for broader temporal coverage.
//!
//! ## Technical scope
//!
//! - [`eval`] — Chebyshev polynomial evaluation and [`eval::SegmentDescriptor`],
//!   which maps a JD interval to a coefficient block and evaluates position /
//!   velocity for one body.
//! - [`bodies`] — generic body-chain resolution that derives Earth, Sun, and Moon
//!   from the natively integrated barycentric states (Earth–Moon Barycenter +
//!   Moon offset).
//! - [`DeData`] — trait abstracting over the per-version coefficient tables.
//! - [`DeEphemeris`] — generic zero-sized ephemeris backend generic over any
//!   `DeData` implementation; implements the shared [`Ephemeris`] trait.
//!
//! ## References
//!
//! - Standish, E. M. (1998). "JPL Planetary and Lunar Ephemerides, DE405/LE405".
//!   *JPL Interoffice Memorandum* 312.F-98-048.
//! - Folkner, W. M., Williams, J. G., Boggs, D. H., Park, R. S., & Kuchynka, P.
//!   (2014). "The Planetary and Lunar Ephemerides DE430 and DE431".
//!   *IPN Progress Report* 42-196, 1–81.
//! - Park, R. S., et al. (2021). "The JPL Planetary and Lunar Ephemerides DE440
//!   and DE441". *The Astronomical Journal* 161, 105.
//!   <https://doi.org/10.3847/1538-3881/abd414>

pub mod bodies;
pub mod eval;

#[cfg(feature = "de440")]
pub mod de440;
#[cfg(feature = "de441")]
pub mod de441;

use eval::SegmentDescriptor;

use crate::calculus::ephemeris::{AuPerDay, Ephemeris, EphemerisError};
use crate::coordinates::{
    cartesian::{Position, Velocity},
    centers::{Barycentric, Geocentric, Heliocentric},
    frames::EclipticMeanJ2000,
};
use crate::qtty::{AstronomicalUnit, Kilometer};
use crate::time::JulianDate;

use core::marker::PhantomData;

/// Trait abstracting over DE4xx coefficient data sources.
///
/// Each JPL DE version (DE440, DE441, …) implements this trait with
/// its own embedded segment descriptors. The descriptors are associated
/// constants so the compiler can fully monomorphize all evaluation paths.
pub trait DeData: 'static {
    /// Segment descriptor for the Sun (NAIF 10 → SSB).
    const SUN: SegmentDescriptor;
    /// Segment descriptor for the Earth-Moon Barycenter (NAIF 3 → SSB).
    const EMB: SegmentDescriptor;
    /// Segment descriptor for the Moon (NAIF 301 → EMB).
    const MOON: SegmentDescriptor;
}

/// Generic zero-sized ephemeris backend for any JPL DE4xx dataset.
///
/// `D` selects the coefficient data at compile time. All dispatch is
/// monomorphized away, there is no runtime cost vs. a hand-written backend.
///
/// # Example
/// ```rust,ignore
/// use siderust::calculus::jpl::DeEphemeris;
/// use siderust::calculus::de440::De440Data;
///
/// type De440Ephemeris = DeEphemeris<De440Data>;
/// ```
#[derive(Debug, Clone, Copy, Default)]
pub struct DeEphemeris<D: DeData>(PhantomData<D>);

impl<D: DeData> Ephemeris for DeEphemeris<D> {
    #[inline]
    fn try_sun_barycentric(
        jd: JulianDate,
    ) -> Result<Position<Barycentric, EclipticMeanJ2000, AstronomicalUnit>, EphemerisError> {
        bodies::try_sun_barycentric(jd, &D::SUN)
    }

    #[inline]
    fn sun_barycentric(
        jd: JulianDate,
    ) -> Position<Barycentric, EclipticMeanJ2000, AstronomicalUnit> {
        bodies::sun_barycentric(jd, &D::SUN)
    }

    #[inline]
    fn try_earth_barycentric(
        jd: JulianDate,
    ) -> Result<Position<Barycentric, EclipticMeanJ2000, AstronomicalUnit>, EphemerisError> {
        bodies::try_earth_barycentric(jd, &D::EMB, &D::MOON)
    }

    #[inline]
    fn earth_barycentric(
        jd: JulianDate,
    ) -> Position<Barycentric, EclipticMeanJ2000, AstronomicalUnit> {
        bodies::earth_barycentric(jd, &D::EMB, &D::MOON)
    }

    #[inline]
    fn try_earth_heliocentric(
        jd: JulianDate,
    ) -> Result<Position<Heliocentric, EclipticMeanJ2000, AstronomicalUnit>, EphemerisError> {
        bodies::try_earth_heliocentric(jd, &D::SUN, &D::EMB, &D::MOON)
    }

    #[inline]
    fn earth_heliocentric(
        jd: JulianDate,
    ) -> Position<Heliocentric, EclipticMeanJ2000, AstronomicalUnit> {
        bodies::earth_heliocentric(jd, &D::SUN, &D::EMB, &D::MOON)
    }

    #[inline]
    fn try_earth_barycentric_velocity(
        jd: JulianDate,
    ) -> Result<Velocity<EclipticMeanJ2000, AuPerDay>, EphemerisError> {
        bodies::try_earth_barycentric_velocity(jd, &D::EMB, &D::MOON)
    }

    #[inline]
    fn earth_barycentric_velocity(jd: JulianDate) -> Velocity<EclipticMeanJ2000, AuPerDay> {
        bodies::earth_barycentric_velocity(jd, &D::EMB, &D::MOON)
    }

    #[inline]
    fn try_moon_geocentric(
        jd: JulianDate,
    ) -> Result<Position<Geocentric, EclipticMeanJ2000, Kilometer>, EphemerisError> {
        bodies::try_moon_geocentric(jd, &D::MOON)
    }

    #[inline]
    fn moon_geocentric(jd: JulianDate) -> Position<Geocentric, EclipticMeanJ2000, Kilometer> {
        bodies::moon_geocentric(jd, &D::MOON)
    }
}
