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

use crate::calculus::ephemeris::{
    AuPerDay, Ephemeris, EphemerisError, MajorPlanet, PlanetEphemerisError, PlanetPoint,
};
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
    /// Segment descriptor for Mercury system barycenter (NAIF 1 -> SSB).
    const MERCURY_BARYCENTER: SegmentDescriptor;
    /// Segment descriptor for Mercury center (NAIF 199 -> Mercury barycenter).
    const MERCURY_CENTER: SegmentDescriptor;
    /// Segment descriptor for Venus system barycenter (NAIF 2 -> SSB).
    const VENUS_BARYCENTER: SegmentDescriptor;
    /// Segment descriptor for Venus center (NAIF 299 -> Venus barycenter).
    const VENUS_CENTER: SegmentDescriptor;
    /// Segment descriptor for Mars system barycenter (NAIF 4 -> SSB).
    const MARS_BARYCENTER: SegmentDescriptor;
    /// Segment descriptor for Jupiter system barycenter (NAIF 5 -> SSB).
    const JUPITER_BARYCENTER: SegmentDescriptor;
    /// Segment descriptor for Saturn system barycenter (NAIF 6 -> SSB).
    const SATURN_BARYCENTER: SegmentDescriptor;
    /// Segment descriptor for Uranus system barycenter (NAIF 7 -> SSB).
    const URANUS_BARYCENTER: SegmentDescriptor;
    /// Segment descriptor for Neptune system barycenter (NAIF 8 -> SSB).
    const NEPTUNE_BARYCENTER: SegmentDescriptor;
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

impl<D: DeData> DeEphemeris<D> {
    /// Fallible embedded JPL barycentric state for a major planet point.
    ///
    /// Generic DE kernels carry system barycenters for all major planets and
    /// direct center offsets for Mercury and Venus. Outer planet center
    /// offsets live in satellite SPK kernels and require the runtime kernel
    /// set API.
    pub fn try_major_planet_barycentric(
        planet: MajorPlanet,
        point: PlanetPoint,
        jd: JulianDate,
    ) -> Result<Position<Barycentric, EclipticMeanJ2000, AstronomicalUnit>, PlanetEphemerisError>
    {
        match (planet, point) {
            (MajorPlanet::Mercury, PlanetPoint::SystemBarycenter) => {
                Ok(bodies::try_planet_barycentric(jd, &D::MERCURY_BARYCENTER)?)
            }
            (MajorPlanet::Mercury, PlanetPoint::Center) => {
                Ok(bodies::try_child_planet_barycentric(
                    jd,
                    &D::MERCURY_BARYCENTER,
                    &D::MERCURY_CENTER,
                )?)
            }
            (MajorPlanet::Venus, PlanetPoint::SystemBarycenter) => {
                Ok(bodies::try_planet_barycentric(jd, &D::VENUS_BARYCENTER)?)
            }
            (MajorPlanet::Venus, PlanetPoint::Center) => Ok(bodies::try_child_planet_barycentric(
                jd,
                &D::VENUS_BARYCENTER,
                &D::VENUS_CENTER,
            )?),
            (MajorPlanet::Mars, PlanetPoint::SystemBarycenter) => {
                Ok(bodies::try_planet_barycentric(jd, &D::MARS_BARYCENTER)?)
            }
            (MajorPlanet::Jupiter, PlanetPoint::SystemBarycenter) => {
                Ok(bodies::try_planet_barycentric(jd, &D::JUPITER_BARYCENTER)?)
            }
            (MajorPlanet::Saturn, PlanetPoint::SystemBarycenter) => {
                Ok(bodies::try_planet_barycentric(jd, &D::SATURN_BARYCENTER)?)
            }
            (MajorPlanet::Uranus, PlanetPoint::SystemBarycenter) => {
                Ok(bodies::try_planet_barycentric(jd, &D::URANUS_BARYCENTER)?)
            }
            (MajorPlanet::Neptune, PlanetPoint::SystemBarycenter) => {
                Ok(bodies::try_planet_barycentric(jd, &D::NEPTUNE_BARYCENTER)?)
            }
            _ => Err(PlanetEphemerisError::UnsupportedPoint { planet, point }),
        }
    }

    /// Fallible embedded JPL Earth-centered state for a major planet point.
    pub fn try_major_planet_geocentric(
        planet: MajorPlanet,
        point: PlanetPoint,
        jd: JulianDate,
    ) -> Result<Position<Geocentric, EclipticMeanJ2000, AstronomicalUnit>, PlanetEphemerisError>
    {
        match (planet, point) {
            (MajorPlanet::Mercury, PlanetPoint::SystemBarycenter) => Ok(
                bodies::try_planet_geocentric(jd, &D::MERCURY_BARYCENTER, &D::EMB, &D::MOON)?,
            ),
            (MajorPlanet::Mercury, PlanetPoint::Center) => Ok(bodies::try_child_planet_geocentric(
                jd,
                &D::MERCURY_BARYCENTER,
                &D::MERCURY_CENTER,
                &D::EMB,
                &D::MOON,
            )?),
            (MajorPlanet::Venus, PlanetPoint::SystemBarycenter) => Ok(
                bodies::try_planet_geocentric(jd, &D::VENUS_BARYCENTER, &D::EMB, &D::MOON)?,
            ),
            (MajorPlanet::Venus, PlanetPoint::Center) => Ok(bodies::try_child_planet_geocentric(
                jd,
                &D::VENUS_BARYCENTER,
                &D::VENUS_CENTER,
                &D::EMB,
                &D::MOON,
            )?),
            (MajorPlanet::Mars, PlanetPoint::SystemBarycenter) => Ok(
                bodies::try_planet_geocentric(jd, &D::MARS_BARYCENTER, &D::EMB, &D::MOON)?,
            ),
            (MajorPlanet::Jupiter, PlanetPoint::SystemBarycenter) => Ok(
                bodies::try_planet_geocentric(jd, &D::JUPITER_BARYCENTER, &D::EMB, &D::MOON)?,
            ),
            (MajorPlanet::Saturn, PlanetPoint::SystemBarycenter) => Ok(
                bodies::try_planet_geocentric(jd, &D::SATURN_BARYCENTER, &D::EMB, &D::MOON)?,
            ),
            (MajorPlanet::Uranus, PlanetPoint::SystemBarycenter) => Ok(
                bodies::try_planet_geocentric(jd, &D::URANUS_BARYCENTER, &D::EMB, &D::MOON)?,
            ),
            (MajorPlanet::Neptune, PlanetPoint::SystemBarycenter) => Ok(
                bodies::try_planet_geocentric(jd, &D::NEPTUNE_BARYCENTER, &D::EMB, &D::MOON)?,
            ),
            _ => Err(PlanetEphemerisError::UnsupportedPoint { planet, point }),
        }
    }
}
