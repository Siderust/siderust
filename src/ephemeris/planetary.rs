// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Body-parameterised planetary heliocentric states
//!
//! ## Scientific scope
//!
//! The core [`Ephemeris`](super::Ephemeris) trait covers Sun, Earth, and Moon —
//! the bodies required by the coordinate-transform and aberration pipelines.
//! Physical orientation of a major planet additionally needs that body's
//! heliocentric ecliptic state so geocentric (or topocentric) directions can be
//! formed with the **same** compile-time backend that supplies Earth.
//!
//! ## Technical scope
//!
//! [`PlanetState<Body>`] extends [`Ephemeris`](super::Ephemeris) with a single
//! heliocentric accessor for a phantom body type. Backends implement the trait
//! per body they can evaluate; orientation entry points require
//! `E: PlanetState<Body>` and therefore cannot silently mix VSOP87 with another
//! planetary source. [`Vsop87Ephemeris`](super::Vsop87Ephemeris) implements
//! `PlanetState` for Mercury–Neptune. [`RuntimeEphemeris`](super::RuntimeEphemeris)
//! currently indexes only Sun/Earth/Moon segments and does **not** implement it.
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

/// Compile-time capability: ephemeris backend `Self` can supply heliocentric
/// ecliptic states for phantom body `Body`.
///
/// Implementors must evaluate `Body` with the **same** theory/data set used by
/// their [`Ephemeris`](super::Ephemeris) Earth methods.
pub trait PlanetState<Body>: Ephemeris {
    /// Fallible heliocentric ecliptic position of `Body` (AU).
    #[inline]
    fn try_heliocentric(
        jd: JulianDate,
    ) -> Result<Position<Heliocentric, EclipticMeanJ2000, AstronomicalUnit>, super::EphemerisError>
    {
        Ok(Self::heliocentric(jd))
    }

    /// Heliocentric ecliptic position of `Body` (AU).
    fn heliocentric(jd: JulianDate) -> Position<Heliocentric, EclipticMeanJ2000, AstronomicalUnit>;
}
