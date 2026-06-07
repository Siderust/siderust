// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! VSOP87 / ELP2000-82B ephemeris backend.
//!
//! This is the default (always-available) backend. It delegates to the
//! existing `Planet::vsop87*` inherent methods and `Moon::get_geo_position`.

use super::{AuPerDay, Ephemeris};
use crate::bodies::solar_system::{Earth, Moon, Sun};
use crate::coordinates::{
    cartesian::{Position, Velocity},
    centers::{Barycentric, Geocentric, Heliocentric},
    frames::EclipticMeanJ2000,
};
use crate::qtty::{AstronomicalUnit, Kilometer};
use crate::time::JulianDate;

/// VSOP87 / ELP2000-82B ephemeris backend (always available).
///
/// Zero-sized marker type. Coefficient tables originate from
/// [`siderust-archive`](https://crates.io/crates/siderust-archive); this crate
/// only consumes the archive-backed static coefficient snapshots. No local
/// build-time generation occurs in Siderust.
#[derive(Debug, Clone, Copy, Default)]
pub struct Vsop87Ephemeris;

impl Ephemeris for Vsop87Ephemeris {
    #[inline]
    fn sun_barycentric(
        jd: JulianDate,
    ) -> Position<Barycentric, EclipticMeanJ2000, AstronomicalUnit> {
        Sun::vsop87e(jd)
    }

    #[inline]
    fn earth_barycentric(
        jd: JulianDate,
    ) -> Position<Barycentric, EclipticMeanJ2000, AstronomicalUnit> {
        Earth::vsop87e(jd)
    }

    #[inline]
    fn earth_heliocentric(
        jd: JulianDate,
    ) -> Position<Heliocentric, EclipticMeanJ2000, AstronomicalUnit> {
        Earth::vsop87a(jd)
    }

    #[inline]
    fn earth_barycentric_velocity(jd: JulianDate) -> Velocity<EclipticMeanJ2000, AuPerDay> {
        Earth::vsop87e_vel(jd)
    }

    #[inline]
    fn moon_geocentric(jd: JulianDate) -> Position<Geocentric, EclipticMeanJ2000, Kilometer> {
        Moon::get_geo_position(jd)
    }
}
