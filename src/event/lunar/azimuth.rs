// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Lunar Azimuth
//!
//! ## Scientific scope
//!
//! Moon‑specific azimuth scalar function. Returns the topocentric
//! apparent azimuth of the Moon (North‑clockwise convention, range
//! `[0, 2π)`) for a geodetic observer at the given instant. The Moon
//! position is derived from `Moon::get_horizontal`, which applies
//! topocentric parallax (~1° at horizon), precession, and the full
//! ecliptic → equatorial → horizontal transform chain. No atmospheric
//! refraction is applied.
//!
//! ## Technical scope
//!
//! Single crate‑internal helper [`moon_azimuth_rad`] used by the unified
//! azimuth calculus API via the
//! [`AzimuthProvider`](crate::event::azimuth::AzimuthProvider) trait.
//!
//! ## References
//! - Meeus, J. (1998). *Astronomical Algorithms*, 2nd ed., Willmann‑Bell.

use crate::bodies::solar_system::Moon;
use crate::coordinates::centers::Geodetic;
use crate::coordinates::frames::ECEF;
use crate::qtty::*;
use crate::time::{JulianDate, ModifiedJulianDate};

/// Computes the Moon's topocentric azimuth in **radians** (North-clockwise,
/// range `[0, 2π)`) at a given Modified Julian Date and observer site.
///
/// Accounts for topocentric parallax (~1° at horizon), precession, and the
/// full ecliptic → equatorial → horizontal transform chain.
///
/// # Arguments
///
/// * `mjd`, instant on the TT axis
/// * `site`, geodetic observer location
///
/// # Returns
///
/// `Radians` in `[0, 2π)`, North‑clockwise (no refraction).
pub(crate) fn moon_azimuth_rad(mjd: ModifiedJulianDate, site: &Geodetic<ECEF>) -> Radians {
    let jd: JulianDate = mjd.to::<crate::JD>();
    Moon::get_horizontal::<Kilometer>(jd, *site)
        .az()
        .to::<Radian>()
}
