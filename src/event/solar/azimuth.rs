// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Solar Azimuth
//!
//! ## Scientific scope
//!
//! Sun‑specific azimuth scalar function. Returns the topocentric apparent
//! azimuth of the Sun (North‑clockwise convention, range `[0, 2π)`) for a
//! geodetic observer at the given instant. The Sun position is derived
//! from `Sun::get_horizontal` (VSOP87 + nutation + aberration). No
//! atmospheric refraction is applied.
//!
//! ## Technical scope
//!
//! Single crate‑internal helper [`sun_azimuth_rad`] used by the unified
//! azimuth calculus API via the
//! [`AzimuthProvider`](crate::event::azimuth::AzimuthProvider) trait.
//!
//! ## References
//! - Meeus, J. (1998). *Astronomical Algorithms*, 2nd ed., Willmann‑Bell.

use crate::bodies::solar_system::Sun;
use crate::coordinates::centers::Geodetic;
use crate::coordinates::frames::ECEF;
use crate::qtty::*;
use crate::time::{JulianDate, ModifiedJulianDate};

/// Computes the Sun's topocentric azimuth in **radians** (North-clockwise,
/// range `[0, 2π)`) at a given Modified Julian Date and observer site.
///
/// # Arguments
///
/// * `mjd`, instant on the TT axis
/// * `site`, geodetic observer location
///
/// # Returns
///
/// `Radians` in `[0, 2π)`, North‑clockwise (no refraction).
pub(crate) fn sun_azimuth_rad(mjd: ModifiedJulianDate, site: &Geodetic<ECEF>) -> Radians {
    let jd: JulianDate = mjd.to::<crate::JD>();
    Sun::get_horizontal::<AstronomicalUnit>(jd, *site)
        .az()
        .to::<Radian>()
}
