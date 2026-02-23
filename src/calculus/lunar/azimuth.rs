// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Lunar Azimuth
//!
//! Moon-specific azimuth scalar function, used by the unified azimuth
//! calculus API via the [`AzimuthProvider`](crate::calculus::azimuth::AzimuthProvider) trait.

use crate::bodies::solar_system::Moon;
use crate::coordinates::centers::Geodetic;
use crate::coordinates::frames::ECEF;
use crate::time::{JulianDate, ModifiedJulianDate};
use qtty::*;

/// Computes the Moon's topocentric azimuth in **radians** (North-clockwise,
/// range `[0, 2π)`) at a given Modified Julian Date and observer site.
///
/// Accounts for topocentric parallax (~1° at horizon), precession, and the
/// full ecliptic → equatorial → horizontal transform chain.
pub(crate) fn moon_azimuth_rad(mjd: ModifiedJulianDate, site: &Geodetic<ECEF>) -> Radians {
    let jd: JulianDate = mjd.into();
    Moon::get_horizontal::<Kilometer>(jd, *site)
        .az()
        .to::<Radian>()
}
