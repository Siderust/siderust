// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Solar Azimuth
//!
//! Sun-specific azimuth scalar function, used by the unified azimuth
//! calculus API via the [`AzimuthProvider`](crate::calculus::azimuth::AzimuthProvider) trait.

use crate::bodies::solar_system::Sun;
use crate::coordinates::centers::Geodetic;
use crate::coordinates::frames::ECEF;
use crate::time::{JulianDate, ModifiedJulianDate};
use qtty::*;

/// Computes the Sun's topocentric azimuth in **radians** (North-clockwise,
/// range `[0, 2π)`) at a given Modified Julian Date and observer site.
pub(crate) fn sun_azimuth_rad(mjd: ModifiedJulianDate, site: &Geodetic<ECEF>) -> Radians {
    let jd: JulianDate = mjd.into();
    Sun::get_horizontal::<AstronomicalUnit>(jd, *site)
        .az()
        .to::<Radian>()
}
