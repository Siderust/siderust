// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Stellar Azimuth
//!
//! Fixed-star azimuth scalar function, used by the unified azimuth calculus
//! API via the [`AzimuthProvider`](crate::calculus::azimuth::AzimuthProvider) trait.

use crate::coordinates::centers::Geodetic;
use crate::coordinates::frames::ECEF;
use crate::time::{JulianDate, ModifiedJulianDate};
use qtty::*;

/// Computes the topocentric azimuth of a fixed star in **radians**
/// (North-clockwise, range `[0, 2π)`) at a given Julian Date and observer site.
///
/// Uses IAU 2006 precession + IAU 2000B nutation (via NPB matrix), ERA-based
/// GAST, and the standard equatorial→horizontal formula.
///
/// # Arguments
/// * `mjd`       — Modified Julian Date (TT)
/// * `site`      — Observer's geographic location
/// * `ra_j2000`  — Right Ascension (J2000), in degrees
/// * `dec_j2000` — Declination (J2000), in degrees
pub(crate) fn fixed_star_azimuth_rad(
    mjd: ModifiedJulianDate,
    site: &Geodetic<ECEF>,
    ra_j2000: Degrees,
    dec_j2000: Degrees,
) -> Radians {
    let jd: JulianDate = mjd.into();
    crate::calculus::horizontal::star_horizontal(ra_j2000, dec_j2000, site, jd)
        .az()
        .to::<Radian>()
}
