// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Precise solar altitude evaluation for event search.
//!
//! Event search windows are `ModifiedJulianDate` values on the crate's TT
//! axis. This module only evaluates altitude at one typed instant; crossing
//! prediction, fallback search, and period assembly live in sibling modules.

use crate::coordinates::centers::Geodetic;
use crate::coordinates::frames::ECEF;
use crate::qtty::*;
use crate::time::{JulianDate, ModifiedJulianDate};

/// Computes the Sun's topocentric altitude in radians at an MJD/TT instant.
pub(crate) fn sun_altitude_rad(mjd: ModifiedJulianDate, site: &Geodetic<ECEF>) -> Radians {
    let jd: JulianDate = mjd.to::<crate::JD>();
    crate::bodies::solar_system::Sun::get_horizontal::<AstronomicalUnit>(jd, *site)
        .alt()
        .to::<Radian>()
}
