// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Solar Altitude Window Periods
//!
//! Sun-specific routines for finding time intervals where the Sun's
//! altitude is above, below, or within a given range.
//!
//! All period-finding is delegated to [`crate::calculus::math_core::intervals`]
//! which provides scan + Brent refinement + crossing classification + interval
//! assembly.  This module supplies the altitude closure and JD↔f64 / Mjd
//! conversions.
//!
//! ## Key Points
//!
//! * A 2-hour scan step safely brackets every sunrise/sunset
//!   (the shortest daylight arc at 65° latitude is ~5 h).
//! * Below-threshold and range variants are derived at negligible cost
//!   via [`crate::time::complement_within`] / set intersection.
//!
//! ## Performance
//!
//! The 2-hour scan uses ~12 VSOP87 evaluations per day, compared to ~72
//! for the previous culmination-based approach (20-min hour-angle scan).
//! This yields a ~6× reduction in ephemeris evaluations.

use crate::bodies::solar_system::Sun;
use crate::calculus::math_core::intervals;
use crate::coordinates::centers::Geodetic;
use crate::coordinates::frames::ECEF;
use crate::time::{complement_within, JulianDate, ModifiedJulianDate, Period, MJD};
use qtty::*;

// =============================================================================
// Constants
// =============================================================================

/// Scan step for Sun altitude threshold detection (2 hours in days).
/// Sunrise/sunset events are separated by ≥5 h even at high latitudes,
/// so 2-hour steps safely detect all crossings (~12 altitude calls per day).
const SCAN_STEP: Days = Quantity::<Hour>::new(2.0).to_const::<Day>();

// =============================================================================
// Core Altitude Function
// =============================================================================

/// Computes the Sun's altitude in **radians** at a given Julian Date and observer site.
/// Positive above the horizon, negative below.
pub(crate) fn sun_altitude_rad(mjd: ModifiedJulianDate, site: &Geodetic<ECEF>) -> Quantity<Radian> {
    let jd: JulianDate = mjd.into();
    Sun::get_horizontal::<AstronomicalUnit>(jd, *site)
        .alt()
        .to::<Radian>()
}

// =============================================================================
// Main API
// =============================================================================

/// Finds day periods (Sun **above** `threshold`) inside `period`.
///
/// Uses a 2-hour scan + Brent refinement via [`math_core::intervals`].
pub(crate) fn find_day_periods(
    site: Geodetic<ECEF>,
    period: Period<MJD>,
    threshold: Degrees,
) -> Vec<Period<MJD>> {
    let thr = threshold.to::<Radian>();

    let f = |t: ModifiedJulianDate| -> Radians { sun_altitude_rad(t, &site) };

    intervals::above_threshold_periods(period, SCAN_STEP, &f, thr)
}

/// Finds night periods (Sun **below** `twilight`) inside `period`.
///
/// Complement of [`find_day_periods`] within `period`.
pub(crate) fn find_night_periods(
    site: Geodetic<ECEF>,
    period: Period<MJD>,
    twilight: Degrees,
) -> Vec<Period<MJD>> {
    let days = find_day_periods(site, period, twilight);
    complement_within(period, &days)
}

/// Finds periods where Sun altitude is within `range` `[min, max]`.
///
/// Computed as `above(min) ∩ complement(above(max))` via
/// [`math_core::intervals::in_range_periods`].
pub(crate) fn find_sun_range_periods(
    site: Geodetic<ECEF>,
    period: Period<MJD>,
    range: (Degrees, Degrees),
) -> Vec<Period<MJD>> {
    let h_min = range.0.to::<Radian>();
    let h_max = range.1.to::<Radian>();

    let f = |t: ModifiedJulianDate| -> Radians { sun_altitude_rad(t, &site) };

    intervals::in_range_periods(period, SCAN_STEP, &f, h_min, h_max)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn greenwich_site() -> Geodetic<ECEF> {
        Geodetic::<ECEF>::new(
            Degrees::new(0.0),
            Degrees::new(51.4769),
            Quantity::<Meter>::new(0.0),
        )
    }

    #[test]
    fn test_sun_altitude_basic() {
        let site = greenwich_site();
        let mjd: ModifiedJulianDate = crate::time::JulianDate::J2000.into();
        let alt = sun_altitude_rad(mjd, &site);
        assert!(
            alt > Radians::new(-std::f64::consts::FRAC_PI_2)
                && alt < Radians::new(std::f64::consts::FRAC_PI_2)
        );
    }

    #[test]
    fn test_find_night_periods() {
        use crate::calculus::solar::twilight;

        let site = greenwich_site();
        let mjd_start = ModifiedJulianDate::new(60000.0);
        let mjd_end = ModifiedJulianDate::new(60007.0);
        let period = Period::new(mjd_start, mjd_end);

        let nights = find_night_periods(site, period, twilight::ASTRONOMICAL);
        assert!(
            !nights.is_empty(),
            "Should find night periods at 51° latitude"
        );

        for night in &nights {
            assert!(
                night.duration_days() > 0.0,
                "Night duration should be positive"
            );
            assert!(
                night.duration_days() < 1.0,
                "Night should be less than 24 hours"
            );
        }
    }

    #[test]
    fn test_find_altitude_range_periods() {
        let site = greenwich_site();
        let mjd_start = ModifiedJulianDate::new(60000.0);
        let mjd_end = ModifiedJulianDate::new(60007.0);

        let period = Period::new(mjd_start, mjd_end);

        let nights =
            find_sun_range_periods(site, period, (Degrees::new(-90.0), Degrees::new(-18.0)));

        assert!(!nights.is_empty(), "Should find night periods using range");

        let nautical =
            find_sun_range_periods(site, period, (Degrees::new(-18.0), Degrees::new(-12.0)));

        assert!(
            !nautical.is_empty(),
            "Should find nautical twilight periods"
        );
    }
}
