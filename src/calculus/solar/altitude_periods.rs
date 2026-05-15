// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Solar Altitude Window Periods
//!
//! ## Scientific scope
//!
//! Sun‑specific routines for finding time intervals where the topocentric
//! altitude of the Sun is above, below, or within a given band — the
//! kinematic basis for day/night classification, twilight detection, and
//! observation‑window planning. The Sun position comes from
//! `Sun::get_horizontal` (VSOP87 + nutation + aberration); refraction is
//! not applied. A 2‑hour scan step safely brackets every sunrise/sunset
//! since the shortest daylight arc at 65° latitude is ≳ 5 h.
//!
//! ## Technical scope
//!
//! Crate‑internal API: `sun_altitude_rad`, [`find_day_periods`],
//! [`find_night_periods`], [`find_sun_range_periods`]. All period‑finding
//! is delegated to [`crate::calculus::math_core::intervals`] which
//! provides scan + Brent refinement + crossing classification + interval
//! assembly. Below‑threshold and range variants are derived at negligible
//! cost via [`crate::time::complement_within`] / set intersection.
//!
//! Performance note: the 2‑hour scan uses ~12 VSOP87 evaluations per day,
//! versus ~72 for a 20‑min hour‑angle scan (~6× reduction).
//!
//! ## References
//! None.

use crate::bodies::solar_system::Sun;
use crate::calculus::math_core::intervals;
use crate::coordinates::centers::Geodetic;
use crate::coordinates::frames::ECEF;
use crate::qtty::*;
use crate::time::{complement_within, JulianDate, ModifiedJulianDate, Period};

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
///
/// # Arguments
///
/// * `mjd`, instant on the TT axis
/// * `site`, geodetic observer location
///
/// # Returns
///
/// Topocentric altitude as `Quantity<Radian>` (no refraction).
pub(crate) fn sun_altitude_rad(mjd: ModifiedJulianDate, site: &Geodetic<ECEF>) -> Quantity<Radian> {
    let jd: JulianDate = mjd.to_time().to::<crate::time::JD>();
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
///
/// # Arguments
///
/// * `site`, geodetic observer location
/// * `period`, MJD/TT search window
/// * `threshold`, altitude threshold
///
/// # Returns
///
/// Sorted, non‑overlapping `Vec<Period<ModifiedJulianDate>>` where the Sun
/// is above `threshold`.
pub(crate) fn find_day_periods(
    site: Geodetic<ECEF>,
    period: Period<ModifiedJulianDate>,
    threshold: Degrees,
) -> Vec<Period<ModifiedJulianDate>> {
    let thr = threshold.to::<Radian>();

    let f = |t: ModifiedJulianDate| -> Radians { sun_altitude_rad(t, &site) };

    intervals::above_threshold_periods(period, SCAN_STEP, &f, thr)
}

/// Finds night periods (Sun **below** `twilight`) inside `period`.
///
/// Complement of [`find_day_periods`] within `period`.
///
/// # Arguments
///
/// * `site`, geodetic observer location
/// * `period`, MJD/TT search window
/// * `twilight`, altitude threshold (e.g. `−18°` for astronomical night)
///
/// # Returns
///
/// Sorted, non‑overlapping `Vec<Period<ModifiedJulianDate>>` where the Sun
/// is at or below `twilight`.
pub(crate) fn find_night_periods(
    site: Geodetic<ECEF>,
    period: Period<ModifiedJulianDate>,
    twilight: Degrees,
) -> Vec<Period<ModifiedJulianDate>> {
    let days = find_day_periods(site, period, twilight);
    complement_within(period, &days)
}

/// Finds periods where Sun altitude is within `range` `[min, max]`.
///
/// Computed as `above(min) ∩ complement(above(max))` via
/// [`math_core::intervals::in_range_periods`].
///
/// # Arguments
///
/// * `site`, geodetic observer location
/// * `period`, MJD/TT search window
/// * `range`, `(min_altitude, max_altitude)` band
///
/// # Returns
///
/// Sorted, non‑overlapping `Vec<Period<ModifiedJulianDate>>` of intervals
/// where `min ≤ altitude(t) ≤ max`.
pub(crate) fn find_sun_range_periods(
    site: Geodetic<ECEF>,
    period: Period<ModifiedJulianDate>,
    range: (Degrees, Degrees),
) -> Vec<Period<ModifiedJulianDate>> {
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
        let mjd: ModifiedJulianDate = crate::time::JulianDate::J2000
            .to_time()
            .to::<crate::time::MJD>();
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
        let mjd_start = ModifiedJulianDate::from_raw_unchecked(qtty::Day::new(60000.0));
        let mjd_end = ModifiedJulianDate::from_raw_unchecked(qtty::Day::new(60007.0));
        let period = Period::new(mjd_start, mjd_end);

        let nights = find_night_periods(site, period, twilight::ASTRONOMICAL);
        assert!(
            !nights.is_empty(),
            "Should find night periods at 51° latitude"
        );

        for night in &nights {
            assert!(
                (night.end.raw() - night.start.raw()) > Days::new(0.0),
                "Night duration should be positive"
            );
            assert!(
                (night.end.raw() - night.start.raw()) < Days::new(1.0),
                "Night should be less than 24 hours"
            );
        }
    }

    #[test]
    fn test_find_altitude_range_periods() {
        let site = greenwich_site();
        let mjd_start = ModifiedJulianDate::from_raw_unchecked(qtty::Day::new(60000.0));
        let mjd_end = ModifiedJulianDate::from_raw_unchecked(qtty::Day::new(60007.0));

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
