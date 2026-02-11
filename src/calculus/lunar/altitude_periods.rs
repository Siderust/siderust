// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Lunar Altitude Window Periods
//!
//! Moon-specific routines for finding time intervals where the Moon's
//! altitude is above, below, or within a given range.
//!
//! All period-finding is delegated to [`crate::calculus::math_core::intervals`]
//! which provides scan + Brent refinement + crossing classification + interval
//! assembly.  This module supplies the altitude closure and JD↔f64 / Mjd
//! conversions.
//!
//! ## Key Points
//!
//! * **Topocentric parallax** (~1° at horizon) is handled by
//!   [`Moon::get_horizontal`] inside [`moon_altitude_rad`].
//! * A 2-hour scan step safely brackets every moonrise/moonset
//!   (the shortest above-horizon arc is ~4 h).
//! * Below-threshold and range variants are derived at negligible cost
//!   via [`crate::time::complement_within`] / set intersection.

use crate::bodies::solar_system::Moon;
use crate::calculus::math_core::intervals;
use crate::coordinates::centers::ObserverSite;
use crate::time::{complement_within, JulianDate, ModifiedJulianDate, Period, MJD};
use qtty::*;

use super::moon_cache::{find_and_label_crossings, MoonAltitudeContext};

// =============================================================================
// Constants
// =============================================================================

/// Scan step for Moon altitude threshold detection (2 hours in days).
/// Moon rises/sets are separated by ~12+ hours, so 2-hour steps
/// safely detect all horizon crossings (~12 altitude calls per day).
const SCAN_STEP: Days = Quantity::new(2.0 / 24.0);

// =============================================================================
// Core Altitude Function
// =============================================================================

/// Computes the Moon's **topocentric** altitude at a given Julian Date.
///
/// Accounts for topocentric parallax (~1° at horizon), precession, and the
/// full ecliptic → equatorial → horizontal transform chain.
///
/// # Arguments
/// * `jd` - Julian Date (TT) for the calculation
/// * `site` - Observer's geographic location
///
/// # Returns
/// Altitude as `Quantity<Radian>` (positive above horizon, negative below)
pub(crate) fn moon_altitude_rad(mjd: ModifiedJulianDate, site: &ObserverSite) -> Quantity<Radian> {
    let jd: JulianDate = mjd.into();
    Moon::get_horizontal::<Kilometer>(jd, *site)
        .alt()
        .to::<Radian>()
}

// =============================================================================
// Main API
// =============================================================================

/// Finds periods when the Moon is above the given altitude threshold.
///
/// Uses a 2-hour scan + Brent refinement via [`math_core::intervals`].
///
/// # Arguments
/// * `site` - Observer's geographic location
/// * `period` - Time period to search
/// * `threshold` - Minimum altitude (e.g., 0° for horizon)
///
/// # Example
/// ```ignore
/// let moonrise_periods = find_moon_above_horizon(site, period, Degrees::new(0.0));
/// ```
pub(crate) fn find_moon_above_horizon(
    site: ObserverSite,
    period: Period<MJD>,
    threshold: Degrees,
) -> Vec<Period<MJD>> {
    let thr = threshold.to::<Radian>();

    // Build Chebyshev + nutation caches for the period
    let ctx = MoonAltitudeContext::new(period.start, period.end, site);

    let f = |t: ModifiedJulianDate| -> Radians { ctx.altitude_rad(t) };

    // Use find_and_label_crossings to avoid probe evaluations
    let (labeled, start_above) = find_and_label_crossings(period, SCAN_STEP, &f, thr);
    intervals::build_above_periods(&labeled, period, start_above, &f, thr)
}

/// Finds periods when the Moon is below the given altitude threshold.
///
/// Complement of [`find_moon_above_horizon`] within `period`.
///
/// # Example
/// ```ignore
/// let moonless_periods = find_moon_below_horizon(site, period, Degrees::new(-0.5));
/// ```
pub(crate) fn find_moon_below_horizon(
    site: ObserverSite,
    period: Period<MJD>,
    threshold: Degrees,
) -> Vec<Period<MJD>> {
    let above = find_moon_above_horizon(site, period, threshold);
    complement_within(period, &above)
}

/// Finds periods when Moon altitude is within a range `[min, max]`.
///
/// Computed as `above(min) ∩ complement(above(max))` via
/// [`math_core::intervals::in_range_periods`].
///
/// # Example
/// ```ignore
/// let low_moon = find_moon_altitude_range(site, period, (Degrees::new(0.0), Degrees::new(30.0)));
/// ```
pub(crate) fn find_moon_altitude_range(
    site: ObserverSite,
    period: Period<MJD>,
    range: (Degrees, Degrees),
) -> Vec<Period<MJD>> {
    let h_min = range.0.to::<Radian>();
    let h_max = range.1.to::<Radian>();

    // Build Chebyshev + nutation caches for the period
    let ctx = MoonAltitudeContext::new(period.start, period.end, site);

    let f = |t: ModifiedJulianDate| -> Radians { ctx.altitude_rad(t) };

    intervals::in_range_periods(period, SCAN_STEP, &f, h_min, h_max)
}

// =============================================================================
// Scan-based variants (10-minute step, for comparison / validation)
// =============================================================================

#[cfg(test)]
/// Scan step for 10-minute scan variants (days).
const SCAN_STEP_10MIN: Days = Quantity::new(10.0 / 1440.0);

#[cfg(test)]
/// Finds periods using the generic scan-based algorithm (above threshold).
///
/// Uses 10-minute steps via the generic scan engine. Slower but useful
/// for comparison / validation.
fn find_moon_above_horizon_scan(
    site: ObserverSite,
    period: Period<MJD>,
    threshold: Degrees,
) -> Vec<Period<MJD>> {
    let thr = threshold.to::<Radian>();

    let f = |t: ModifiedJulianDate| -> Radians { moon_altitude_rad(t, &site) };

    intervals::above_threshold_periods(period, SCAN_STEP_10MIN, &f, thr)
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{observatories::ROQUE_DE_LOS_MUCHACHOS, time::JulianDate};

    fn greenwich_site() -> ObserverSite {
        ObserverSite::new(Degrees::new(0.0), Degrees::new(51.4769), Meters::new(0.0))
    }

    #[test]
    fn test_moon_altitude_basic() {
        let site = greenwich_site();
        let mjd: ModifiedJulianDate = JulianDate::J2000.into();
        let alt = moon_altitude_rad(mjd, &site);
        assert!(
            alt > -std::f64::consts::FRAC_PI_2 * RAD && alt < std::f64::consts::FRAC_PI_2 * RAD
        );
    }

    #[test]
    fn test_find_moon_above_horizon() {
        let site = greenwich_site();
        let mjd_start = ModifiedJulianDate::new(60000.0);
        let mjd_end = ModifiedJulianDate::new(60007.0);
        let period = Period::new(mjd_start, mjd_end);

        let periods = find_moon_above_horizon(site, period, Degrees::new(0.0));
        assert!(
            !periods.is_empty(),
            "Should find moon-up periods over 7 days"
        );

        for p in &periods {
            assert!(
                p.duration_days() > 0.0,
                "Period duration should be positive"
            );
        }
    }

    #[test]
    fn test_find_moon_below_horizon() {
        let site = ObserverSite::from_geographic(&ROQUE_DE_LOS_MUCHACHOS);
        let mjd_start = ModifiedJulianDate::new(60676.0);
        let mjd_end = ModifiedJulianDate::new(60683.0);
        let period = Period::new(mjd_start, mjd_end);

        let periods = find_moon_below_horizon(site, period, Degrees::new(-0.5));
        assert!(!periods.is_empty(), "Should find moon-down periods");
    }

    #[test]
    fn test_above_vs_scan_consistency() {
        // Verify that math_core-based and scan-based algorithms produce similar results
        let site = greenwich_site();
        let mjd_start = ModifiedJulianDate::new(60000.0);
        let mjd_end = ModifiedJulianDate::new(60003.0);
        let period = Period::new(mjd_start, mjd_end);

        let main_result = find_moon_above_horizon(site, period, Degrees::new(0.0));
        let scan_result = find_moon_above_horizon_scan(site, period, Degrees::new(0.0));

        assert!(!main_result.is_empty());
        assert!(!scan_result.is_empty());

        assert_eq!(
            main_result.len(),
            scan_result.len(),
            "Should find same number of periods"
        );

        // Boundaries should match within tolerance (~1 minute)
        let tolerance = Days::new(1.0 / 1440.0);
        for (m, s) in main_result.iter().zip(scan_result.iter()) {
            assert!(
                (m.start - s.start).abs() < tolerance,
                "Start times should match within 1 minute"
            );
            assert!(
                (m.end - s.end).abs() < tolerance,
                "End times should match within 1 minute"
            );
        }
    }
}
