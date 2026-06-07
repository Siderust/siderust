// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Lunar Altitude Window Periods
//!
//! ## Scientific scope
//!
//! Moon‑specific routines for finding time intervals where the
//! topocentric Moon altitude is above, below, or within a given band —
//! the kinematic basis for moonrise/moonset, dark‑sky planning, and
//! lunar visibility windows. Topocentric parallax (~1° at horizon) is
//! handled by [`Moon::get_horizontal`] inside [`moon_altitude_rad`]; no
//! atmospheric refraction is applied. Position accuracy is bounded by the
//! underlying lunar engine (ELP truncation / DE chebyshev cache).
//!
//! ## Technical scope
//!
//! Crate‑internal API: `moon_altitude_rad`, [`find_moon_above_horizon`],
//! [`find_moon_below_horizon`], [`find_moon_altitude_range`]. All
//! period‑finding uses Chebyshev-first labelled crossing discovery with
//! precise validation and local scan+Brent fallback via
//! [`crate::event::search::crossings`]. `find_and_label_crossings` remains
//! available as a scan+Brent validation baseline.
//!
//! ## References
//! None.

use crate::bodies::solar_system::Moon;
use crate::coordinates::centers::Geodetic;
use crate::coordinates::frames::ECEF;
use crate::event::altitude::search::SearchOptsV2;
use crate::event::altitude::{CrossingDirection, CrossingEvent};
use crate::event::search::crossings;
use crate::event::search::intervals;
use crate::qtty::*;
use crate::time::{complement_within, Interval, JulianDate, ModifiedJulianDate};

use super::moon_cache::MoonAltitudeContext;

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
pub(crate) fn moon_altitude_rad(
    mjd: ModifiedJulianDate,
    site: &Geodetic<ECEF>,
) -> Quantity<Radian> {
    let jd: JulianDate = mjd.to::<crate::JD>();
    Moon::get_horizontal::<Kilometer>(jd, *site)
        .alt()
        .to::<Radian>()
}

// =============================================================================
// Main API
// =============================================================================

/// Finds periods when the Moon is above the given altitude threshold.
///
/// Uses Chebyshev-first crossing discovery with scan+Brent fallback.
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
///
/// # Returns
///
/// Sorted, non‑overlapping `Vec<Interval<ModifiedJulianDate>>` where the
/// Moon altitude is at or above `threshold`.
pub(crate) fn find_moon_above_horizon_with_search_opts(
    site: Geodetic<ECEF>,
    period: Interval<ModifiedJulianDate>,
    threshold: Degrees,
    opts: SearchOptsV2,
) -> Vec<Interval<ModifiedJulianDate>> {
    let thr = threshold.to::<Radian>();

    // Build Chebyshev + nutation caches for the period
    let ctx = MoonAltitudeContext::new(period.start, period.end, site);

    let (labeled, start_above) = find_moon_labeled_crossings_with_context(&ctx, period, thr, opts);
    intervals::build_above_periods_directed(&labeled, period, start_above)
}

/// Finds periods when the Moon is below the given altitude threshold.
///
/// Complement of [`find_moon_above_horizon`] within `period`.
///
/// # Arguments
///
/// * `site`, geodetic observer location
/// * `period`, MJD/TT search window
/// * `threshold`, altitude upper bound (e.g. `−0.5°` for "Moon down")
///
/// # Example
/// ```ignore
/// let moonless_periods = find_moon_below_horizon(site, period, Degrees::new(-0.5));
/// ```
///
/// # Returns
///
/// Sorted, non‑overlapping `Vec<Interval<ModifiedJulianDate>>` where the
/// Moon altitude is at or below `threshold`.
pub(crate) fn find_moon_below_horizon_with_search_opts(
    site: Geodetic<ECEF>,
    period: Interval<ModifiedJulianDate>,
    threshold: Degrees,
    opts: SearchOptsV2,
) -> Vec<Interval<ModifiedJulianDate>> {
    let above = find_moon_above_horizon_with_search_opts(site, period, threshold, opts);
    complement_within(period, &above)
}

/// Finds periods when Moon altitude is within a range `[min, max]`.
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
/// # Example
/// ```ignore
/// let low_moon = find_moon_altitude_range(site, period, (Degrees::new(0.0), Degrees::new(30.0)));
/// ```
///
/// # Returns
///
/// Sorted, non‑overlapping `Vec<Interval<ModifiedJulianDate>>` of intervals
/// where `min ≤ altitude(t) ≤ max`.
pub(crate) fn find_moon_altitude_range_with_search_opts(
    site: Geodetic<ECEF>,
    period: Interval<ModifiedJulianDate>,
    range: (Degrees, Degrees),
    opts: SearchOptsV2,
) -> Vec<Interval<ModifiedJulianDate>> {
    let h_min = range.0.to::<Radian>();
    let h_max = range.1.to::<Radian>();

    // Build Chebyshev + nutation caches for the period
    let ctx = MoonAltitudeContext::new(period.start, period.end, site);

    let (min_crossings, start_above_min) =
        find_moon_labeled_crossings_with_context(&ctx, period, h_min, opts);
    let above_min =
        intervals::build_above_periods_directed(&min_crossings, period, start_above_min);

    let (max_crossings, start_above_max) =
        find_moon_labeled_crossings_with_context(&ctx, period, h_max, opts);
    let above_max =
        intervals::build_above_periods_directed(&max_crossings, period, start_above_max);
    let below_max = intervals::complement(period, &above_max);
    intervals::intersect(&above_min, &below_max)
}

pub(crate) fn find_moon_crossings_with_search_opts(
    site: Geodetic<ECEF>,
    period: Interval<ModifiedJulianDate>,
    threshold: Degrees,
    opts: SearchOptsV2,
) -> Vec<CrossingEvent> {
    let thr = threshold.to::<Radian>();
    let ctx = MoonAltitudeContext::new(period.start, period.end, site);
    let (labeled, _) = find_moon_labeled_crossings_with_context(&ctx, period, thr, opts);
    labeled
        .iter()
        .map(|crossing| CrossingEvent {
            mjd: crossing.t,
            direction: if crossing.direction > 0 {
                CrossingDirection::Rising
            } else {
                CrossingDirection::Setting
            },
        })
        .collect()
}

fn find_moon_labeled_crossings_with_context(
    ctx: &MoonAltitudeContext,
    period: Interval<ModifiedJulianDate>,
    threshold: Radians,
    opts: SearchOptsV2,
) -> (Vec<intervals::LabeledCrossing>, bool) {
    let signal = |t: ModifiedJulianDate| -> f64 { ctx.altitude_rad(t).sin() };
    let (labeled, start_above, _) =
        crossings::find_labelled_crossings(period, SCAN_STEP, &signal, threshold.sin(), opts);
    (labeled, start_above)
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::catalogs::observatories::ROQUE_DE_LOS_MUCHACHOS;

    fn greenwich_site() -> Geodetic<ECEF> {
        Geodetic::<ECEF>::new(Degrees::new(0.0), Degrees::new(51.4769), Meters::new(0.0))
    }

    #[test]
    fn test_moon_altitude_basic() {
        let site = greenwich_site();
        let mjd: ModifiedJulianDate = crate::J2000.to::<crate::MJD>();
        let alt = moon_altitude_rad(mjd, &site);
        assert!(
            alt > -std::f64::consts::FRAC_PI_2 * RAD && alt < std::f64::consts::FRAC_PI_2 * RAD
        );
    }

    #[test]
    fn test_find_moon_above_horizon() {
        let site = greenwich_site();
        let mjd_start = crate::time::ModifiedJulianDate::new(60000.0);
        let mjd_end = crate::time::ModifiedJulianDate::new(60007.0);
        let period = Interval::new(mjd_start, mjd_end);

        let periods = find_moon_above_horizon_with_search_opts(
            site,
            period,
            Degrees::new(0.0),
            SearchOptsV2::default(),
        );
        assert!(
            !periods.is_empty(),
            "Should find moon-up periods over 7 days"
        );

        for p in &periods {
            assert!(
                p.length() > Days::new(0.0),
                "Period duration should be positive"
            );
        }
    }

    #[test]
    fn test_find_moon_below_horizon() {
        let site = ROQUE_DE_LOS_MUCHACHOS.geodetic();
        let mjd_start = crate::time::ModifiedJulianDate::new(60676.0);
        let mjd_end = crate::time::ModifiedJulianDate::new(60683.0);
        let period = Interval::new(mjd_start, mjd_end);

        let periods = find_moon_below_horizon_with_search_opts(
            site,
            period,
            Degrees::new(-0.5),
            SearchOptsV2::default(),
        );
        assert!(!periods.is_empty(), "Should find moon-down periods");
    }

    fn scan_baseline_above_horizon(
        ctx: &MoonAltitudeContext,
        period: Interval<ModifiedJulianDate>,
        threshold: Degrees,
    ) -> Vec<Interval<ModifiedJulianDate>> {
        let f = |t: ModifiedJulianDate| -> Radians { ctx.altitude_rad(t) };
        let (labeled, start_above) = crate::event::lunar::moon_cache::find_and_label_crossings(
            period,
            SCAN_STEP,
            &f,
            threshold.to::<Radian>(),
        );
        intervals::build_above_periods_directed(&labeled, period, start_above)
    }

    #[test]
    fn chebyshev_moon_above_horizon_matches_scan_baseline() {
        // Verify that math_core-based and scan-based algorithms produce similar results
        let site = greenwich_site();
        let mjd_start = crate::time::ModifiedJulianDate::new(60000.0);
        let mjd_end = crate::time::ModifiedJulianDate::new(60003.0);
        let period = Interval::new(mjd_start, mjd_end);

        let main_result = find_moon_above_horizon_with_search_opts(
            site,
            period,
            Degrees::new(0.0),
            SearchOptsV2::default(),
        );
        let ctx = MoonAltitudeContext::new(period.start, period.end, site);
        let scan_result = scan_baseline_above_horizon(&ctx, period, Degrees::new(0.0));

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
                (m.start.raw() - s.start.raw()).abs() < tolerance,
                "Start times should match within 1 minute"
            );
            assert!(
                (m.end.raw() - s.end.raw()).abs() < tolerance,
                "End times should match within 1 minute"
            );
        }
    }
}
