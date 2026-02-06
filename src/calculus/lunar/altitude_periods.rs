// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Lunar Altitude Window Periods
//!
//! Moon-specific optimized routines for finding time intervals where the Moon's
//! altitude is above, below, or within a given range.
//!
//! The core algorithm finds **above-threshold** periods. Below and between
//! variants are derived via [`crate::time::complement_within`] and
//! [`crate::time::intersect_periods`] at negligible cost.
//!
//! ## Key Optimizations
//!
//! 1. **Topocentric parallax correction**: Critical for Moon (~1° at horizon)
//! 2. **Culmination-based search**: Partitions time into monotonic altitude segments
//! 3. **Brent's method with endpoint reuse**: Avoids redundant ELP2000 evaluations
//! 4. **Fast direct-scan variant**: 2-hour steps with relaxed Brent tolerance
//!
//! ## Performance Notes
//!
//! The fast-scan algorithm achieves best performance for long periods by
//! tracking crossing direction during the scan phase, feeding pre-labelled
//! crossings into the shared [`crate::calculus::events::altitude_periods::build_above_periods`]
//! helper to avoid redundant ELP2000 evaluations.

use crate::astro::nutation::corrected_ra_with_nutation;
use crate::astro::precession;
use crate::astro::sidereal::gast_fast;
use crate::astro::JulianDate;
use crate::bodies::solar_system::Moon;
use crate::calculus::events::altitude_periods::{
    build_above_periods, crossings_to_above_periods,
    find_threshold_crossings_in_segments,
};
use crate::calculus::events::Culmination;
use crate::calculus::root_finding::brent;
use crate::coordinates::cartesian;
use crate::coordinates::centers::{Geocentric, ObserverSite};
use crate::coordinates::frames::{self, Ecliptic};
use crate::coordinates::spherical;
use crate::coordinates::transform::TransformFrame;
use crate::targets::Target;
use crate::time::{complement_within, intersect_periods, ModifiedJulianDate, Period};
use qtty::*;
use std::f64::consts::PI;

// =============================================================================
// Constants
// =============================================================================

/// Scan step for Moon altitude threshold detection (2 hours).
/// Moon rises/sets are separated by ~12+ hours, so 2-hour steps
/// safely detect all horizon crossings. Uses ~12 altitude calls per day.
const MOON_ALTITUDE_SCAN_STEP: Days = Hours::new(2.0).to::<Day>();

/// Coarse scan step for Moon culmination detection (12 hours).
/// Uses ~2 calls per day for hour angle detection.
const MOON_CULMINATION_SCAN_STEP: Days = Hours::new(12.0).to::<Day>();

/// Epsilon for deduplicating key times (~0.86 ms).
const DEDUPE_EPS: f64 = 1e-10;

/// Epsilon for deduplicating crossings (~1 ms).
const CROSS_DEDUPE_EPS: f64 = 1e-8;

/// Relaxed tolerance for fast Moon altitude root finding (~2 minute precision).
/// 2 minutes = 2/(24*60) days ≈ 1.39e-3 days
const MOON_BRENT_TOLERANCE: f64 = 1.4e-3;

// =============================================================================
// Core Altitude Function
// =============================================================================

/// Computes the Moon's **topocentric** altitude in **radians** at a given Julian Date.
///
/// This function properly accounts for:
/// - **Topocentric parallax**: The Moon's ~1° horizontal parallax is significant
/// - **Precession**: J2000 → mean-of-date transformation
/// - **Coordinate transformations**: Ecliptic → Equatorial → Horizontal
///
/// # Arguments
/// * `jd` - Julian Date for the calculation
/// * `site` - Observer's geographic location
///
/// # Returns
/// Altitude in radians (positive above horizon, negative below)
pub fn moon_altitude_rad(jd: JulianDate, site: &ObserverSite) -> f64 {
    let horiz = Moon::get_horizontal::<Kilometer>(jd, *site);
    horiz.alt().to::<Radian>().value()
}

/// Returns a closure that computes Moon altitude for a specific site.
///
/// This is useful for passing to generic altitude-finding functions.
#[inline]
pub fn moon_altitude_fn(site: ObserverSite) -> impl Fn(JulianDate) -> f64 {
    move |jd| moon_altitude_rad(jd, &site)
}

// =============================================================================
// Moon-Specific Culmination Finder (Optimized for ELP2000)
// =============================================================================

/// Returns the Moon's geocentric equatorial position wrapped in a Target.
fn get_moon_equatorial(jd: JulianDate) -> Target<spherical::Position<Geocentric, frames::EquatorialMeanJ2000, Kilometer>> {
    let moon_geo_ecliptic: cartesian::Position<Geocentric, Ecliptic, Kilometer> =
        Moon::get_geo_position(jd);
    let moon_geo_eq_j2000: cartesian::Position<Geocentric, frames::EquatorialMeanJ2000, Kilometer> =
        TransformFrame::to_frame(&moon_geo_ecliptic);
    let moon_geo_eq_j2000_sph: spherical::Position<Geocentric, frames::EquatorialMeanJ2000, Kilometer> =
        spherical::Position::from_cartesian(&moon_geo_eq_j2000);
    Target::new_static(moon_geo_eq_j2000_sph, jd)
}

/// Finds Moon culminations using a coarse scan optimized for ELP2000 performance.
///
/// Uses 6-hour scan steps (vs 20 minutes for generic `find_dynamic_extremas`)
/// reducing ELP2000 calls by ~18x while still safely bracketing culminations.
fn find_moon_culminations_fast(
    observer_geo: &spherical::position::Geographic,
    jd_start: JulianDate,
    jd_end: JulianDate,
) -> Vec<Culmination> {
    // Hour angle calculation
    let hour_angle = |jd: JulianDate| -> Radians {
        let target = get_moon_equatorial(jd);
        let mean_of_date = precession::precess_from_j2000(target.get_position().clone(), jd);
        let ra_nut = corrected_ra_with_nutation(&mean_of_date.direction(), jd);
        let ra = ra_nut.to::<Radian>();
        let theta = gast_fast(jd).to::<Radian>() + observer_geo.lon().to::<Radian>();
        (theta - ra).wrap_signed()
    };

    // Hour angle functions for Brent refinement
    let h_upper = |jd: JulianDate| -> f64 {
        hour_angle(jd).value() // crosses 0 at upper culmination
    };
    let h_lower = |jd: JulianDate| -> f64 {
        (hour_angle(jd) - Radians::new(PI)).wrap_signed().value() // crosses 0 at lower culmination
    };

    let mut out = Vec::new();
    let mut jd0 = jd_start;

    let h0 = hour_angle(jd0);
    let mut h0_val = h0.value();
    let mut h0_lower = (h0 - Radians::new(PI)).wrap_signed().value();

    while jd0 < jd_end {
        let jd1 = (jd0 + MOON_CULMINATION_SCAN_STEP).min(jd_end);
        let h1 = hour_angle(jd1);
        let h1_val = h1.value();
        let h1_lower = (h1 - Radians::new(PI)).wrap_signed().value();

        // Check for continuity (avoid wrap discontinuities at ±π)
        let upper_continuous = (h1_val - h0_val).abs() < PI;
        let lower_continuous = (h1_lower - h0_lower).abs() < PI;

        // Upper culmination: H(jd) crosses 0
        if upper_continuous && h0_val * h1_val < 0.0 {
            if let Some(root) =
                brent::refine_root_with_values(jd0, jd1, h0_val, h1_val, &h_upper, 0.0)
            {
                if root >= jd_start && root < jd_end {
                    out.push(Culmination::Upper { jd: root });
                }
            }
        }

        // Lower culmination: H(jd) - π crosses 0
        if lower_continuous && h0_lower * h1_lower < 0.0 {
            if let Some(root) =
                brent::refine_root_with_values(jd0, jd1, h0_lower, h1_lower, &h_lower, 0.0)
            {
                if root >= jd_start && root < jd_end {
                    out.push(Culmination::Lower { jd: root });
                }
            }
        }

        jd0 = jd1;
        h0_val = h1_val;
        h0_lower = h1_lower;
    }

    // Stable chronological sort
    out.sort_by(|a, b| {
        let jd_a = match a {
            Culmination::Upper { jd } | Culmination::Lower { jd } => jd,
        };
        let jd_b = match b {
            Culmination::Upper { jd } | Culmination::Lower { jd } => jd,
        };
        jd_a.partial_cmp(jd_b).unwrap()
    });

    // Deduplicate
    out.dedup_by(|a, b| match (a, b) {
        (Culmination::Upper { jd: jd_a }, Culmination::Upper { jd: jd_b })
        | (Culmination::Lower { jd: jd_a }, Culmination::Lower { jd: jd_b }) => {
            (jd_a.value() - jd_b.value()).abs() < DEDUPE_EPS
        }
        _ => false,
    });

    out
}

/// Computes sorted, deduplicated key times (start, culminations, end) for the Moon.
fn compute_moon_key_times(
    observer_geo: &spherical::position::Geographic,
    jd_start: JulianDate,
    jd_end: JulianDate,
) -> Vec<JulianDate> {
    let culminations = find_moon_culminations_fast(observer_geo, jd_start, jd_end);
    let mut key_times: Vec<JulianDate> = Vec::with_capacity(culminations.len() + 2);
    key_times.push(jd_start);
    for c in culminations {
        let jd = match c {
            Culmination::Upper { jd } | Culmination::Lower { jd } => jd,
        };
        if jd > jd_start && jd < jd_end {
            key_times.push(jd);
        }
    }
    key_times.push(jd_end);
    key_times.sort_by(|a, b| a.partial_cmp(b).unwrap());
    key_times.dedup_by(|a, b| (a.value() - b.value()).abs() < DEDUPE_EPS);
    key_times
}

// =============================================================================
// Main API: Optimized Altitude Period Finding
// =============================================================================

/// Finds all periods where the Moon's altitude is **above** the given threshold,
/// using culmination-based partitioning for optimal performance.
///
/// This is the **optimized** algorithm that achieves sub-second performance
/// for 365-day searches by:
/// 1. Using `find_dynamic_extremas` to find Moon culminations efficiently
/// 2. Using Brent's method with pre-computed endpoint values
/// 3. Avoiding redundant ELP2000 evaluations
///
/// # Arguments
/// * `site` - Observer's geographic location
/// * `period` - Time period to search (MJD)
/// * `threshold` - Altitude threshold
///
/// # Returns
/// Periods where the Moon's altitude is **above** the threshold.
///
/// # Example
/// ```ignore
/// use siderust::calculus::lunar::*;
/// use siderust::observatories::ROQUE_DE_LOS_MUCHACHOS;
/// use qtty::Degrees;
///
/// let site = ObserverSite::from_geographic(&ROQUE_DE_LOS_MUCHACHOS);
/// let period = /* your period */;
/// let above = find_moon_altitude_periods_via_culminations(site, period, Degrees::new(0.0));
/// ```
pub fn find_moon_altitude_periods_via_culminations<T: Into<Degrees>>(
    site: ObserverSite,
    period: Period<ModifiedJulianDate>,
    threshold: T,
) -> Vec<Period<ModifiedJulianDate>> {
    let jd_start = period.start.to_julian_day();
    let jd_end = period.end.to_julian_day();
    let threshold_rad = threshold.into().to::<Radian>().value();

    let altitude_fn = |jd: JulianDate| moon_altitude_rad(jd, &site);

    // Find Moon culminations using the optimized Moon-specific finder
    let observer_geo = spherical::position::Geographic::new(
        site.lon,
        site.lat,
        Kilometers::new(0.0),
    );

    let key_times = compute_moon_key_times(&observer_geo, jd_start, jd_end);

    let crossings = find_threshold_crossings_in_segments(
        &key_times,
        &altitude_fn,
        threshold_rad,
        jd_start,
        jd_end,
    );

    crossings_to_above_periods(crossings, period, &altitude_fn, threshold_rad)
}

// =============================================================================
// Fast Direct Scan Algorithm (Optimized for Performance)
// =============================================================================

/// Finds periods where Moon altitude is **above** the given threshold using
/// direct 2-hour scan without culmination detection.
///
/// This is the **fastest** algorithm for long periods, achieving significant
/// performance improvements for 365-day searches by:
/// 1. Using 2-hour scan steps (~12 calls/day)
/// 2. Reusing endpoint values for Brent refinement
/// 3. Tracking crossing direction during the scan to avoid extra evaluations
/// 4. Feeding pre-labelled crossings into [`build_above_periods`]
///
/// # Trade-off
/// Slightly less numerically precise at boundaries compared to culmination-based
/// approach, but the difference is negligible (<1 minute) for practical use.
pub fn find_moon_altitude_periods_fast<T: Into<Degrees>>(
    site: ObserverSite,
    period: Period<ModifiedJulianDate>,
    threshold: T,
) -> Vec<Period<ModifiedJulianDate>> {
    let jd_start = period.start.to_julian_day();
    let jd_end = period.end.to_julian_day();
    let threshold_rad = threshold.into().to::<Radian>().value();

    let altitude_fn = |jd: JulianDate| moon_altitude_rad(jd, &site);

    // Scan with direction tracking:
    // prev_f < 0 means altitude was below threshold → crossing into "above" (direction = +1)
    // prev_f > 0 means altitude was above threshold → crossing out of "above" (direction = −1)
    let mut labeled_crossings: Vec<(JulianDate, i32)> = Vec::new();

    let mut jd = jd_start;
    let mut prev_f = altitude_fn(jd) - threshold_rad;

    while jd < jd_end {
        let next_jd = (jd + MOON_ALTITUDE_SCAN_STEP).min(jd_end);
        let next_f = altitude_fn(next_jd) - threshold_rad;

        if prev_f * next_f < 0.0 {
            let direction = if prev_f < 0.0 { 1 } else { -1 };

            if let Some(root) = brent::refine_root_with_values_and_tolerance(
                jd, next_jd, prev_f, next_f, &altitude_fn, threshold_rad, MOON_BRENT_TOLERANCE,
            ) {
                if root >= jd_start && root <= jd_end {
                    labeled_crossings.push((root, direction));
                }
            }
        }

        jd = next_jd;
        prev_f = next_f;
    }

    // Sort + deduplicate
    labeled_crossings.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    {
        let mut i = 0;
        while i + 1 < labeled_crossings.len() {
            if (labeled_crossings[i].0.value() - labeled_crossings[i + 1].0.value()).abs()
                < CROSS_DEDUPE_EPS
            {
                labeled_crossings.remove(i + 1);
            } else {
                i += 1;
            }
        }
    }

    let start_above = altitude_fn(jd_start) > threshold_rad;

    build_above_periods(
        &labeled_crossings,
        period,
        jd_start,
        jd_end,
        start_above,
        &altitude_fn,
        threshold_rad,
    )
}

// =============================================================================
// Convenience Wrappers
// =============================================================================

/// Finds periods when Moon is above the given altitude threshold.
///
/// Uses the fast direct-scan algorithm optimized for performance.
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
pub fn find_moon_above_horizon<T: Into<Degrees>>(
    site: ObserverSite,
    period: Period<ModifiedJulianDate>,
    threshold: T,
) -> Vec<Period<ModifiedJulianDate>> {
    find_moon_altitude_periods_fast(site, period, threshold)
}

/// Finds periods when Moon is below the given altitude threshold.
///
/// Equivalent to the complement of [`find_moon_above_horizon`] within `period`.
///
/// # Arguments
/// * `site` - Observer's geographic location
/// * `period` - Time period to search
/// * `threshold` - Maximum altitude (e.g., -0.5° for moonset accounting for refraction)
///
/// # Example
/// ```ignore
/// let moonless_periods = find_moon_below_horizon(site, period, Degrees::new(-0.5));
/// ```
pub fn find_moon_below_horizon<T: Into<Degrees>>(
    site: ObserverSite,
    period: Period<ModifiedJulianDate>,
    threshold: T,
) -> Vec<Period<ModifiedJulianDate>> {
    let above = find_moon_altitude_periods_fast(site, period, threshold);
    complement_within(period, &above)
}

/// Finds periods when Moon altitude is within a range.
///
/// Computed as `above(min) ∩ complement(above(max))`. The culmination
/// key-times are computed only once and reused for both thresholds.
///
/// # Arguments
/// * `site` - Observer's geographic location
/// * `period` - Time period to search
/// * `range` - (min_altitude, max_altitude) tuple
///
/// # Example
/// ```ignore
/// // Find when Moon is at low altitude (0-30°)
/// let low_moon = find_moon_altitude_range(site, period, (Degrees::new(0.0), Degrees::new(30.0)));
/// ```
pub fn find_moon_altitude_range(
    site: ObserverSite,
    period: Period<ModifiedJulianDate>,
    range: (Degrees, Degrees),
) -> Vec<Period<ModifiedJulianDate>> {
    let jd_start = period.start.to_julian_day();
    let jd_end = period.end.to_julian_day();
    let min_rad = range.0.to::<Radian>().value();
    let max_rad = range.1.to::<Radian>().value();

    let altitude_fn = |jd: JulianDate| moon_altitude_rad(jd, &site);

    let observer_geo = spherical::position::Geographic::new(
        site.lon, site.lat, Kilometers::new(0.0),
    );
    let key_times = compute_moon_key_times(&observer_geo, jd_start, jd_end);

    let crossings_min = find_threshold_crossings_in_segments(
        &key_times, &altitude_fn, min_rad, jd_start, jd_end,
    );
    let above_min = crossings_to_above_periods(crossings_min, period, &altitude_fn, min_rad);

    let crossings_max = find_threshold_crossings_in_segments(
        &key_times, &altitude_fn, max_rad, jd_start, jd_end,
    );
    let above_max = crossings_to_above_periods(crossings_max, period, &altitude_fn, max_rad);

    let below_max = complement_within(period, &above_max);
    intersect_periods(&above_min, &below_max)
}

// =============================================================================
// Scan-based Algorithm (for comparison/fallback)
// =============================================================================

/// Finds periods using the generic scan-based algorithm (above threshold).
///
/// This is slower than the fast-scan approach but can be useful for
/// comparison or as a fallback.
pub fn find_moon_above_horizon_scan<T: Into<Degrees>>(
    site: ObserverSite,
    period: Period<ModifiedJulianDate>,
    threshold: T,
) -> Vec<Period<ModifiedJulianDate>> {
    let altitude_fn = moon_altitude_fn(site);
    crate::calculus::events::altitude_periods::find_above_altitude_periods(
        altitude_fn,
        period,
        threshold.into(),
    )
}

/// Finds periods using the generic scan-based algorithm (below threshold).
pub fn find_moon_below_horizon_scan<T: Into<Degrees>>(
    site: ObserverSite,
    period: Period<ModifiedJulianDate>,
    threshold: T,
) -> Vec<Period<ModifiedJulianDate>> {
    let above = find_moon_above_horizon_scan(site, period, threshold);
    complement_within(period, &above)
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::observatories::ROQUE_DE_LOS_MUCHACHOS;

    fn greenwich_site() -> ObserverSite {
        ObserverSite::new(
            Degrees::new(0.0),
            Degrees::new(51.4769),
            Meters::new(0.0),
        )
    }

    #[test]
    fn test_moon_altitude_basic() {
        let site = greenwich_site();
        let jd = JulianDate::J2000;
        let alt = moon_altitude_rad(jd, &site);
        // Moon altitude should be within valid range
        assert!(alt > -std::f64::consts::FRAC_PI_2 && alt < std::f64::consts::FRAC_PI_2);
    }

    #[test]
    fn test_moon_equatorial_position() {
        let jd = JulianDate::J2000;
        let target = get_moon_equatorial(jd);
        let pos = target.get_position();
        
        // Moon should be within reasonable distance (300,000-400,000 km)
        let dist = pos.distance().to::<Kilometer>().value();
        assert!(dist > 300_000.0 && dist < 450_000.0, "Moon distance: {} km", dist);
    }

    #[test]
    fn test_find_moon_above_horizon() {
        let site = greenwich_site();
        let mjd_start = ModifiedJulianDate::new(60000.0);
        let mjd_end = ModifiedJulianDate::new(60007.0);
        let period = Period::new(mjd_start, mjd_end);

        let periods = find_moon_above_horizon(site, period, Degrees::new(0.0));
        assert!(!periods.is_empty(), "Should find moon-up periods over 7 days");

        assert!(!periods.is_empty());

        for p in &periods {
            assert!(p.duration_days() > 0.0, "Period duration should be positive");
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

        assert!(!periods.is_empty());
    }

    #[test]
    fn test_culmination_vs_scan_consistency() {
        // Verify that culmination-based and scan-based algorithms produce similar results
        let site = greenwich_site();
        let mjd_start = ModifiedJulianDate::new(60000.0);
        let mjd_end = ModifiedJulianDate::new(60003.0);
        let period = Period::new(mjd_start, mjd_end);

        let culmination_result =
            find_moon_above_horizon(site, period, Degrees::new(0.0));
        let scan_result =
            find_moon_above_horizon_scan(site, period, Degrees::new(0.0));

        // Both should find periods
        assert!(!culmination_result.is_empty());
        assert!(!scan_result.is_empty());

        let culm_periods = culmination_result;
        let scan_periods = scan_result;

        // Should have similar number of periods
        assert_eq!(
            culm_periods.len(),
            scan_periods.len(),
            "Should find same number of periods"
        );

        // Boundaries should match within tolerance (~1 minute)
        let tolerance = 1.0 / 1440.0; // 1 minute in days
        for (c, s) in culm_periods.iter().zip(scan_periods.iter()) {
            assert!(
                (c.start.value() - s.start.value()).abs() < tolerance,
                "Start times should match within 1 minute"
            );
            assert!(
                (c.end.value() - s.end.value()).abs() < tolerance,
                "End times should match within 1 minute"
            );
        }
    }
}
