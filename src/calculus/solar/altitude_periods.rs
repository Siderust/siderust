// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Solar Altitude Window Periods
//!
//! Sun-specific convenience wrappers around the generic altitude-window routines
//! in [`crate::calculus::events::altitude_periods`].
//!
//! The core algorithm finds **above-threshold** periods. Below and between
//! variants are derived via [`crate::time::complement_within`] and
//! [`crate::time::intersect_periods`] at negligible cost.
//!
//! ## Performance
//!
//! The culmination-based variants (`find_day_periods`, `find_night_periods`,
//! `find_sun_range_periods`) partition the time axis using solar culminations,
//! then bracket threshold crossings within each monotonic segment.
//! This avoids the coarse 10-minute scan and reduces the number of VSOP87
//! evaluations significantly.
//!
//! The `_scan` variants fall back to the generic scan+refine algorithm from
//! [`crate::calculus::events::altitude_periods`].

use crate::astro::JulianDate;
use crate::bodies::solar_system::Sun;
use crate::calculus::events::altitude_periods::{
    crossings_to_above_periods, find_above_altitude_periods,
    find_threshold_crossings_in_segments,
};
use crate::calculus::events::{find_dynamic_extremas, Culmination};
use crate::coordinates::centers::ObserverSite;
use crate::coordinates::transform::Transform;
use crate::coordinates::{cartesian, spherical};
use crate::targets::Target;
use crate::time::{complement_within, intersect_periods, ModifiedJulianDate, Period};
use qtty::{AstronomicalUnit, Degrees, Kilometers, Quantity, Radian};

/// Computes the Sun's altitude in **radians** at a given Julian Date and observer site.
/// Positive above the horizon, negative below.
pub fn sun_altitude_rad(jd: JulianDate, site: &ObserverSite) -> Quantity<Radian> {
    let horiz = Sun::get_horizontal::<AstronomicalUnit>(jd, *site);
    horiz.alt().to::<Radian>()
}

// =============================================================================
// Internal: Culmination-based key-times computation
// =============================================================================

/// Computes the sorted, deduplicated list of key times (start, culminations, end)
/// that partition the time axis into quasi-monotonic altitude segments for the Sun.
fn compute_solar_key_times(
    site: &ObserverSite,
    jd_start: JulianDate,
    jd_end: JulianDate,
) -> Vec<JulianDate> {
    let observer_geo =
        spherical::position::Geographic::new(site.lon, site.lat, Kilometers::new(0.0));

    let get_equatorial = |jd: JulianDate| {
        let helio = cartesian::position::Ecliptic::<AstronomicalUnit>::CENTER;
        let geo_cart: cartesian::position::EquatorialMeanJ2000<AstronomicalUnit> =
            helio.transform(jd);
        let geo_sph: spherical::position::EquatorialMeanJ2000<AstronomicalUnit> =
            spherical::Position::from_cartesian(&geo_cart);
        Target::new_static(geo_sph, jd)
    };

    let culminations = find_dynamic_extremas(get_equatorial, &observer_geo, jd_start, jd_end);

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

    const DEDUPE_EPS: f64 = 1e-10; // ~0.86 ms
    key_times.dedup_by(|a, b| (a.value() - b.value()).abs() < DEDUPE_EPS);

    key_times
}

// =============================================================================
// Core: Culmination-based above-threshold finder
// =============================================================================

/// Finds periods where Sun altitude is **above** `threshold` using
/// culmination-based partitioning.
///
/// This is the primary building block. All other solar altitude finders
/// (day, night, range, and their scan variants) are thin wrappers around
/// this function and the set operations on period vectors.
pub fn find_sun_altitude_periods_via_culminations<T: Into<Degrees>>(
    site: ObserverSite,
    period: Period<ModifiedJulianDate>,
    threshold: T,
) -> Vec<Period<ModifiedJulianDate>> {
    let jd_start = period.start.to_julian_day();
    let jd_end = period.end.to_julian_day();
    let threshold_rad = threshold.into().to::<Radian>().value();

    let altitude_fn = |jd: JulianDate| sun_altitude_rad(jd, &site).value();

    let key_times = compute_solar_key_times(&site, jd_start, jd_end);

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
// Public convenience wrappers
// =============================================================================

/// Finds day periods (Sun **above** `twilight`) inside `period`.
///
/// Uses culmination-based partitioning for optimal performance.
pub fn find_day_periods<T: Into<Degrees>>(
    site: ObserverSite,
    period: Period<ModifiedJulianDate>,
    twilight: T,
) -> Vec<Period<ModifiedJulianDate>> {
    find_sun_altitude_periods_via_culminations(site, period, twilight)
}

/// Finds night periods (Sun **below** `twilight`) inside `period`.
///
/// Equivalent to the complement of [`find_day_periods`] within `period`.
pub fn find_night_periods<T: Into<Degrees>>(
    site: ObserverSite,
    period: Period<ModifiedJulianDate>,
    twilight: T,
) -> Vec<Period<ModifiedJulianDate>> {
    let days = find_day_periods(site, period, twilight);
    complement_within(period, &days)
}

/// Finds periods where Sun altitude is within `range` `(min, max)`.
///
/// Computed as `above(min) ∩ complement(above(max))`. The culmination
/// key-times are computed only once and reused for both thresholds.
pub fn find_sun_range_periods(
    site: ObserverSite,
    period: Period<ModifiedJulianDate>,
    range: (Degrees, Degrees),
) -> Vec<Period<ModifiedJulianDate>> {
    let jd_start = period.start.to_julian_day();
    let jd_end = period.end.to_julian_day();
    let min_rad = range.0.to::<Radian>().value();
    let max_rad = range.1.to::<Radian>().value();

    let altitude_fn = |jd: JulianDate| sun_altitude_rad(jd, &site).value();

    // Compute key_times ONCE
    let key_times = compute_solar_key_times(&site, jd_start, jd_end);

    // Find above-min periods
    let crossings_min = find_threshold_crossings_in_segments(
        &key_times, &altitude_fn, min_rad, jd_start, jd_end,
    );
    let above_min = crossings_to_above_periods(crossings_min, period, &altitude_fn, min_rad);

    // Find above-max periods
    let crossings_max = find_threshold_crossings_in_segments(
        &key_times, &altitude_fn, max_rad, jd_start, jd_end,
    );
    let above_max = crossings_to_above_periods(crossings_max, period, &altitude_fn, max_rad);

    // Between(min, max) = above(min) ∩ complement(above(max))
    let below_max = complement_within(period, &above_max);
    intersect_periods(&above_min, &below_max)
}

// =============================================================================
// Scan-based variants (for comparison / fallback)
// =============================================================================

/// Finds day periods using the generic scan+refine algorithm.
///
/// Prefer [`find_day_periods`] unless you specifically want to compare behavior.
pub fn find_day_periods_scan<T: Into<Degrees>>(
    site: ObserverSite,
    period: Period<ModifiedJulianDate>,
    twilight: T,
) -> Vec<Period<ModifiedJulianDate>> {
    let altitude_fn = |jd: JulianDate| sun_altitude_rad(jd, &site).value();
    find_above_altitude_periods(altitude_fn, period, twilight.into())
}

/// Finds night periods using the generic scan+refine algorithm.
///
/// Prefer [`find_night_periods`] unless you specifically want to compare behavior.
pub fn find_night_periods_scan<T: Into<Degrees>>(
    site: ObserverSite,
    period: Period<ModifiedJulianDate>,
    twilight: T,
) -> Vec<Period<ModifiedJulianDate>> {
    let days = find_day_periods_scan(site, period, twilight);
    complement_within(period, &days)
}

/// Finds periods where Sun altitude is within `range` using the generic scan+refine algorithm.
///
/// Prefer [`find_sun_range_periods`] unless you specifically want to compare behavior.
pub fn find_sun_range_periods_scan(
    site: ObserverSite,
    period: Period<ModifiedJulianDate>,
    range: (Degrees, Degrees),
) -> Vec<Period<ModifiedJulianDate>> {
    let altitude_fn = |jd: JulianDate| sun_altitude_rad(jd, &site).value();
    let above_min = find_above_altitude_periods(&altitude_fn, period, range.0);
    let above_max = find_above_altitude_periods(&altitude_fn, period, range.1);
    let below_max = complement_within(period, &above_max);
    intersect_periods(&above_min, &below_max)
}

#[cfg(test)]
mod tests {
    use super::*;
    use qtty::*;

    fn greenwich_site() -> ObserverSite {
        ObserverSite::new(
            Degrees::new(0.0),
            Degrees::new(51.4769),
            Quantity::<Meter>::new(0.0),
        )
    }

    #[test]
    fn test_sun_altitude_basic() {
        let site = greenwich_site();
        let jd = JulianDate::J2000;
        let alt = sun_altitude_rad(jd, &site);
        assert!(alt.value() > -std::f64::consts::FRAC_PI_2 && alt.value() < std::f64::consts::FRAC_PI_2);
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

        assert!(!nights.is_empty(), "Should have at least one night period");

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
        assert!(!nights.is_empty(), "Should have at least one night period");

        let nautical =
            find_sun_range_periods(site, period, (Degrees::new(-18.0), Degrees::new(-12.0)));

        assert!(!nautical.is_empty(), "Should find nautical twilight periods");
    }
}
