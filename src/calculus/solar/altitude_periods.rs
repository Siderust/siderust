// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Solar Altitude Window Periods
//!
//! Sun-specific convenience wrappers around the generic altitude-window routines
//! in [`crate::calculus::events::altitude_periods`].
//!
//! This module exists to keep `calculus::events` generic (body-agnostic) while
//! still providing ergonomic helpers for common solar concepts like twilight.

use crate::astro::JulianDate;
use crate::bodies::solar_system::Sun;
use crate::calculus::events::altitude_periods::{
    find_altitude_periods, AltitudeCondition,
};
use crate::calculus::events::{find_dynamic_extremas, Culmination};
use crate::coordinates::centers::ObserverSite;
use crate::coordinates::transform::Transform;
use crate::coordinates::{cartesian, spherical};
use crate::targets::Target;
use crate::time::{ModifiedJulianDate, Period};
use qtty::{AstronomicalUnit, Degrees, Kilometers, Radian};

/// Computes the Sun's altitude in **radians** at a given Julian Date and observer site.
/// Positive above the horizon, negative below.
pub fn sun_altitude_rad(jd: JulianDate, site: &ObserverSite) -> f64 {
    let horiz = Sun::get_horizontal::<AstronomicalUnit>(jd, *site);
    horiz.alt().to::<Radian>().value()
}

pub fn find_sun_altitude_periods_via_culminations(
    site: ObserverSite,
    period: Period<ModifiedJulianDate>,
    condition: AltitudeCondition,
) -> Option<Vec<Period<ModifiedJulianDate>>> {
    let jd_start = period.start.to_julian_day();
    let jd_end = period.end.to_julian_day();

    let altitude_fn = |jd: JulianDate| sun_altitude_rad(jd, &site);

    // Collect all boundary crossings (may be 1 or 2 boundaries depending on condition)
    let boundaries = match condition {
        AltitudeCondition::Below(threshold) | AltitudeCondition::Above(threshold) => {
            vec![threshold.to::<Radian>().value()]
        }
        AltitudeCondition::Between { min, max } => {
            vec![min.to::<Radian>().value(), max.to::<Radian>().value()]
        }
    };

    // Find upper/lower culminations (meridian crossings) for the Sun across the interval.
    //
    // We use `find_dynamic_extremas` to get a stable partition of the time axis into
    // quasi-monotonic segments, then bracket threshold crossings within each segment.
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

    let mut all_crossings: Vec<JulianDate> = Vec::new();

    for &boundary_rad in &boundaries {
        for window in key_times.windows(2) {
            let a = window[0];
            let b = window[1];
            if a.partial_cmp(&b) != Some(std::cmp::Ordering::Less) {
                continue;
            }

            let f_a = altitude_fn(a) - boundary_rad;
            let f_b = altitude_fn(b) - boundary_rad;

            const ROOT_EPS: f64 = 1e-12;
            if f_a.abs() < ROOT_EPS {
                all_crossings.push(a);
                continue;
            }
            if f_b.abs() < ROOT_EPS {
                all_crossings.push(b);
                continue;
            }

            if f_a * f_b < 0.0 {
                // Use Brent's method with pre-computed values: avoids 2 redundant VSOP evaluations
                if let Some(root) =
                    crate::calculus::root_finding::find_crossing_brent_with_values(a, b, f_a, f_b, &altitude_fn, boundary_rad)
                {
                    if root >= jd_start && root <= jd_end {
                        all_crossings.push(root);
                    }
                }
            }
        }
    }

    // Sort crossings chronologically
    all_crossings.sort_by(|a, b| a.partial_cmp(b).unwrap());

    // Deduplicate crossings that are very close (can happen at interval boundaries)
    const CROSS_DEDUPE_EPS: f64 = 1e-8; // ~1 ms
    all_crossings.dedup_by(|a, b| (a.value() - b.value()).abs() < CROSS_DEDUPE_EPS);

    // Classify each crossing: +1 = entering valid range, -1 = exiting valid range
    let mut labeled: Vec<(JulianDate, i32)> = Vec::new();
    for &root in &all_crossings {
        let dt = qtty::Days::new(10.0 * crate::calculus::root_finding::DERIVATIVE_STEP.value());
        let alt_before = altitude_fn(root - dt);
        let alt_after = altitude_fn(root + dt);

        let inside_before = condition.is_inside(alt_before);
        let inside_after = condition.is_inside(alt_after);

        if !inside_before && inside_after {
            labeled.push((root, 1)); // entering
        } else if inside_before && !inside_after {
            labeled.push((root, -1)); // exiting
        }
    }

    // Check if we start inside the valid range
    let start_altitude = altitude_fn(jd_start);
    let start_inside = condition.is_inside(start_altitude);

    // Build intervals by pairing enter/exit crossings
    let mut periods: Vec<Period<ModifiedJulianDate>> = Vec::new();

    if labeled.is_empty() {
        if start_inside {
            return Some(vec![period]);
        }
        return None;
    }

    let mut i = 0;

    // If we start inside and first crossing is an exit, add initial interval
    if start_inside && labeled[0].1 == -1 {
        let exit_mjd = ModifiedJulianDate::new(labeled[0].0.value() - 2400000.5);
        let mid = JulianDate::new((jd_start.value() + labeled[0].0.value()) * 0.5);
        if condition.is_inside(altitude_fn(mid)) {
            periods.push(Period::<ModifiedJulianDate>::new(period.start, exit_mjd));
        }
        i = 1;
    }

    // Process remaining crossings as enter/exit pairs
    while i < labeled.len() {
        if labeled[i].1 == 1 {
            let enter_jd = labeled[i].0;
            let enter_mjd = ModifiedJulianDate::new(enter_jd.value() - 2400000.5);

            let exit_mjd = if i + 1 < labeled.len() && labeled[i + 1].1 == -1 {
                let exit_jd = labeled[i + 1].0;
                i += 2;
                ModifiedJulianDate::new(exit_jd.value() - 2400000.5)
            } else {
                i += 1;
                period.end
            };

            let mid = JulianDate::new((enter_jd.value() + exit_mjd.to_julian_day().value()) * 0.5);
            if condition.is_inside(altitude_fn(mid)) {
                periods.push(Period::<ModifiedJulianDate>::new(enter_mjd, exit_mjd));
            }
        } else {
            i += 1;
        }
    }

    if periods.is_empty() {
        None
    } else {
        Some(periods)
    }
}

/// Common twilight types.
///
/// Re-exported from [`crate::calculus::solar::night_types`].
pub use super::night_types::Twilight;

/// Finds night periods (Sun below `twilight`) inside `period`.
pub fn find_night_periods<T: Into<Degrees>>(
    site: ObserverSite,
    period: Period<ModifiedJulianDate>,
    twilight: T,
) -> Option<Vec<Period<ModifiedJulianDate>>> {
    let tw: Degrees = twilight.into();
    find_sun_altitude_periods_via_culminations(site, period, AltitudeCondition::below(tw))
}

/// Finds night periods (Sun below `twilight`) using the generic scan+refine algorithm.
///
/// Prefer [`find_night_periods`] unless you specifically want to compare behavior.
pub fn find_night_periods_scan<T: Into<Degrees>>(
    site: ObserverSite,
    period: Period<ModifiedJulianDate>,
    twilight: T,
) -> Option<Vec<Period<ModifiedJulianDate>>> {
    let altitude_fn = |jd: JulianDate| sun_altitude_rad(jd, &site);
    let tw: Degrees = twilight.into();
    find_altitude_periods(altitude_fn, period, AltitudeCondition::below(tw))
}

/// Finds day periods (Sun above `twilight`) inside `period`.
pub fn find_day_periods<T: Into<Degrees>>(
    site: ObserverSite,
    period: Period<ModifiedJulianDate>,
    twilight: T,
) -> Option<Vec<Period<ModifiedJulianDate>>> {
    let tw: Degrees = twilight.into();
    find_sun_altitude_periods_via_culminations(site, period, AltitudeCondition::above(tw))
}

/// Finds day periods (Sun above `twilight`) using the generic scan+refine algorithm.
///
/// Prefer [`find_day_periods`] unless you specifically want to compare behavior.
pub fn find_day_periods_scan<T: Into<Degrees>>(
    site: ObserverSite,
    period: Period<ModifiedJulianDate>,
    twilight: T,
) -> Option<Vec<Period<ModifiedJulianDate>>> {
    let altitude_fn = |jd: JulianDate| sun_altitude_rad(jd, &site);
    let tw: Degrees = twilight.into();
    find_altitude_periods(altitude_fn, period, AltitudeCondition::above(tw))
}

/// Finds periods where Sun altitude is within `range` (min, max) inside `period`.
pub fn find_sun_range_periods(
    site: ObserverSite,
    period: Period<ModifiedJulianDate>,
    range: (Degrees, Degrees),
) -> Option<Vec<Period<ModifiedJulianDate>>> {
    find_sun_altitude_periods_via_culminations(
        site,
        period,
        AltitudeCondition::between(range.0, range.1),
    )
}

/// Finds periods where Sun altitude is within `range` (min, max) using the generic scan+refine algorithm.
///
/// Prefer [`find_sun_range_periods`] unless you specifically want to compare behavior.
pub fn find_sun_range_periods_scan(
    site: ObserverSite,
    period: Period<ModifiedJulianDate>,
    range: (Degrees, Degrees),
) -> Option<Vec<Period<ModifiedJulianDate>>> {
    let altitude_fn = |jd: JulianDate| sun_altitude_rad(jd, &site);
    find_altitude_periods(
        altitude_fn,
        period,
        AltitudeCondition::between(range.0, range.1),
    )
}

/// Standard twilight threshold definitions (Sun center altitude).
pub mod twilight {
    use qtty::Degrees;

    /// Civil twilight: Sun center 6° below horizon (-6°)
    pub const CIVIL: Degrees = Degrees::new(-6.0);

    /// Nautical twilight: Sun center 12° below horizon (-12°)
    pub const NAUTICAL: Degrees = Degrees::new(-12.0);

    /// Astronomical twilight: Sun center 18° below horizon (-18°)
    pub const ASTRONOMICAL: Degrees = Degrees::new(-18.0);

    /// Sunrise/sunset: Sun center at geometric horizon (0°)
    /// Note: For apparent sunrise/sunset, use -0.833° to account for refraction
    pub const HORIZON: Degrees = Degrees::new(0.0);

    /// Apparent sunrise/sunset accounting for atmospheric refraction (-0.833°)
    pub const APPARENT_HORIZON: Degrees = Degrees::new(-0.833);
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
        assert!(alt > -std::f64::consts::FRAC_PI_2 && alt < std::f64::consts::FRAC_PI_2);
    }

    #[test]
    fn test_find_night_periods() {
        let site = greenwich_site();
        let mjd_start = ModifiedJulianDate::new(60000.0);
        let mjd_end = ModifiedJulianDate::new(60007.0);
        let period = Period::new(mjd_start, mjd_end);

        let nights = find_night_periods(site, period, twilight::ASTRONOMICAL);
        assert!(
            nights.is_some(),
            "Should find night periods at 51° latitude"
        );

        let nights = nights.unwrap();
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

        assert!(nights.is_some(), "Should find night periods using range");
        let nights = nights.unwrap();
        assert!(!nights.is_empty(), "Should have at least one night period");

        let nautical =
            find_sun_range_periods(site, period, (Degrees::new(-18.0), Degrees::new(-12.0)));

        assert!(nautical.is_some(), "Should find nautical twilight periods");
    }
}
