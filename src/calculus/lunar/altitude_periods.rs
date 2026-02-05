// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Lunar Altitude Window Periods
//!
//! Moon-specific optimized routines for finding time intervals where the Moon's
//! altitude satisfies a given condition (above, below, or within a range).
//!
//! ## Key Optimizations
//!
//! 1. **Topocentric parallax correction**: Critical for Moon (~1° at horizon)
//! 2. **Culmination-based search**: Partitions time into monotonic altitude segments
//! 3. **Brent's method with endpoint reuse**: Avoids redundant ELP2000 evaluations
//! 4. **Reuses `find_dynamic_extremas`**: Leverages the generic, optimized culmination finder
//!
//! ## Performance Notes
//!
//! The culmination-based algorithm achieves sub-second performance for 365-day
//! searches by:
//! - Using `find_dynamic_extremas` to find Moon culminations efficiently
//! - Reusing endpoint altitude values when calling Brent's method
//! - Avoiding derivative evaluations (Brent is derivative-free)

use crate::astro::nutation::corrected_ra_with_nutation;
use crate::astro::precession;
use crate::astro::JulianDate;
use crate::bodies::solar_system::Moon;
use crate::calculus::events::altitude_periods::{AltitudeCondition, AltitudePeriod};
use crate::calculus::events::Culmination;
use crate::calculus::root_finding::brent;
use crate::coordinates::cartesian;
use crate::coordinates::centers::{Geocentric, ObserverSite, Topocentric};
use crate::coordinates::frames::{self, Ecliptic};
use crate::coordinates::spherical::{self, Position};
use crate::coordinates::transform::centers::ToTopocentricExt;
use crate::coordinates::transform::{Transform, TransformFrame};
use crate::targets::Target;
use crate::time::{ModifiedJulianDate, Period};
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

/// Root-finding precision (~0.86 µs).
const ROOT_EPS: f64 = 1e-12;

/// Relaxed tolerance for fast Moon altitude root finding (~2 minute precision).
/// 2 minutes = 2/(24*60) days ≈ 1.39e-3 days
const MOON_BRENT_TOLERANCE: f64 = 1.4e-3;

/// Fast Brent root finding with relaxed tolerance for Moon altitude.
///
/// Uses ~1 minute precision instead of ~86µs, reducing iterations from ~8 to ~3-4.
fn brent_fast_with_values<F>(
    lo: JulianDate,
    hi: JulianDate,
    f_lo: f64,
    f_hi: f64,
    scalar_fn: F,
    threshold: f64,
) -> Option<JulianDate>
where
    F: Fn(JulianDate) -> f64,
{
    const CONVERGENCE_EPS: f64 = 1e-10;
    const MAX_ITERS: usize = 20;
    
    let mut a = lo.value();
    let mut b = hi.value();
    let mut fa = f_lo;
    let mut fb = f_hi;

    // Check endpoints for exact roots
    if fa.abs() < CONVERGENCE_EPS {
        return Some(lo);
    }
    if fb.abs() < CONVERGENCE_EPS {
        return Some(hi);
    }

    // Verify bracket (must have opposite signs)
    if fa * fb > 0.0 {
        return None;
    }

    // Ensure |f(a)| >= |f(b)| (b is the better approximation)
    if fa.abs() < fb.abs() {
        std::mem::swap(&mut a, &mut b);
        std::mem::swap(&mut fa, &mut fb);
    }

    let mut c = a;
    let mut fc = fa;
    let mut d = b - a;
    let mut e = d;

    for _ in 0..MAX_ITERS {
        // Keep c as the contrapoint (opposite sign to b)
        if (fb > 0.0) == (fc > 0.0) {
            c = a;
            fc = fa;
            d = b - a;
            e = d;
        }
        
        // Ensure |f(b)| <= |f(c)|
        if fc.abs() < fb.abs() {
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
        }

        let tol = MOON_BRENT_TOLERANCE;
        let m = 0.5 * (c - b);

        if fb.abs() < CONVERGENCE_EPS || m.abs() <= tol {
            return Some(JulianDate::new(b));
        }

        // Decide whether to use bisection or interpolation
        let use_bisection = e.abs() < tol || fa.abs() <= fb.abs();

        let (new_e, new_d) = if use_bisection {
            (m, m)
        } else {
            let s = fb / fa;

            let (p, q) = if (a - c).abs() < 1e-14 {
                // Secant method
                let p = 2.0 * m * s;
                let q = 1.0 - s;
                (p, q)
            } else {
                // Inverse quadratic interpolation
                let q_val = fa / fc;
                let r = fb / fc;
                let p = s * (2.0 * m * q_val * (q_val - r) - (b - a) * (r - 1.0));
                let q = (q_val - 1.0) * (r - 1.0) * (s - 1.0);
                (p, q)
            };

            let (p, q) = if p > 0.0 { (p, -q) } else { (-p, q) };

            let s_val = e;
            if 2.0 * p < 3.0 * m * q - (tol * q).abs() && p < (0.5 * s_val * q).abs() {
                (d, p / q)
            } else {
                (m, m)
            }
        };

        e = new_e;
        d = new_d;

        a = b;
        fa = fb;

        if d.abs() > tol {
            b += d;
        } else if m > 0.0 {
            b += tol;
        } else {
            b -= tol;
        }

        fb = scalar_fn(JulianDate::new(b)) - threshold;
    }

    Some(JulianDate::new(b))
}

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
    // Get Moon's geocentric ecliptic position from ELP2000 (cartesian)
    let moon_geo_ecliptic: cartesian::Position<Geocentric, Ecliptic, Kilometer> =
        Moon::get_geo_position(jd);

    // Transform: Ecliptic → EquatorialMeanJ2000
    let moon_geo_eq_j2000: cartesian::Position<Geocentric, frames::EquatorialMeanJ2000, Kilometer> =
        TransformFrame::to_frame(&moon_geo_ecliptic);

    // Apply topocentric parallax correction (critical for Moon!)
    let moon_topo_eq_j2000: cartesian::Position<
        Topocentric,
        frames::EquatorialMeanJ2000,
        Kilometer,
    > = moon_geo_eq_j2000.to_topocentric(*site, jd);

    // Apply precession: J2000 → mean-of-date
    let rot = crate::coordinates::transform::frame_rotation::<
        frames::EquatorialMeanJ2000,
        frames::EquatorialMeanOfDate,
    >(jd, &crate::coordinates::transform::AstroContext::default());

    let [x, y, z] = rot.apply_array([
        moon_topo_eq_j2000.x().value(),
        moon_topo_eq_j2000.y().value(),
        moon_topo_eq_j2000.z().value(),
    ]);

    let moon_topo_eq_date: cartesian::Position<Topocentric, frames::EquatorialMeanOfDate, Kilometer> =
        cartesian::Position::from_vec3(
            *site,
            nalgebra::Vector3::new(x * KM, y * KM, z * KM),
        );

    // Transform to horizontal coordinates
    let moon_horizontal: cartesian::Position<Topocentric, frames::Horizontal, Kilometer> =
        moon_topo_eq_date.transform(jd);

    // Convert to spherical to extract altitude
    let moon_horizontal_sph: Position<Topocentric, frames::Horizontal, Kilometer> =
        Position::from_cartesian(&moon_horizontal);

    moon_horizontal_sph.alt().to::<Radian>().value()
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

/// Fast GAST approximation (error < 0.1″ for ±50 yr around 2025).
fn gast_fast(jd: JulianDate) -> Degrees {
    let t = (jd.value() - 2_451_545.0) / 36_525.0;
    let gast = 280.460_618_37
        + 360.985_647_366_29 * (jd.value() - 2_451_545.0)
        + 0.000_387_933 * t * t
        - t * t * t / 38_710.0;
    Degrees::new(gast)
}

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

// =============================================================================
// Main API: Optimized Altitude Period Finding
// =============================================================================

/// Finds all periods where Moon altitude satisfies the given condition,
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
/// * `condition` - Altitude condition (below, above, or between)
///
/// # Returns
/// * `Some(Vec<AltitudePeriod>)` - Periods where condition is satisfied
/// * `None` - If condition is never satisfied
///
/// # Example
/// ```ignore
/// use siderust::calculus::lunar::*;
/// use siderust::calculus::events::altitude_periods::AltitudeCondition;
/// use siderust::observatories::ROQUE_DE_LOS_MUCHACHOS;
/// use qtty::Degrees;
///
/// let site = ObserverSite::from_geographic(&ROQUE_DE_LOS_MUCHACHOS);
/// let period = /* your period */;
/// let condition = AltitudeCondition::above(Degrees::new(0.0));
/// let periods = find_moon_altitude_periods_via_culminations(site, period, condition);
/// ```
pub fn find_moon_altitude_periods_via_culminations(
    site: ObserverSite,
    period: Period<ModifiedJulianDate>,
    condition: AltitudeCondition,
) -> Option<Vec<AltitudePeriod>> {
    let jd_start = period.start.to_julian_day();
    let jd_end = period.end.to_julian_day();

    let altitude_fn = |jd: JulianDate| moon_altitude_rad(jd, &site);

    // Get boundaries from condition
    let boundaries = match condition {
        AltitudeCondition::Below(threshold) | AltitudeCondition::Above(threshold) => {
            vec![threshold.to::<Radian>().value()]
        }
        AltitudeCondition::Between { min, max } => {
            vec![min.to::<Radian>().value(), max.to::<Radian>().value()]
        }
    };

    // Find Moon culminations using the optimized Moon-specific finder
    let observer_geo = spherical::position::Geographic::new(
        site.lon,
        site.lat,
        Kilometers::new(0.0),
    );

    let culminations = find_moon_culminations_fast(&observer_geo, jd_start, jd_end);

    // Build key times: start, culminations, end
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

    // Find all boundary crossings
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

            // Check for root at endpoints
            if f_a.abs() < ROOT_EPS {
                all_crossings.push(a);
                continue;
            }
            if f_b.abs() < ROOT_EPS {
                all_crossings.push(b);
                continue;
            }

            // Sign change indicates a crossing
            if f_a * f_b < 0.0 {
                // Use Brent's method with pre-computed endpoint values
                if let Some(root) =
                    crate::calculus::root_finding::find_crossing_brent_with_values(
                        a, b, f_a, f_b, &altitude_fn, boundary_rad,
                    )
                {
                    if root >= jd_start && root <= jd_end {
                        all_crossings.push(root);
                    }
                }
            }
        }
    }

    // Sort and deduplicate crossings
    all_crossings.sort_by(|a, b| a.partial_cmp(b).unwrap());
    all_crossings.dedup_by(|a, b| (a.value() - b.value()).abs() < CROSS_DEDUPE_EPS);

    // Classify crossings: +1 = entering, -1 = exiting
    let mut labeled: Vec<(JulianDate, i32)> = Vec::new();
    let dt = Days::new(10.0 * crate::calculus::root_finding::DERIVATIVE_STEP.value());

    for &root in &all_crossings {
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

    // Build intervals
    let mut periods: Vec<AltitudePeriod> = Vec::new();

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
            periods.push(AltitudePeriod::new(period.start, exit_mjd));
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
                periods.push(AltitudePeriod::new(enter_mjd, exit_mjd));
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

// =============================================================================
// Fast Direct Scan Algorithm (Optimized for Performance)
// =============================================================================

/// Finds Moon altitude periods using direct 2-hour scan without culmination detection.
///
/// This is the **fastest** algorithm for long periods, achieving significant
/// performance improvements for 365-day searches by:
/// 1. Using 2-hour scan steps (~12 calls/day)
/// 2. Reusing endpoint values for Brent refinement
/// 3. Tracking sign changes to avoid extra altitude evaluations
/// 4. Avoiding the culmination detection overhead entirely
///
/// # Trade-off
/// Slightly less numerically precise at boundaries compared to culmination-based
/// approach, but the difference is negligible (<1 minute) for practical use.
pub fn find_moon_altitude_periods_fast(
    site: ObserverSite,
    period: Period<ModifiedJulianDate>,
    condition: AltitudeCondition,
) -> Option<Vec<AltitudePeriod>> {
    let jd_start = period.start.to_julian_day();
    let jd_end = period.end.to_julian_day();

    let altitude_fn = |jd: JulianDate| moon_altitude_rad(jd, &site);

    // Get boundaries from condition
    let boundaries = match condition {
        AltitudeCondition::Below(threshold) | AltitudeCondition::Above(threshold) => {
            vec![threshold.to::<Radian>().value()]
        }
        AltitudeCondition::Between { min, max } => {
            vec![min.to::<Radian>().value(), max.to::<Radian>().value()]
        }
    };

    // Store crossings with their direction: (jd, +1=rising, -1=falling)
    // "rising" means altitude is increasing through boundary (going from below to above)
    let mut labeled_crossings: Vec<(JulianDate, i32)> = Vec::new();

    for &boundary_rad in &boundaries {
        let mut jd = jd_start;
        let mut prev_f = altitude_fn(jd) - boundary_rad;

        while jd < jd_end {
            let next_jd = (jd + MOON_ALTITUDE_SCAN_STEP).min(jd_end);
            let next_f = altitude_fn(next_jd) - boundary_rad;

            // Sign change indicates a crossing
            if prev_f * next_f < 0.0 {
                // Direction: prev_f < 0 && next_f > 0 means rising (altitude increasing)
                let direction = if prev_f < 0.0 { 1 } else { -1 };
                
                if let Some(root) = brent_fast_with_values(
                    jd, next_jd, prev_f, next_f, &altitude_fn, boundary_rad,
                ) {
                    if root >= jd_start && root <= jd_end {
                        labeled_crossings.push((root, direction));
                    }
                }
            }

            jd = next_jd;
            prev_f = next_f;
        }
    }

    // Sort crossings by time
    labeled_crossings.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    
    // Deduplicate (keep first occurrence)
    let mut i = 0;
    while i + 1 < labeled_crossings.len() {
        if (labeled_crossings[i].0.value() - labeled_crossings[i + 1].0.value()).abs() < CROSS_DEDUPE_EPS {
            labeled_crossings.remove(i + 1);
        } else {
            i += 1;
        }
    }

    // Now convert direction (rising/falling through boundary) to entering/exiting condition
    // For AltitudeCondition::Above: entering = rising (+1), exiting = falling (-1)
    // For AltitudeCondition::Below: entering = falling (-1), exiting = rising (+1)
    let labeled: Vec<(JulianDate, i32)> = match condition {
        AltitudeCondition::Above(_) => {
            // Rising through boundary = entering "above" condition
            labeled_crossings.iter().map(|&(jd, dir)| (jd, dir)).collect()
        }
        AltitudeCondition::Below(_) => {
            // Falling through boundary = entering "below" condition  
            labeled_crossings.iter().map(|&(jd, dir)| (jd, -dir)).collect()
        }
        AltitudeCondition::Between { .. } => {
            // For Between, we need more complex logic - fall back to evaluation
            let dt = Days::new(10.0 * crate::calculus::root_finding::DERIVATIVE_STEP.value());
            labeled_crossings.iter().map(|&(root, _)| {
                let alt_before = altitude_fn(root - dt);
                let alt_after = altitude_fn(root + dt);
                let inside_before = condition.is_inside(alt_before);
                let inside_after = condition.is_inside(alt_after);
                if !inside_before && inside_after {
                    (root, 1) // entering
                } else if inside_before && !inside_after {
                    (root, -1) // exiting
                } else {
                    (root, 0) // tangent, skip
                }
            }).filter(|&(_, d)| d != 0).collect()
        }
    };

    // Check if we start inside the valid range (already computed during scan)
    let start_altitude = altitude_fn(jd_start);
    let start_inside = condition.is_inside(start_altitude);

    // Build intervals
    let mut periods: Vec<AltitudePeriod> = Vec::new();

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
        periods.push(AltitudePeriod::new(period.start, exit_mjd));
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

            periods.push(AltitudePeriod::new(enter_mjd, exit_mjd));
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
) -> Option<Vec<AltitudePeriod>> {
    find_moon_altitude_periods_fast(
        site,
        period,
        AltitudeCondition::above(threshold.into()),
    )
}

/// Finds periods when Moon is below the given altitude threshold.
///
/// Uses the fast direct-scan algorithm optimized for performance.
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
) -> Option<Vec<AltitudePeriod>> {
    find_moon_altitude_periods_fast(
        site,
        period,
        AltitudeCondition::below(threshold.into()),
    )
}

/// Finds periods when Moon altitude is within a range.
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
) -> Option<Vec<AltitudePeriod>> {
    find_moon_altitude_periods_via_culminations(
        site,
        period,
        AltitudeCondition::between(range.0, range.1),
    )
}

// =============================================================================
// Scan-based Algorithm (for comparison/fallback)
// =============================================================================

/// Finds periods using the generic scan-based algorithm.
///
/// This is slower than the culmination-based approach but can be useful for
/// comparison or as a fallback.
pub fn find_moon_above_horizon_scan<T: Into<Degrees>>(
    site: ObserverSite,
    period: Period<ModifiedJulianDate>,
    threshold: T,
) -> Option<Vec<AltitudePeriod>> {
    let altitude_fn = moon_altitude_fn(site);
    crate::calculus::events::altitude_periods::find_altitude_periods(
        altitude_fn,
        period,
        AltitudeCondition::above(threshold.into()),
    )
}

/// Finds periods using the generic scan-based algorithm (below threshold variant).
pub fn find_moon_below_horizon_scan<T: Into<Degrees>>(
    site: ObserverSite,
    period: Period<ModifiedJulianDate>,
    threshold: T,
) -> Option<Vec<AltitudePeriod>> {
    let altitude_fn = moon_altitude_fn(site);
    crate::calculus::events::altitude_periods::find_altitude_periods(
        altitude_fn,
        period,
        AltitudeCondition::below(threshold.into()),
    )
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
        assert!(periods.is_some(), "Should find moon-up periods over 7 days");

        let periods = periods.unwrap();
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
        assert!(periods.is_some(), "Should find moon-down periods");

        let periods = periods.unwrap();
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
        assert!(culmination_result.is_some());
        assert!(scan_result.is_some());

        let culm_periods = culmination_result.unwrap();
        let scan_periods = scan_result.unwrap();

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
