// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Altitude Window Periods
//!
//! Generic, body-agnostic tools for finding time intervals where a celestial
//! body's altitude is above a given threshold.
//!
//! ## Core Idea
//!
//! Every algorithm in this module (and in the solar / lunar modules) reduces to
//! a single primitive: **find periods where `altitude > threshold`**.
//!
//! * *Below*-threshold periods (e.g. night) are the **complement** of
//!   above-threshold periods — see [`crate::time::complement_within`].
//! * *Between*-range periods (e.g. twilight bands) are the **intersection** of
//!   `above(min)` with the complement of `above(max)` — see
//!   [`crate::time::intersect_periods`].
//!
//! The set operations run in O(n) on the resulting period vectors, adding zero
//! ephemeris cost.
//!
//! ## Algorithm
//!
//! The public [`find_above_altitude_periods`] function uses a **scan + refine**
//! approach:
//!
//! 1. Coarse scan at 10-minute steps to detect sign changes of
//!    `altitude(t) − threshold`.
//! 2. Brent's method refinement for each bracket (derivative-free, ~1–2
//!    evaluations / iteration).
//! 3. Classification of each crossing as *entering* or *exiting* the
//!    above-threshold region.
//! 4. Pairing of crossings into contiguous intervals with midpoint validation.
//!
//! Body-specific modules (solar, lunar) provide faster variants that replace
//! step 1 with culmination-based partitioning, sharing steps 2–4 through the
//! [`crossings_to_above_periods`] and [`find_threshold_crossings_in_segments`]
//! helpers exported here at `pub(crate)` visibility.
//!
//! ## Accuracy
//!
//! Timing precision is ~1 µs (Brent convergence criterion `1e-11` days ≈ 0.86 µs).

use crate::astro::JulianDate;
use crate::calculus::math_core;
use crate::time::{ModifiedJulianDate, Period};
use qtty::{Day, Days, Degrees, Minutes, Radian};

// =============================================================================
// Constants
// =============================================================================

/// Scan step for coarse bracket detection (10 minutes).
const SCAN_STEP: Minutes = Minutes::new(10.0);

// =============================================================================
// Shared Helpers (pub(crate))
// =============================================================================

/// Finds all crossings of a threshold within pre-computed time segments.
///
/// Each consecutive pair in `key_times` defines a quasi-monotonic segment
/// (e.g. between successive culminations), so at most one crossing per segment
/// is expected. Crossings are found via Brent's method with pre-computed
/// endpoint values, avoiding two redundant altitude evaluations per segment.
///
/// Returns the raw (unsorted) list of crossing JDs.
pub(crate) fn find_threshold_crossings_in_segments<F>(
    key_times: &[JulianDate],
    altitude_fn: &F,
    threshold_rad: f64,
    jd_start: JulianDate,
    jd_end: JulianDate,
) -> Vec<JulianDate>
where
    F: Fn(JulianDate) -> f64,
{
    // Convert JulianDate key_times to f64 for math_core
    let key_f64: Vec<f64> = key_times.iter().map(|jd| jd.value()).collect();

    // Wrap the altitude_fn to operate on plain f64
    let f = |t: f64| altitude_fn(JulianDate::new(t));

    let raw = math_core::intervals::find_crossings_in_segments(
        &key_f64,
        &f,
        threshold_rad,
        jd_start.value(),
        jd_end.value(),
    );

    raw.into_iter().map(JulianDate::new).collect()
}

/// Converts a list of threshold crossings into "above threshold" periods.
///
/// This is the shared pipeline used by all altitude-finding algorithms:
///
/// 1. Sort and deduplicate crossings chronologically.
/// 2. Classify each crossing as *entering above* (+1) or *exiting above* (−1)
///    by probing the altitude just before and after the root.
/// 3. Pair crossings into contiguous `Period`s, validating each with a
///    midpoint check.
pub(crate) fn crossings_to_above_periods<F>(
    mut crossings: Vec<JulianDate>,
    period: Period<ModifiedJulianDate>,
    altitude_fn: &F,
    threshold_rad: f64,
) -> Vec<Period<ModifiedJulianDate>>
where
    F: Fn(JulianDate) -> f64,
{
    let jd_start = period.start.to_julian_day();
    let jd_end = period.end.to_julian_day();

    // Sort crossings chronologically
    crossings.sort_by(|a, b| a.partial_cmp(b).unwrap());

    // Deduplicate crossings that are very close (~1 ms)
    const DEDUPE_EPS: f64 = 1e-8;
    crossings.dedup_by(|a, b| (a.value() - b.value()).abs() < DEDUPE_EPS);

    let is_above = |alt: f64| alt > threshold_rad;

    // Classify each crossing: +1 = entering above, -1 = exiting above
    // Probe ±10 seconds away from the crossing to determine direction.
    const PROBE_DT: Days = Days::new(10.0 / 86_400.0);
    let dt = PROBE_DT;
    let mut labeled: Vec<(JulianDate, i32)> = Vec::with_capacity(crossings.len());
    for &root in &crossings {
        let above_before = is_above(altitude_fn(root - dt));
        let above_after = is_above(altitude_fn(root + dt));

        if !above_before && above_after {
            labeled.push((root, 1)); // entering above
        } else if above_before && !above_after {
            labeled.push((root, -1)); // exiting above
        }
        // tangent (same side) → skip
    }

    let start_above = is_above(altitude_fn(jd_start));

    build_above_periods(&labeled, period, jd_start, jd_end, start_above, altitude_fn, threshold_rad)
}

/// Builds "above threshold" periods from pre-labelled crossings.
///
/// Separated from [`crossings_to_above_periods`] so that algorithms which
/// already know the crossing direction (e.g. lunar fast-scan) can skip the
/// classification step and avoid extra altitude evaluations.
///
/// # Arguments
///
/// * `labeled` – `(jd, direction)` pairs where `+1` = entering above,
///   `−1` = exiting above. Must be sorted chronologically.
/// * `period` – bounding search interval.
/// * `jd_start`, `jd_end` – Julian-Day equivalents of `period.start` / `.end`.
/// * `start_above` – whether the altitude at `jd_start` is above the threshold.
/// * `altitude_fn` – altitude evaluator (used only for midpoint validation).
/// * `threshold_rad` – threshold in radians.
pub(crate) fn build_above_periods<F>(
    labeled: &[(JulianDate, i32)],
    period: Period<ModifiedJulianDate>,
    jd_start: JulianDate,
    jd_end: JulianDate,
    start_above: bool,
    altitude_fn: &F,
    threshold_rad: f64,
) -> Vec<Period<ModifiedJulianDate>>
where
    F: Fn(JulianDate) -> f64,
{
    let is_above = |alt: f64| alt > threshold_rad;

    let mut periods: Vec<Period<ModifiedJulianDate>> = Vec::new();

    if labeled.is_empty() {
        if start_above {
            return vec![period];
        }
        return Vec::new();
    }

    let mut i = 0;

    // If we start above and first crossing is an exit, add initial interval
    if start_above && labeled[0].1 == -1 {
        let exit_mjd = ModifiedJulianDate::new(labeled[0].0.value() - 2400000.5);
        let mid = JulianDate::new((jd_start.value() + labeled[0].0.value()) * 0.5);
        if is_above(altitude_fn(mid)) {
            periods.push(Period::new(period.start, exit_mjd));
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

            let mid_jd = (enter_jd.value() + exit_mjd.to_julian_day().value()) * 0.5;
            if mid_jd >= jd_start.value() && mid_jd <= jd_end.value()
                && is_above(altitude_fn(JulianDate::new(mid_jd)))
            {
                periods.push(Period::new(enter_mjd, exit_mjd));
            }
        } else {
            i += 1;
        }
    }

    periods
}

// =============================================================================
// Public API
// =============================================================================

/// Finds all time periods where altitude is **above** the given threshold.
///
/// This is a generic, body-agnostic function that works with any altitude
/// evaluator closure. For Sun / Moon specific helpers, see the
/// [`crate::calculus::solar`] and [`crate::calculus::lunar`] modules.
///
/// To obtain *below*-threshold or *between*-range periods, compose with
/// [`crate::time::complement_within`] and [`crate::time::intersect_periods`].
///
/// # Arguments
///
/// * `altitude_fn` – returns altitude in **radians** at a given JD.
/// * `period` – the search interval.
/// * `threshold` – the altitude threshold.
///
/// # Algorithm
///
/// 1. Coarse scan at 10-minute intervals.
/// 2. Brent refinement for each sign-change bracket.
/// 3. Crossing classification and interval construction via
///    [`crossings_to_above_periods`].
pub fn find_above_altitude_periods<F>(
    altitude_fn: F,
    period: Period<ModifiedJulianDate>,
    threshold: Degrees,
) -> Vec<Period<ModifiedJulianDate>>
where
    F: Fn(JulianDate) -> f64,
{
    let jd_start = period.start.to_julian_day().value();
    let jd_end = period.end.to_julian_day().value();
    let threshold_rad = threshold.to::<Radian>().value();

    let step: Days = SCAN_STEP.to::<Day>();
    let f = |t: f64| altitude_fn(JulianDate::new(t));

    let raw = math_core::intervals::above_threshold_periods(
        jd_start, jd_end, step.value(), &f, threshold_rad,
    );

    raw.into_iter()
        .map(|iv| {
            Period::new(
                ModifiedJulianDate::new(iv.start - 2_400_000.5),
                ModifiedJulianDate::new(iv.end - 2_400_000.5),
            )
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_find_above_altitude_periods_synthetic() {
        // Synthetic altitude: a smooth sine wave over 1 day, shifted so that
        // the zero-crossings fall well inside the interval (not at endpoints).
        let altitude_fn = |jd: JulianDate| {
            let t = jd.value() - JulianDate::J2000.value();
            let phase = (t + 0.25) * core::f64::consts::TAU; // shift by ¼ day
            phase.sin()
        };

        let mjd_start = ModifiedJulianDate::new(JulianDate::J2000.value() - 2400000.5);
        let mjd_end = ModifiedJulianDate::new(mjd_start.value() + 1.0);

        let above = find_above_altitude_periods(
            altitude_fn,
            Period::new(mjd_start, mjd_end),
            Degrees::new(0.0),
        );
        assert!(!above.is_empty());
        for p in &above {
            assert!(p.duration_days() > 0.0);
        }

        // Complement should give "below 0" periods
        let below = crate::time::complement_within(
            Period::new(mjd_start, mjd_end),
            &above,
        );
        assert!(!below.is_empty());
        for p in &below {
            assert!(p.duration_days() > 0.0);
        }

        // Total coverage should sum to ~1.0 day
        let total: f64 = above.iter().chain(below.iter()).map(|p| p.duration_days()).sum();
        assert!((total - 1.0).abs() < 1e-6);
    }
}
