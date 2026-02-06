// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Altitude Window Periods
//!
//! Generic, body-agnostic tool for finding time intervals where a celestial
//! body's altitude is above a given threshold.
//!
//! ## Core Idea
//!
//! Every algorithm in the solar / lunar modules reduces to a single primitive:
//! **find periods where `altitude > threshold`**.
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
//! The public [`find_above_altitude_periods`] function delegates to
//! [`crate::calculus::math_core::intervals::above_threshold_periods`] with a
//! 10-minute coarse scan step, producing a "slow but safe" fallback.  The
//! body-specific modules (solar, lunar) use a 2-hour step for better
//! performance.
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
