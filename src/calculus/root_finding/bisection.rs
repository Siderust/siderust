// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Bisection method for root finding.
//!
//! The bisection method is a simple, robust bracketing algorithm that repeatedly
//! halves an interval known to contain a root. It guarantees convergence but has
//! only linear convergence (interval width decreases by factor of 2 per iteration).
//!
//! **Use this method when:**
//! - You need guaranteed convergence regardless of function behavior
//! - The function is discontinuous or has other pathologies
//! - You want a simple, predictable fallback
//!
//! **Consider alternatives:**
//! - [`super::brent::refine_root`] — Usually faster while maintaining robustness
//! - [`super::newton::refine_root`] — Much faster when you have a good initial guess

use crate::astro::JulianDate;
use qtty::Days;

// =============================================================================
// Constants
// =============================================================================

/// Bisection tolerance (≈ 0.086 ms)
const TOLERANCE: Days = Days::new(1e-9);

/// Maximum bisection iterations
const MAX_ITERS: usize = 80;

/// Function value convergence threshold
const CONVERGENCE_EPS: f64 = 1e-12;

// =============================================================================
// Implementation
// =============================================================================

/// Bisection fallback for a scalar root between two Julian dates.
/// 
/// Requires `scalar_fn(lo) - threshold` and `scalar_fn(hi) - threshold` to have
/// opposite signs (valid bracket).
///
/// # Algorithm
///
/// Repeatedly bisects the interval [lo, hi] and updates the bracket based on
/// the sign of the function at the midpoint:
///
/// ```text
/// mid = (lo + hi) / 2
/// if sign(f(mid)) == sign(f(lo)):
///     lo = mid
/// else:
///     hi = mid
/// ```
///
/// Terminates when the interval width is below tolerance or the function value
/// is sufficiently close to zero.
pub fn refine_root<F>(
    mut lo: JulianDate,
    mut hi: JulianDate,
    scalar_fn: F,
    threshold: f64,
) -> Option<JulianDate>
where
    F: Fn(JulianDate) -> f64,
{
    let mut f_lo = scalar_fn(lo) - threshold;
    let f_hi = scalar_fn(hi) - threshold;

    // Check endpoints for exact roots
    if f_lo.abs() < CONVERGENCE_EPS {
        return Some(lo);
    }
    if f_hi.abs() < CONVERGENCE_EPS {
        return Some(hi);
    }

    // Verify bracket (must have opposite signs)
    if f_lo * f_hi > 0.0 {
        return None;
    }

    for _ in 0..MAX_ITERS {
        let mid = JulianDate::new((lo.value() + hi.value()) * 0.5);
        let f_mid = scalar_fn(mid) - threshold;
        let interval_width = Days::new((hi.value() - lo.value()).abs());

        // Check convergence: function value or interval width
        if f_mid.abs() < CONVERGENCE_EPS || interval_width < TOLERANCE {
            return Some(mid);
        }

        // Update bracket based on sign of f_mid
        if f_lo * f_mid <= 0.0 {
            hi = mid;
        } else {
            lo = mid;
            f_lo = f_mid;
        }
    }

    // Return midpoint if max iterations reached
    Some(JulianDate::new((lo.value() + hi.value()) * 0.5))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::cell::Cell;

    #[test]
    fn returns_lo_when_endpoint_matches_root() {
        let threshold = 12_345.678;
        let lo = JulianDate::new(threshold);
        let hi = JulianDate::new(threshold + 4.0);
        let scalar = |jd: JulianDate| jd.value();

        assert_eq!(refine_root(lo, hi, scalar, threshold), Some(lo));
    }

    #[test]
    fn returns_hi_when_endpoint_matches_root() {
        let threshold = 54_321.0;
        let lo = JulianDate::new(threshold - 4.0);
        let hi = JulianDate::new(threshold);
        let scalar = |jd: JulianDate| jd.value();

        assert_eq!(refine_root(lo, hi, scalar, threshold), Some(hi));
    }

    #[test]
    fn finds_root_between_brackets() {
        let threshold = 42.0;
        let lo = JulianDate::new(threshold - 1.0);
        let hi = JulianDate::new(threshold + 1.0);
        let scalar = |jd: JulianDate| jd.value();

        let found = refine_root(lo, hi, scalar, threshold).expect("bisection failed");
        assert!((found.value() - threshold).abs() < 1e-9);
    }

    #[test]
    fn returns_midpoint_after_max_iterations() {
        let calls = Cell::new(0usize);
        let threshold = 1.0e15 + 1.0;
        let lo = JulianDate::new(0.0);
        let hi = JulianDate::new(2.0e15);
        let scalar = |jd: JulianDate| {
            calls.set(calls.get() + 1);
            jd.value()
        };

        let found = refine_root(lo, hi, scalar, threshold).expect("bisection failed");
        // Ensure we iterated at least once and terminated early due to tolerance
        // being tighter than the floating precision at this magnitude.
        assert!(calls.get() > 2 && calls.get() <= MAX_ITERS + 2);
        assert!(found.value().is_finite());
    }

    #[test]
    fn returns_none_for_invalid_bracket() {
        let scalar = |_jd: JulianDate| 99.0;
        let lo = JulianDate::new(0.0);
        let hi = JulianDate::new(1.0);

        assert!(refine_root(lo, hi, scalar, 0.0).is_none());
    }

    #[test]
    fn handles_step_function() {
        let lo = JulianDate::new(-1.0);
        let hi = JulianDate::new(1.0);
        let scalar = |jd: JulianDate| {
            if jd.value() < 0.0 {
                -5.0
            } else {
                5.0
            }
        };

        let result = refine_root(lo, hi, scalar, 0.0).expect("should handle step function");
        assert!(result.value().abs() < 1e-6);
    }
}
