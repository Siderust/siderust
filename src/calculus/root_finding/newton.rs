// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Newton-Raphson method with finite-difference derivatives.
//!
//! Newton-Raphson provides quadratic convergence near the root, but requires
//! computing derivatives. Since analytical derivatives are not available for
//! general astronomical functions, we use symmetric finite differences:
//!
//! ```text
//! f'(x) ≈ (f(x+h) - f(x-h)) / (2h)
//! ```
//!
//! This requires **3 function evaluations per iteration** (f(x-h), f(x), f(x+h)),
//! which can be expensive for VSOP87/ELP2000 calculations. Consider using
//! Brent's method when brackets are available.

use crate::astro::JulianDate;
use qtty::Days;

use super::DERIVATIVE_STEP;

// =============================================================================
// Constants
// =============================================================================

/// Finite-difference step for derivative estimation (1 second)
const FD_STEP: Days = DERIVATIVE_STEP;

/// Newton-Raphson tolerance (≈ 0.86 µs)
const TOLERANCE: Days = Days::new(1e-11);

/// Maximum Newton iterations
const MAX_ITERS: usize = 20;

/// Function value convergence threshold
const CONVERGENCE_EPS: f64 = 1e-12;

/// Minimum acceptable derivative value
const MIN_DERIVATIVE: f64 = 1e-14;

// =============================================================================
// Implementation
// =============================================================================

/// Newton-Raphson refinement for a scalar root of `f(x) - threshold = 0`.
///
/// This is a generic routine that takes an initial Julian date guess `jd`, a
/// scalar function `scalar_fn` evaluated at a given `JulianDate`, and the
/// `threshold` to subtract. Returns `Some(jd)` if a root is found, `None` on
/// failure.
///
/// # Algorithm
///
/// Iteratively refines the guess using:
/// ```text
/// x_{n+1} = x_n - f(x_n) / f'(x_n)
/// ```
///
/// where `f'(x_n)` is estimated via symmetric finite differences.
///
/// # Performance
///
/// Each iteration requires **3 function evaluations**: f(x-h), f(x), f(x+h).
/// For expensive ephemeris calculations, this can dominate runtime. Consider
/// [`super::brent::refine_root`] when valid brackets are available.
pub fn refine_root<F>(mut jd: JulianDate, scalar_fn: F, threshold: f64) -> Option<JulianDate>
where
    F: Fn(JulianDate) -> f64,
{
    for _ in 0..MAX_ITERS {
        let f = scalar_fn(jd) - threshold;

        // Check if already converged
        if f.abs() < CONVERGENCE_EPS {
            return Some(jd);
        }

        // Symmetric finite-difference derivative: f'(x) ≈ (f(x+h) - f(x-h)) / 2h
        let f_plus = scalar_fn(jd + FD_STEP);
        let f_minus = scalar_fn(jd - FD_STEP);
        let deriv = (f_plus - f_minus) / (2.0 * FD_STEP.value()); // units/day

        // Protect against near-zero derivative
        if deriv.abs() < MIN_DERIVATIVE {
            return None;
        }

        // Newton step: x_new = x - f(x)/f'(x)
        let delta = Days::new(f / deriv);
        jd -= delta;

        if delta.abs() < TOLERANCE {
            return Some(jd);
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::cell::Cell;

    #[test]
    fn returns_guess_when_root_already_exact() {
        let threshold = 25_000.0;
        let guess = JulianDate::new(threshold);

        let scalar = |jd: JulianDate| jd.value();
        let found = refine_root(guess, scalar, threshold);

        assert_eq!(found, Some(guess));
    }

    #[test]
    fn refines_linear_root() {
        let threshold = 2_458_850.0;
        let guess = JulianDate::new(threshold + 0.5);
        let scalar = |jd: JulianDate| jd.value();

        let found = refine_root(guess, scalar, threshold).expect("Newton should converge");
        assert!((found.value() - threshold).abs() < 1e-9);
    }

    #[test]
    fn returns_none_for_constant_function() {
        let threshold = 0.0;
        let guess = JulianDate::J2000;
        let scalar = |_jd: JulianDate| std::f64::consts::PI;

        assert!(refine_root(guess, scalar, threshold).is_none());
    }

    #[test]
    fn returns_none_after_max_iterations() {
        let calls = Cell::new(0usize);
        let scalar = |_jd: JulianDate| {
            let count = calls.get() + 1;
            calls.set(count);
            count as f64
        };

        let guess = JulianDate::J2000;
        assert!(refine_root(guess, scalar, 0.0).is_none());
    }

    #[test]
    fn finds_sine_root() {
        let guess = JulianDate::new(3.5);
        let scalar = |jd: JulianDate| jd.value().sin();

        let result = refine_root(guess, scalar, 0.0).expect("should find pi");
        assert!((result.value() - std::f64::consts::PI).abs() < 1e-9);
    }
}
