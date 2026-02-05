// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Brent's method for root finding.
//!
//! Brent's method is a hybrid root-finding algorithm that combines the robustness
//! of bisection with the speed of inverse quadratic interpolation and the secant
//! method. It provides:
//!
//! - **Guaranteed convergence** when given a valid bracket (like bisection)
//! - **Superlinear convergence** in practice (faster than bisection)
//! - **Derivative-free**: only 1-2 function evaluations per iteration
//!
//! This is the **recommended method** for most astronomical applications where
//! brackets are available from sign-change scanning.
//!
//! ## Performance Advantage
//!
//! For expensive ephemeris calculations (VSOP87, ELP2000):
//! - Newton with finite-diff: ~3 evaluations/iteration
//! - Brent: ~1-2 evaluations/iteration
//! - **Typical speedup: 50-70% fewer total evaluations**

use crate::astro::JulianDate;
use qtty::Days;

// =============================================================================
// Constants
// =============================================================================

/// Brent's method tolerance — ~86 µs precision, balances speed and accuracy
const TOLERANCE: Days = Days::new(1e-9); // ~86 µs

/// Maximum Brent iterations (typically needs fewer than bisection)
const MAX_ITERS: usize = 100;

/// Function value convergence threshold
const CONVERGENCE_EPS: f64 = 1e-12;

// =============================================================================
// Implementation
// =============================================================================

/// Brent's method for finding a root in a bracketed interval.
///
/// This is a hybrid algorithm combining bisection, secant, and inverse quadratic
/// interpolation. It provides guaranteed convergence (like bisection) with
/// superlinear convergence in practice.
///
/// # Arguments
/// - `lo`, `hi`: Bracket endpoints where `scalar_fn(lo) - threshold` and
///   `scalar_fn(hi) - threshold` have opposite signs
/// - `scalar_fn`: Function to find root of
/// - `threshold`: Value to subtract from `scalar_fn` (finds where `scalar_fn(x) == threshold`)
///
/// # Returns
/// - `Some(jd)`: The root within tolerance `TOLERANCE` (≈ 0.86 µs)
/// - `None`: If bracket is invalid (same sign at both endpoints)
///
/// # Algorithm
///
/// At each iteration, Brent's method chooses between:
/// 1. **Inverse quadratic interpolation** — uses three points to fit a parabola
/// 2. **Secant method** — linear interpolation between two points
/// 3. **Bisection** — guaranteed progress when interpolation fails
///
/// The choice is made to balance speed with robustness, ensuring convergence
/// while maintaining superlinear convergence rates when possible.
pub fn refine_root<F>(
    lo: JulianDate,
    hi: JulianDate,
    scalar_fn: F,
    threshold: f64,
) -> Option<JulianDate>
where
    F: Fn(JulianDate) -> f64,
{
    // Evaluate at bracket endpoints
    let mut a = lo.value();
    let mut b = hi.value();
    let mut fa = scalar_fn(lo) - threshold;
    let mut fb = scalar_fn(hi) - threshold;

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

    let mut c = a; // Previous iterate
    let mut fc = fa;
    let mut d = b - a; // Step before last
    let mut e = d; // Last step taken

    for _ in 0..MAX_ITERS {
        // Convergence check
        let tol = TOLERANCE.value();
        let m = 0.5 * (c - b); // Midpoint offset

        if fb.abs() < CONVERGENCE_EPS || m.abs() <= tol {
            return Some(JulianDate::new(b));
        }

        // Decide whether to use bisection or interpolation
        let use_bisection = e.abs() < tol || fa.abs() <= fb.abs();

        let (new_e, new_d) = if use_bisection {
            // Bisection step
            (m, m)
        } else {
            // Try interpolation
            let s = fb / fa;

            let (p, q) = if (a - c).abs() < 1e-14 {
                // Secant method (linear interpolation)
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

            // Ensure p/q has correct sign
            let (p, q) = if p > 0.0 { (p, -q) } else { (-p, q) };

            // Accept interpolation only if it's better than bisection
            let s_val = e;
            if 2.0 * p < 3.0 * m * q - (tol * q).abs() && p < (0.5 * s_val * q).abs() {
                (d, p / q)
            } else {
                // Reject interpolation, use bisection
                (m, m)
            }
        };

        e = new_e;
        d = new_d;

        // Update previous point
        a = b;
        fa = fb;

        // Compute new point
        if d.abs() > tol {
            b += d;
        } else {
            // Ensure minimum step toward root
            b += if m > 0.0 { tol } else { -tol };
        }
        fb = scalar_fn(JulianDate::new(b)) - threshold;

        // Keep bracket: ensure f(b) and f(c) have opposite signs
        if fb * fc > 0.0 {
            c = a;
            fc = fa;
            e = b - a;
            d = e;
        }

        // Ensure |f(c)| >= |f(b)|
        if fc.abs() < fb.abs() {
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
        }
    }

    // Return best estimate even if max iterations reached
    Some(JulianDate::new(b))
}

/// Brent's method variant that accepts pre-computed function values at bracket endpoints.
///
/// This is **more efficient** when the caller has already evaluated the function at the
/// bracket endpoints (e.g., during a sign-change scan). Avoids 2 redundant function evaluations.
///
/// # Arguments
/// - `lo`, `hi`: Bracket endpoints
/// - `f_lo`, `f_hi`: Pre-computed `scalar_fn(lo) - threshold` and `scalar_fn(hi) - threshold`
/// - `scalar_fn`: Function to find root of
/// - `threshold`: Value to subtract from `scalar_fn` (finds where `scalar_fn(x) == threshold`)
///
/// # Returns
/// - `Some(jd)`: The root within tolerance `TOLERANCE` (≈ 0.86 µs)
/// - `None`: If bracket is invalid (same sign at both endpoints)
pub fn refine_root_with_values<F>(
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
        let tol = TOLERANCE.value();
        let m = 0.5 * (c - b);

        if fb.abs() < CONVERGENCE_EPS || m.abs() <= tol {
            return Some(JulianDate::new(b));
        }

        let use_bisection = e.abs() < tol || fa.abs() <= fb.abs();

        let (new_e, new_d) = if use_bisection {
            (m, m)
        } else {
            let s = fb / fa;

            let (p, q) = if (a - c).abs() < 1e-14 {
                let p = 2.0 * m * s;
                let q = 1.0 - s;
                (p, q)
            } else {
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
        } else {
            b += if m > 0.0 { tol } else { -tol };
        }
        fb = scalar_fn(JulianDate::new(b)) - threshold;

        if fb * fc > 0.0 {
            c = a;
            fc = fa;
            e = b - a;
            d = e;
        }

        if fc.abs() < fb.abs() {
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
        }
    }

    Some(JulianDate::new(b))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::cell::Cell;

    #[test]
    fn finds_sine_root_near_pi() {
        let lo = JulianDate::new(3.0);
        let hi = JulianDate::new(4.0);
        let scalar = |jd: JulianDate| jd.value().sin();

        let result = refine_root(lo, hi, scalar, 0.0).expect("should find pi");
        assert!(
            (result.value() - std::f64::consts::PI).abs() < 1e-10,
            "expected root near π, got {}",
            result.value()
        );
    }

    #[test]
    fn refines_linear_root() {
        let threshold = 2_458_850.0;
        let lo = JulianDate::new(threshold - 1.0);
        let hi = JulianDate::new(threshold + 1.0);
        let scalar = |jd: JulianDate| jd.value();

        let found = refine_root(lo, hi, scalar, threshold).expect("Brent should converge");
        assert!((found.value() - threshold).abs() < 1e-10);
    }

    #[test]
    fn returns_lo_when_endpoint_is_root() {
        let threshold = 12_345.678;
        let lo = JulianDate::new(threshold);
        let hi = JulianDate::new(threshold + 4.0);
        let scalar = |jd: JulianDate| jd.value();

        let result = refine_root(lo, hi, scalar, threshold);
        assert_eq!(result, Some(lo));
    }

    #[test]
    fn returns_hi_when_endpoint_is_root() {
        let threshold = 54_321.0;
        let lo = JulianDate::new(threshold - 4.0);
        let hi = JulianDate::new(threshold);
        let scalar = |jd: JulianDate| jd.value();

        let result = refine_root(lo, hi, scalar, threshold);
        assert_eq!(result, Some(hi));
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

    #[test]
    fn uses_fewer_evaluations_than_newton_for_trig() {
        // Brent should need ~1-2 evals/iter vs Newton's ~3 evals/iter
        let brent_calls = Cell::new(0usize);
        let newton_calls = Cell::new(0usize);

        let lo = JulianDate::new(3.0);
        let hi = JulianDate::new(4.0);

        // Brent
        let scalar_brent = |jd: JulianDate| {
            brent_calls.set(brent_calls.get() + 1);
            jd.value().sin()
        };
        let _ = refine_root(lo, hi, scalar_brent, 0.0);

        // Newton
        let guess = JulianDate::new(3.5);
        let scalar_newton = |jd: JulianDate| {
            newton_calls.set(newton_calls.get() + 1);
            jd.value().sin()
        };
        let _ = super::super::newton::refine_root(guess, scalar_newton, 0.0);

        // Brent should use fewer function evaluations
        assert!(
            brent_calls.get() < newton_calls.get(),
            "Brent used {} calls, Newton used {} calls",
            brent_calls.get(),
            newton_calls.get()
        );
    }

    #[test]
    fn handles_quadratic_with_two_roots() {
        // f(x) = x^2 - 4, roots at x = ±2
        let lo = JulianDate::new(1.0);
        let hi = JulianDate::new(3.0);
        let scalar = |jd: JulianDate| jd.value().powi(2);

        let result = refine_root(lo, hi, scalar, 4.0).expect("should find x=2");
        assert!((result.value() - 2.0).abs() < 1e-9);
    }

    #[test]
    fn converges_for_cubic() {
        // f(x) = x^3 - 2, root at x = 2^(1/3) ≈ 1.2599
        let lo = JulianDate::new(1.0);
        let hi = JulianDate::new(2.0);
        let scalar = |jd: JulianDate| jd.value().powi(3);

        let result = refine_root(lo, hi, scalar, 2.0).expect("should find cube root of 2");
        let expected = 2.0_f64.powf(1.0 / 3.0);
        assert!((result.value() - expected).abs() < 1e-9);
    }
}
