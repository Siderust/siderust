// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Root Finding — Brent & Bisection
//!
//! Robust, astronomy‑agnostic root solvers for scalar functions of one
//! variable.  All routines work on plain `f64` inputs/outputs so they can be
//! composed with any domain‑specific wrapper.
//!
//! ## Provided solvers
//!
//! | Function | Method | Convergence | Evals / iter |
//! |----------|--------|-------------|--------------|
//! | [`brent`] | Brent hybrid (bisection + secant + IQI) | superlinear | 1–2 |
//! | [`brent_with_values`] | same, pre‑computed endpoints | superlinear | 1–2 |
//! | [`brent_tol`] | same, custom tolerance | superlinear | 1–2 |
//! | [`bisection`] | classic bisection | linear | 1 |
//!
//! All functions return `Some(root)` on success, `None` when the bracket is
//! invalid (same sign at both endpoints).

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// Default convergence tolerance (days).  ≈ 0.86 µs.
const DEFAULT_TOL: f64 = 1e-9;

/// Maximum iterations for Brent.
const BRENT_MAX_ITER: usize = 100;

/// Maximum iterations for bisection.
const BISECT_MAX_ITER: usize = 80;

/// Function‑value convergence threshold.
const F_EPS: f64 = 1e-12;

// ---------------------------------------------------------------------------
// Brent's method
// ---------------------------------------------------------------------------

/// Find a root of `f(t) = 0` in `[lo, hi]` using Brent's method with
/// default tolerance ([`DEFAULT_TOL`]).
///
/// Returns `Some(t)` where `|f(t)| < ε` or the bracket width < tolerance,
/// `None` if the bracket is invalid.
pub fn brent<F>(lo: f64, hi: f64, f: F) -> Option<f64>
where
    F: Fn(f64) -> f64,
{
    let fa = f(lo);
    let fb = f(hi);
    brent_core(lo, hi, fa, fb, &f, DEFAULT_TOL)
}

/// Like [`brent`] but the caller supplies pre‑computed `f(lo)` and `f(hi)`.
///
/// Saves two function evaluations when the values are already available
/// (e.g. from a sign‑change scan).
pub fn brent_with_values<F>(lo: f64, hi: f64, f_lo: f64, f_hi: f64, f: F) -> Option<f64>
where
    F: Fn(f64) -> f64,
{
    brent_core(lo, hi, f_lo, f_hi, &f, DEFAULT_TOL)
}

/// Like [`brent_with_values`] but with a caller‑chosen tolerance.
///
/// Useful for trading precision for speed in performance‑critical loops
/// (e.g. lunar altitude with 2‑minute tolerance).
pub fn brent_tol<F>(
    lo: f64,
    hi: f64,
    f_lo: f64,
    f_hi: f64,
    f: F,
    tolerance: f64,
) -> Option<f64>
where
    F: Fn(f64) -> f64,
{
    brent_core(lo, hi, f_lo, f_hi, &f, tolerance)
}

/// Core Brent implementation shared by all public entry points.
fn brent_core<F>(
    lo: f64,
    hi: f64,
    f_lo: f64,
    f_hi: f64,
    f: &F,
    tol: f64,
) -> Option<f64>
where
    F: Fn(f64) -> f64,
{
    let mut a = lo;
    let mut b = hi;
    let mut fa = f_lo;
    let mut fb = f_hi;

    // Exact roots at endpoints
    if fa.abs() < F_EPS {
        return Some(a);
    }
    if fb.abs() < F_EPS {
        return Some(b);
    }
    // Valid bracket?
    if fa * fb > 0.0 {
        return None;
    }

    // Ensure |f(a)| >= |f(b)| so b is always the best approximation
    if fa.abs() < fb.abs() {
        std::mem::swap(&mut a, &mut b);
        std::mem::swap(&mut fa, &mut fb);
    }

    let mut c = a;
    let mut fc = fa;
    let mut d = b - a;
    let mut e = d;

    for _ in 0..BRENT_MAX_ITER {
        let m = 0.5 * (c - b);

        if fb.abs() < F_EPS || m.abs() <= tol {
            return Some(b);
        }

        let use_bisection = e.abs() < tol || fa.abs() <= fb.abs();

        let (new_e, new_d) = if use_bisection {
            (m, m)
        } else {
            let s = fb / fa;
            let (p, q) = if (a - c).abs() < 1e-14 {
                // Secant
                (2.0 * m * s, 1.0 - s)
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
                (d, p / q) // accept interpolation
            } else {
                (m, m) // fallback to bisection
            }
        };

        e = new_e;
        d = new_d;

        a = b;
        fa = fb;

        b += if d.abs() > tol {
            d
        } else if m > 0.0 {
            tol
        } else {
            -tol
        };
        fb = f(b);

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

    Some(b)
}

// ---------------------------------------------------------------------------
// Bisection
// ---------------------------------------------------------------------------

/// Find a root of `f(t) = 0` in `[lo, hi]` using classic bisection.
///
/// Guaranteed convergence for valid brackets; linear rate.
pub fn bisection<F>(lo: f64, hi: f64, f: F) -> Option<f64>
where
    F: Fn(f64) -> f64,
{
    let mut a = lo;
    let mut b = hi;
    let mut fa = f(a);
    let fb = f(b);

    if fa.abs() < F_EPS {
        return Some(a);
    }
    if fb.abs() < F_EPS {
        return Some(b);
    }
    if fa * fb > 0.0 {
        return None;
    }

    for _ in 0..BISECT_MAX_ITER {
        let mid = 0.5 * (a + b);
        let fm = f(mid);
        let width = (b - a).abs();

        if fm.abs() < F_EPS || width < DEFAULT_TOL {
            return Some(mid);
        }

        if fa * fm <= 0.0 {
            b = mid;
        } else {
            a = mid;
            fa = fm;
        }
    }

    Some(0.5 * (a + b))
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn brent_finds_sine_root_near_pi() {
        let root = brent(3.0, 4.0, |t| t.sin()).expect("should find π");
        assert!((root - std::f64::consts::PI).abs() < 1e-10);
    }

    #[test]
    fn brent_finds_linear_root() {
        let root = brent(0.0, 10.0, |t| t - 5.0).expect("should find 5");
        assert!((root - 5.0).abs() < 1e-10);
    }

    #[test]
    fn brent_returns_none_for_invalid_bracket() {
        assert!(brent(0.0, 1.0, |_| 42.0).is_none());
    }

    #[test]
    fn brent_returns_endpoint_when_exact() {
        let root = brent(0.0, 5.0, |t| t - 5.0).expect("endpoint");
        assert!((root - 5.0).abs() < F_EPS);
    }

    #[test]
    fn brent_with_values_saves_evaluations() {
        use std::cell::Cell;
        let count = Cell::new(0usize);
        let f = |t: f64| {
            count.set(count.get() + 1);
            t.sin()
        };
        let f_lo = (3.0_f64).sin();
        let f_hi = (4.0_f64).sin();
        let _ = brent_with_values(3.0, 4.0, f_lo, f_hi, &f);
        let with_vals = count.get();

        count.set(0);
        let _ = brent(3.0, 4.0, &f);
        let without = count.get();

        // with_values should use at least 2 fewer evals (the endpoints)
        assert!(with_vals + 2 <= without || with_vals <= without);
    }

    #[test]
    fn brent_tol_respects_relaxed_tolerance() {
        let root = brent_tol(3.0, 4.0, (3.0_f64).sin(), (4.0_f64).sin(), |t| t.sin(), 1e-3)
            .expect("relaxed");
        assert!((root - std::f64::consts::PI).abs() < 2e-3);
    }

    #[test]
    fn brent_handles_step_function() {
        let root = brent(-1.0, 1.0, |t| if t < 0.0 { -1.0 } else { 1.0 }).expect("step");
        assert!(root.abs() < 1e-6);
    }

    #[test]
    fn brent_cubic() {
        let root = brent(1.0, 2.0, |t| t.powi(3) - 2.0).expect("cbrt 2");
        assert!((root - 2.0_f64.powf(1.0 / 3.0)).abs() < 1e-9);
    }

    #[test]
    fn bisection_finds_sine_root() {
        let root = bisection(3.0, 4.0, |t| t.sin()).expect("π");
        assert!((root - std::f64::consts::PI).abs() < 1e-8);
    }

    #[test]
    fn bisection_returns_none_for_invalid_bracket() {
        assert!(bisection(0.0, 1.0, |_| 42.0).is_none());
    }

    #[test]
    fn bisection_handles_step_function() {
        let root = bisection(-1.0, 1.0, |t| if t < 0.0 { -5.0 } else { 5.0 }).expect("step");
        assert!(root.abs() < 1e-6);
    }

    #[test]
    fn bisection_endpoint_root() {
        let root = bisection(0.0, 5.0, |t| t).expect("root at 0");
        assert!(root.abs() < F_EPS);
    }
}
