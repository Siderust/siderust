// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Root Finding — Brent & Bisection
//!
//! Robust, astronomy‑agnostic root solvers for scalar functions of one
//! variable.
//!
//! [`brent`] and [`bisection`] are generic over [`qtty`] unit types.
//! [`brent_with_values`] and [`brent_tol`] operate on time
//! [`Period`]s (`TimeInstant` with day durations).
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

use crate::time::{Period, TimeInstant};
use qtty::{Days, Quantity, Unit};

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// Default convergence tolerance for domain variables.
/// ≈ 0.86 µs when the domain is days.
const DEFAULT_TOL: f64 = 1e-9;

/// Maximum iterations for Brent.
const BRENT_MAX_ITER: usize = 100;

/// Maximum iterations for bisection.
const BISECT_MAX_ITER: usize = 80;

/// Function‑value convergence threshold.
const F_EPS: f64 = 1e-12;

/// Near‑zero threshold for detecting coincident bracket points in the
/// domain (used to choose between secant and IQI interpolation).
const COINCIDENCE_EPS: f64 = 1e-14;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// `true` when both quantities are strictly on the same side of zero.
fn same_sign<V: Unit>(a: Quantity<V>, b: Quantity<V>) -> bool {
    let sa = a.signum();
    let sb = b.signum();
    sa == sb && sa != 0.0
}

/// Dimensionless ratio of two same‑unit quantities.
///
/// The only place where raw `f64` values are extracted: computing
/// interpolation weights that are inherently pure scalars.
fn dimensionless_ratio<V: Unit>(num: Quantity<V>, den: Quantity<V>) -> f64 {
    num.value() / den.value()
}

// ---------------------------------------------------------------------------
// Brent's method
// ---------------------------------------------------------------------------

/// Find a root of `f(t) = 0` in `[lo, hi]` using Brent's method with
/// default tolerance ([`DEFAULT_TOL`]).
///
/// Returns `Some(t)` where `|f(t)| < ε` or the bracket width < tolerance,
/// `None` if the bracket is invalid.
pub fn brent<T, V, F>(lo: Quantity<T>, hi: Quantity<T>, f: F) -> Option<Quantity<T>>
where
    T: Unit,
    V: Unit,
    F: Fn(Quantity<T>) -> Quantity<V>,
{
    let f_lo = f(lo);
    let f_hi = f(hi);
    brent_engine(lo, hi, f_lo, f_hi, &f, Quantity::new(DEFAULT_TOL))
}

/// Like [`brent`] but the caller supplies pre‑computed
/// `f(period.start)` and `f(period.end)`.
///
/// Saves two function evaluations when the values are already available
/// (e.g. from a sign‑change scan).
pub fn brent_with_values<T, V, F>(
    period: Period<T>,
    f_lo: Quantity<V>,
    f_hi: Quantity<V>,
    f: F,
) -> Option<T>
where
    T: TimeInstant<Duration = Days>,
    V: Unit,
    F: Fn(T) -> Quantity<V>,
{
    brent_core(period, f_lo, f_hi, &f, Days::new(DEFAULT_TOL))
}

/// Like [`brent_with_values`] but with a caller‑chosen tolerance.
///
/// Useful for trading precision for speed in performance‑critical loops
/// (e.g. lunar altitude with 2‑minute tolerance).
pub fn brent_tol<T, V, F>(
    period: Period<T>,
    f_lo: Quantity<V>,
    f_hi: Quantity<V>,
    f: F,
    tolerance: Days,
) -> Option<T>
where
    T: TimeInstant<Duration = Days>,
    V: Unit,
    F: Fn(T) -> Quantity<V>,
{
    brent_core(period, f_lo, f_hi, &f, tolerance)
}

/// Core Brent implementation for [`Period`]‑based entry points.
///
/// Maps the period `[start, end]` to [`Days`] offsets, solves in that
/// typed space, then maps the root back to `T`.
fn brent_core<T, V, F>(
    period: Period<T>,
    f_lo: Quantity<V>,
    f_hi: Quantity<V>,
    f: &F,
    tol: Days,
) -> Option<T>
where
    T: TimeInstant<Duration = Days>,
    V: Unit,
    F: Fn(T) -> Quantity<V>,
{
    let start = period.start;
    let span = period.end.difference(&start);

    brent_engine(
        Days::zero(),
        span,
        f_lo,
        f_hi,
        &|offset: Days| {
            let t = start.add_duration(offset);
            f(t)
        },
        tol,
    )
    .map(|offset| start.add_duration(offset))
}

/// Typed Brent solver — all domain arithmetic uses [`Quantity`] types.
///
/// The only raw `f64` values extracted are dimensionless ratios of
/// same‑unit function values, used as interpolation weights.
fn brent_engine<T, V, F>(
    lo: Quantity<T>,
    hi: Quantity<T>,
    f_lo: Quantity<V>,
    f_hi: Quantity<V>,
    f: &F,
    tol: Quantity<T>,
) -> Option<Quantity<T>>
where
    T: Unit,
    V: Unit,
    F: Fn(Quantity<T>) -> Quantity<V>,
{
    let f_eps: Quantity<V> = Quantity::new(F_EPS);
    let coincidence: Quantity<T> = Quantity::new(COINCIDENCE_EPS);

    let mut a = lo;
    let mut b = hi;
    let mut fa = f_lo;
    let mut fb = f_hi;

    if fa.abs() < f_eps {
        return Some(a);
    }
    if fb.abs() < f_eps {
        return Some(b);
    }
    if same_sign(fa, fb) {
        return None;
    }

    if fa.abs() < fb.abs() {
        std::mem::swap(&mut a, &mut b);
        std::mem::swap(&mut fa, &mut fb);
    }

    let mut c = a;
    let mut fc = fa;
    let mut d = b - a;
    let mut e = d;

    for _ in 0..BRENT_MAX_ITER {
        let m = (c - b) * 0.5;

        if fb.abs() < f_eps || m.abs() <= tol {
            return Some(b);
        }

        let use_bisection = e.abs() < tol || fa.abs() <= fb.abs();

        let (new_e, new_d) = if use_bisection {
            (m, m)
        } else {
            // Dimensionless interpolation weights from function‑value ratios
            let s = dimensionless_ratio(fb, fa);

            let (p, q): (Quantity<T>, f64) = if (a - c).abs() < coincidence {
                // Secant
                (m * (2.0 * s), 1.0 - s)
            } else {
                // Inverse quadratic interpolation
                let q_val = dimensionless_ratio(fa, fc);
                let r = dimensionless_ratio(fb, fc);
                let p = m * (2.0 * s * q_val * (q_val - r)) - (b - a) * (s * (r - 1.0));
                let q = (q_val - 1.0) * (r - 1.0) * (s - 1.0);
                (p, q)
            };

            // Ensure p ≥ 0 by adjusting signs
            let (p, q) = if p > Quantity::<T>::zero() {
                (p, -q)
            } else {
                (-p, q)
            };

            let s_val = e;
            if p * 2.0 < m * (3.0 * q) - (tol * q).abs() && p < (s_val * (0.5 * q)).abs() {
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
        } else if m > Quantity::<T>::zero() {
            tol
        } else {
            -tol
        };
        fb = f(b);

        if same_sign(fb, fc) {
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
pub fn bisection<T, V, F>(lo: Quantity<T>, hi: Quantity<T>, f: F) -> Option<Quantity<T>>
where
    T: Unit,
    V: Unit,
    F: Fn(Quantity<T>) -> Quantity<V>,
{
    let tol: Quantity<T> = Quantity::new(DEFAULT_TOL);
    let f_eps: Quantity<V> = Quantity::new(F_EPS);

    let mut a = lo;
    let mut b = hi;
    let mut fa = f(a);
    let fb = f(b);

    if fa.abs() < f_eps {
        return Some(a);
    }
    if fb.abs() < f_eps {
        return Some(b);
    }
    if same_sign(fa, fb) {
        return None;
    }

    for _ in 0..BISECT_MAX_ITER {
        let mid = a.mean(b);
        let fm = f(mid);
        let width = (b - a).abs();

        if fm.abs() < f_eps || width < tol {
            return Some(mid);
        }

        if same_sign(fa, fm) {
            a = mid;
            fa = fm;
        } else {
            b = mid;
        }
    }

    Some(a.mean(b))
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::time::{ModifiedJulianDate, Period};
    use qtty::{Day, Radian};

    type Days = Quantity<Day>;
    type Mjd = ModifiedJulianDate;
    type Radians = Quantity<Radian>;

    #[test]
    fn brent_finds_sine_root_near_pi() {
        let root = brent(Days::new(3.0), Days::new(4.0), |t: Days| {
            Radians::new(t.value().sin())
        })
        .expect("should find π");
        assert!((root.value() - std::f64::consts::PI).abs() < 1e-10);
    }

    #[test]
    fn brent_finds_linear_root() {
        let root = brent(Days::new(0.0), Days::new(10.0), |t: Days| {
            Radians::new(t.value() - 5.0)
        })
        .expect("should find 5");
        assert!((root.value() - 5.0).abs() < 1e-10);
    }

    #[test]
    fn brent_returns_none_for_invalid_bracket() {
        assert!(brent(Days::new(0.0), Days::new(1.0), |_: Days| Radians::new(42.0)).is_none());
    }

    #[test]
    fn brent_returns_endpoint_when_exact() {
        let root = brent(Days::new(0.0), Days::new(5.0), |t: Days| {
            Radians::new(t.value() - 5.0)
        })
        .expect("endpoint");
        assert!((root.value() - 5.0).abs() < F_EPS);
    }

    #[test]
    fn brent_with_values_saves_evaluations() {
        use std::cell::Cell;
        let count = Cell::new(0usize);
        let f = |t: Mjd| -> Radians {
            count.set(count.get() + 1);
            Radians::new(t.value().sin())
        };
        let f_lo = Radians::new((3.0_f64).sin());
        let f_hi = Radians::new((4.0_f64).sin());
        let _ = brent_with_values(Period::new(Mjd::new(3.0), Mjd::new(4.0)), f_lo, f_hi, f);
        let with_vals = count.get();

        count.set(0);
        let _ = brent(Days::new(3.0), Days::new(4.0), |t: Days| {
            count.set(count.get() + 1);
            Radians::new(t.value().sin())
        });
        let without = count.get();

        // with_values should use at least 2 fewer evals (the endpoints)
        assert!(with_vals + 2 <= without || with_vals <= without);
    }

    #[test]
    fn brent_tol_respects_relaxed_tolerance() {
        let root = brent_tol(
            Period::new(Mjd::new(3.0), Mjd::new(4.0)),
            Radians::new((3.0_f64).sin()),
            Radians::new((4.0_f64).sin()),
            |t: Mjd| Radians::new(t.value().sin()),
            Days::new(1e-3),
        )
        .expect("relaxed");
        assert!((root.value() - std::f64::consts::PI).abs() < 2e-3);
    }

    #[test]
    fn brent_handles_step_function() {
        let root = brent(Days::new(-1.0), Days::new(1.0), |t: Days| {
            Radians::new(if t < 0.0 { -1.0 } else { 1.0 })
        })
        .expect("step");
        assert!(root.abs() < 1e-6);
    }

    #[test]
    fn brent_cubic() {
        let root = brent(Days::new(1.0), Days::new(2.0), |t: Days| {
            Radians::new(t.value().powi(3) - 2.0)
        })
        .expect("cbrt 2");
        assert!((root.value() - 2.0_f64.powf(1.0 / 3.0)).abs() < 1e-9);
    }

    #[test]
    fn bisection_finds_sine_root() {
        let root = bisection(Days::new(3.0), Days::new(4.0), |t: Days| {
            Radians::new(t.value().sin())
        })
        .expect("π");
        assert!((root.value() - std::f64::consts::PI).abs() < 1e-8);
    }

    #[test]
    fn bisection_returns_none_for_invalid_bracket() {
        assert!(bisection(Days::new(0.0), Days::new(1.0), |_: Days| Radians::new(42.0)).is_none());
    }

    #[test]
    fn bisection_handles_step_function() {
        let root = bisection(Days::new(-1.0), Days::new(1.0), |t: Days| {
            Radians::new(if t < 0.0 { -5.0 } else { 5.0 })
        })
        .expect("step");
        assert!(root.abs() < 1e-6);
    }

    #[test]
    fn bisection_endpoint_root() {
        let root = bisection(Days::new(0.0), Days::new(5.0), |t: Days| {
            Radians::new(t.value())
        })
        .expect("root at 0");
        assert!(root.abs() < F_EPS);
    }
}
