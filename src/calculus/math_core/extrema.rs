// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Extrema Finding — Golden‑Section Search & Classification
//!
//! Derivative‑free optimisation for scalar functions of one variable.
//! Finds local minima, maxima, and classifies stationary points.
//!
//! All routines operate on `Period<ModifiedJulianDate>` time windows and
//! closures `Fn(ModifiedJulianDate) → Quantity<V>`.
//!
//! ## Provided routines
//!
//! | Function | Purpose |
//! |----------|---------|
//! | [`minimize`] | Find the *t* that minimises *f(t)* in a bracket |
//! | [`maximize`] | Find the *t* that maximises *f(t)* in a bracket |
//! | [`find_extrema`] | Scan a period, find all local min/max |
//! | [`classify`] | Determine if a point is a local max, min, or neither |

use crate::time::{ModifiedJulianDate, Period};
use qtty::*;

use super::root_finding;

type Mjd = ModifiedJulianDate;

#[inline]
fn opposite_sign<V: Unit>(a: Quantity<V>, b: Quantity<V>) -> bool {
    a.signum() * b.signum() < 0.0
}

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// Golden ratio conjugate φ = (√5 − 1) / 2 ≈ 0.6180339887
const PHI: f64 = 0.618_033_988_749_895;

/// Default tolerance for golden‑section convergence.
const DEFAULT_TOL: Days = Days::new(1e-9);

/// Maximum iterations for golden‑section.
const MAX_ITER: usize = 100;

/// Small probe offset for classification (in the caller's time unit).
const PROBE_EPS: Days = Days::new(1e-7);

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

/// The kind of extremum.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExtremumKind {
    Maximum,
    Minimum,
}

/// A located extremum: time, function value, and kind.
///
/// Generic over `V` (value unit).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Extremum<V: Unit> {
    pub t: ModifiedJulianDate,
    pub value: Quantity<V>,
    pub kind: ExtremumKind,
}

// ---------------------------------------------------------------------------
// Golden‑section search
// ---------------------------------------------------------------------------

/// Find the value of *t* in `period` that **minimises** `f(t)`.
///
/// Uses golden‑section search (derivative‑free, guaranteed convergence).
/// Returns `(t_min, f(t_min))`.
pub fn minimize<V, F>(period: Period<Mjd>, f: &F) -> (ModifiedJulianDate, Quantity<V>)
where
    V: Unit,
    F: Fn(ModifiedJulianDate) -> Quantity<V>,
{
    minimize_tol(period, f, DEFAULT_TOL)
}

/// Like [`minimize`] but with a caller‑chosen tolerance.
pub fn minimize_tol<V, F>(
    period: Period<Mjd>,
    f: &F,
    tol: Days,
) -> (ModifiedJulianDate, Quantity<V>)
where
    V: Unit,
    F: Fn(ModifiedJulianDate) -> Quantity<V>,
{
    let mut p = period;
    let mut x1 = p.end - PHI * p.duration();
    let mut x2 = p.start + PHI * p.duration();
    let mut f1: Quantity<V> = f(x1);
    let mut f2: Quantity<V> = f(x2);

    for _ in 0..MAX_ITER {
        if p.duration() < tol {
            break;
        }
        if f1 < f2 {
            p.end = x2;
            x2 = x1;
            f2 = f1;
            x1 = p.end - PHI * p.duration();
            f1 = f(x1);
        } else {
            p.start = x1;
            x1 = x2;
            f1 = f2;
            x2 = p.start + PHI * p.duration();
            f2 = f(x2);
        }
    }

    let t = p.start.mean(p.end);
    (t, f(t))
}

/// Find the value of *t* in `period` that **maximises** `f(t)`.
///
/// Implemented as `minimize(-f)`.  Returns `(t_max, f(t_max))`.
pub fn maximize<V, F>(period: Period<Mjd>, f: &F) -> (ModifiedJulianDate, Quantity<V>)
where
    V: Unit,
    F: Fn(ModifiedJulianDate) -> Quantity<V>,
{
    maximize_tol(period, f, DEFAULT_TOL)
}

/// Like [`maximize`] but with a caller‑chosen tolerance.
pub fn maximize_tol<V, F>(
    period: Period<Mjd>,
    f: &F,
    tol: Days,
) -> (ModifiedJulianDate, Quantity<V>)
where
    V: Unit,
    F: Fn(ModifiedJulianDate) -> Quantity<V>,
{
    let neg_f = |t: Mjd| -f(t);
    let (t, _neg_val) = minimize_tol(period, &neg_f, tol);
    (t, f(t))
}

// ---------------------------------------------------------------------------
// Classification
// ---------------------------------------------------------------------------

/// Classify a candidate extremum by probing f(t ± ε).
///
/// Returns `Some(ExtremumKind)` if the point is a local max or min,
/// `None` if inconclusive (e.g. inflection point or saddle).
pub fn classify<V, F>(t: ModifiedJulianDate, f: &F) -> Option<ExtremumKind>
where
    V: Unit,
    F: Fn(ModifiedJulianDate) -> Quantity<V>,
{
    let tv = t;
    let fc = f(t);
    let fl = f(tv - PROBE_EPS);
    let fr = f(tv + PROBE_EPS);

    if fc >= fl && fc >= fr {
        Some(ExtremumKind::Maximum)
    } else if fc <= fl && fc <= fr {
        Some(ExtremumKind::Minimum)
    } else {
        None
    }
}

// ---------------------------------------------------------------------------
// Scan‑based extrema finder
// ---------------------------------------------------------------------------

/// Scan `period` at `step` intervals, detect monotonicity changes,
/// and refine each extremum with golden‑section.
///
/// Returns a chronologically sorted list of [`Extremum`] values.
pub fn find_extrema<V, F>(period: Period<Mjd>, step: Days, f: &F) -> Vec<Extremum<V>>
where
    V: Unit,
    F: Fn(ModifiedJulianDate) -> Quantity<V>,
{
    find_extrema_tol(period, step, f, DEFAULT_TOL)
}

/// Like [`find_extrema`] but with a caller‑chosen tolerance.
pub fn find_extrema_tol<V, F>(period: Period<Mjd>, step: Days, f: &F, tol: Days) -> Vec<Extremum<V>>
where
    V: Unit,
    F: Fn(ModifiedJulianDate) -> Quantity<V>,
{
    let step_v = step;
    let t_start_v = period.start;
    let t_end_v = period.end;

    let mut result = Vec::new();
    if step_v <= Days::new(0.0) || t_start_v >= t_end_v {
        return result;
    }

    let mut t0 = t_start_v;
    let mut f0 = f(t0);
    let mut t1 = (t0 + step_v).min(t_end_v);
    let mut f1 = f(t1);
    let mut prev_rising = f1 > f0;

    loop {
        let t2 = (t1 + step_v).min(t_end_v);
        if t2 <= t1 {
            break;
        }
        let f2 = f(t2);
        let now_rising = f2 > f1;

        if prev_rising && !now_rising {
            // Was rising, now falling → local maximum in [t0, t2]
            let (t_max, v_max) = maximize_tol(Period::new(t0, t2), f, tol);
            result.push(Extremum {
                t: t_max,
                value: v_max,
                kind: ExtremumKind::Maximum,
            });
        } else if !prev_rising && now_rising {
            // Was falling, now rising → local minimum in [t0, t2]
            let (t_min, v_min) = minimize_tol(Period::new(t0, t2), f, tol);
            result.push(Extremum {
                t: t_min,
                value: v_min,
                kind: ExtremumKind::Minimum,
            });
        }

        t0 = t1;
        f0 = f1;
        t1 = t2;
        f1 = f2;
        prev_rising = now_rising;
    }

    // Suppress the unused-variable warning for f0
    let _ = f0;

    result.sort_by(|a, b| a.t.partial_cmp(&b.t).unwrap());
    result
}

/// Find roots of the numerical derivative f'(t) ≈ 0 using finite differences
/// and Brent's method.  This is an alternative extrema finder that scans for
/// sign changes in the finite‑difference derivative.
///
/// The derivative is estimated as: f'(t) ≈ (f(t+h) − f(t−h)) / (2h)
pub fn find_extrema_via_derivative<V, F>(
    period: Period<Mjd>,
    step: Days,
    f: &F,
    fd_step: Days,
) -> Vec<Extremum<V>>
where
    V: Unit,
    F: Fn(ModifiedJulianDate) -> Quantity<V>,
{
    let step_v = step;
    let t_start_v = period.start;
    let t_end_v = period.end;
    let fd_v = fd_step;

    let deriv = |t: Mjd| -> Quantity<V> {
        let tv = t;
        let fwd = f(tv + fd_v);
        let bwd = f(tv - fd_v);
        (fwd - bwd) / (fd_v + fd_v).value()
    };

    let mut result = Vec::new();
    let mut t = t_start_v;
    let mut prev_d = deriv(t);

    while t < t_end_v {
        let next_t = (t + step_v).min(t_end_v);
        let next_d = deriv(next_t);
        if opposite_sign(prev_d, next_d) {
            if let Some(root_mjd) =
                root_finding::brent_with_values(Period::new(t, next_t), prev_d, next_d, deriv)
            {
                let val = f(root_mjd);
                let kind = if let Some(k) = classify::<V, _>(root_mjd, f) {
                    k
                } else {
                    // Default based on derivative sign change direction
                    if prev_d > 0.0 {
                        ExtremumKind::Maximum
                    } else {
                        ExtremumKind::Minimum
                    }
                };
                result.push(Extremum {
                    t: root_mjd,
                    value: val,
                    kind,
                });
            }
        }

        t = next_t;
        prev_d = next_d;
    }

    result.sort_by(|a, b| a.t.partial_cmp(&b.t).unwrap());
    result
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use qtty::Radian;

    type Radians = Quantity<Radian>;

    fn mjd(v: f64) -> Mjd {
        Mjd::new(v)
    }
    fn period(a: f64, b: f64) -> Period<Mjd> {
        Period::new(mjd(a), mjd(b))
    }

    #[test]
    fn minimize_finds_parabola_vertex() {
        let (t, v) = minimize(period(-5.0, 5.0), &|t: Mjd| {
            Radians::new((t.value() - 2.0).powi(2))
        });
        assert!((t.value() - 2.0).abs() < 1e-7, "t = {}", t.value());
        assert!(v < 1e-12, "v = {}", v.value());
    }

    #[test]
    fn maximize_finds_negative_parabola_peak() {
        let (t, v) = maximize(period(-5.0, 5.0), &|t: Mjd| {
            Radians::new(-(t.value() - 3.0).powi(2) + 10.0)
        });
        assert!((t.value() - 3.0).abs() < 1e-7, "t = {}", t.value());
        assert!((v.value() - 10.0).abs() < 1e-6, "v = {}", v.value());
    }

    #[test]
    fn classify_detects_maximum() {
        let f = |t: Mjd| Radians::new(-(t.value() * t.value()));
        assert_eq!(
            classify::<Radian, _>(mjd(0.0), &f),
            Some(ExtremumKind::Maximum)
        );
    }

    #[test]
    fn classify_detects_minimum() {
        let f = |t: Mjd| Radians::new(t.value() * t.value());
        assert_eq!(
            classify::<Radian, _>(mjd(0.0), &f),
            Some(ExtremumKind::Minimum)
        );
    }

    #[test]
    fn classify_returns_none_for_inflection() {
        let f = |t: Mjd| Radians::new(t.value() * t.value() * t.value());
        // At t=0 the cubic has an inflection point, not a max/min
        assert_eq!(classify::<Radian, _>(mjd(0.0), &f), None);
    }

    #[test]
    fn find_extrema_sine_wave() {
        // sin(2πt) over [0, 1] has max at t=0.25, min at t=0.75
        let f = |t: Mjd| Radians::new((2.0 * std::f64::consts::PI * t.value()).sin());
        let extrema: Vec<Extremum<Radian>> = find_extrema(period(0.0, 1.0), Days::new(0.05), &f);

        assert_eq!(extrema.len(), 2, "expected 2 extrema, got {:?}", extrema);

        let max_ext = extrema
            .iter()
            .find(|e| e.kind == ExtremumKind::Maximum)
            .unwrap();
        let min_ext = extrema
            .iter()
            .find(|e| e.kind == ExtremumKind::Minimum)
            .unwrap();

        assert!(
            (max_ext.t.value() - 0.25).abs() < 1e-6,
            "max at {}",
            max_ext.t.value()
        );
        assert!((max_ext.value.value() - 1.0).abs() < 1e-6);
        assert!(
            (min_ext.t.value() - 0.75).abs() < 1e-6,
            "min at {}",
            min_ext.t.value()
        );
        assert!((min_ext.value.value() + 1.0).abs() < 1e-6);
    }

    #[test]
    fn find_extrema_via_derivative_sine() {
        // Shift slightly to avoid derivative sign-change at endpoints
        let f = |t: Mjd| Radians::new((2.0 * std::f64::consts::PI * t.value()).sin());
        // Use a step that doesn't align with extrema at t=0.25, 0.75
        let extrema: Vec<Extremum<Radian>> =
            find_extrema_via_derivative(period(0.01, 0.99), Days::new(0.035), &f, Days::new(1e-5));

        assert!(extrema.len() >= 2, "expected ≥2 extrema, got {:?}", extrema);

        let max_ext = extrema
            .iter()
            .find(|e| e.kind == ExtremumKind::Maximum)
            .unwrap();
        let min_ext = extrema
            .iter()
            .find(|e| e.kind == ExtremumKind::Minimum)
            .unwrap();

        assert!(
            (max_ext.t.value() - 0.25).abs() < 1e-3,
            "max at {}",
            max_ext.t.value()
        );
        assert!(
            (min_ext.t.value() - 0.75).abs() < 1e-3,
            "min at {}",
            min_ext.t.value()
        );
    }

    #[test]
    fn find_extrema_constant_returns_empty() {
        let extrema: Vec<Extremum<Radian>> =
            find_extrema(period(0.0, 10.0), Days::new(1.0), &|_: Mjd| {
                Radians::new(42.0)
            });
        assert!(extrema.is_empty());
    }

    #[test]
    fn find_extrema_monotone_returns_empty() {
        let extrema: Vec<Extremum<Radian>> =
            find_extrema(period(0.0, 10.0), Days::new(0.5), &|t: Mjd| {
                Radians::new(t.value())
            });
        assert!(extrema.is_empty());
    }

    #[test]
    fn find_extrema_multiple_oscillations() {
        // sin(6πt) over [0,1] → 3 maxima, 2–3 minima
        let f = |t: Mjd| Radians::new((6.0 * std::f64::consts::PI * t.value()).sin());
        let extrema: Vec<Extremum<Radian>> = find_extrema(period(0.0, 1.0), Days::new(0.02), &f);

        let n_max = extrema
            .iter()
            .filter(|e| e.kind == ExtremumKind::Maximum)
            .count();
        let n_min = extrema
            .iter()
            .filter(|e| e.kind == ExtremumKind::Minimum)
            .count();

        assert!(n_max >= 2, "expected ≥2 maxima, got {n_max}");
        assert!(n_min >= 2, "expected ≥2 minima, got {n_min}");
    }
}
