// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Extrema Finding — Golden‑Section Search & Classification
//!
//! Derivative‑free optimisation for scalar functions of one variable.
//! Finds local minima, maxima, and classifies stationary points.
//!
//! ## Provided routines
//!
//! | Function | Purpose |
//! |----------|---------|
//! | [`minimize`] | Find the *t* that minimises *f(t)* in a bracket |
//! | [`maximize`] | Find the *t* that maximises *f(t)* in a bracket |
//! | [`find_extrema`] | Scan an interval, find all local min/max |
//! | [`classify`] | Determine if a point is a local max, min, or neither |

use super::root_finding;

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// Golden ratio conjugate φ = (√5 − 1) / 2 ≈ 0.6180339887
const PHI: f64 = 0.618_033_988_749_895;

/// Default tolerance for golden‑section convergence.
const DEFAULT_TOL: f64 = 1e-9;

/// Maximum iterations for golden‑section.
const MAX_ITER: usize = 100;

/// Small probe offset for classification (in the caller's time unit).
const PROBE_EPS: f64 = 1e-7;

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
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Extremum {
    pub t: f64,
    pub value: f64,
    pub kind: ExtremumKind,
}

// ---------------------------------------------------------------------------
// Golden‑section search
// ---------------------------------------------------------------------------

/// Find the value of *t* in `[a, b]` that **minimises** `f(t)`.
///
/// Uses golden‑section search (derivative‑free, guaranteed convergence).
/// Returns `(t_min, f(t_min))`.
pub fn minimize<F: Fn(f64) -> f64>(a: f64, b: f64, f: &F) -> (f64, f64) {
    minimize_tol(a, b, f, DEFAULT_TOL)
}

/// Like [`minimize`] but with a caller‑chosen tolerance.
pub fn minimize_tol<F: Fn(f64) -> f64>(mut a: f64, mut b: f64, f: &F, tol: f64) -> (f64, f64) {
    let mut x1 = b - PHI * (b - a);
    let mut x2 = a + PHI * (b - a);
    let mut f1 = f(x1);
    let mut f2 = f(x2);

    for _ in 0..MAX_ITER {
        if (b - a).abs() < tol {
            break;
        }
        if f1 < f2 {
            b = x2;
            x2 = x1;
            f2 = f1;
            x1 = b - PHI * (b - a);
            f1 = f(x1);
        } else {
            a = x1;
            x1 = x2;
            f1 = f2;
            x2 = a + PHI * (b - a);
            f2 = f(x2);
        }
    }

    let t = 0.5 * (a + b);
    (t, f(t))
}

/// Find the value of *t* in `[a, b]` that **maximises** `f(t)`.
///
/// Implemented as `minimize(-f)`.  Returns `(t_max, f(t_max))`.
pub fn maximize<F: Fn(f64) -> f64>(a: f64, b: f64, f: &F) -> (f64, f64) {
    maximize_tol(a, b, f, DEFAULT_TOL)
}

/// Like [`maximize`] but with a caller‑chosen tolerance.
pub fn maximize_tol<F: Fn(f64) -> f64>(a: f64, b: f64, f: &F, tol: f64) -> (f64, f64) {
    let neg_f = |t: f64| -f(t);
    let (t, _neg_val) = minimize_tol(a, b, &neg_f, tol);
    (t, f(t))
}

// ---------------------------------------------------------------------------
// Classification
// ---------------------------------------------------------------------------

/// Classify a candidate extremum by probing f(t ± ε).
///
/// Returns `Some(ExtremumKind)` if the point is a local max or min,
/// `None` if inconclusive (e.g. inflection point or saddle).
pub fn classify<F: Fn(f64) -> f64>(t: f64, f: &F) -> Option<ExtremumKind> {
    let fc = f(t);
    let fl = f(t - PROBE_EPS);
    let fr = f(t + PROBE_EPS);

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

/// Scan `[t_start, t_end]` at `step` intervals, detect monotonicity changes,
/// and refine each extremum with golden‑section.
///
/// Returns a chronologically sorted list of [`Extremum`] values.
pub fn find_extrema<F: Fn(f64) -> f64>(
    t_start: f64,
    t_end: f64,
    step: f64,
    f: &F,
) -> Vec<Extremum> {
    find_extrema_tol(t_start, t_end, step, f, DEFAULT_TOL)
}

/// Like [`find_extrema`] but with a caller‑chosen tolerance.
pub fn find_extrema_tol<F: Fn(f64) -> f64>(
    t_start: f64,
    t_end: f64,
    step: f64,
    f: &F,
    tol: f64,
) -> Vec<Extremum> {
    let mut result = Vec::new();
    if step <= 0.0 || t_start >= t_end {
        return result;
    }

    let mut t0 = t_start;
    let mut f0 = f(t0);
    let mut t1 = (t0 + step).min(t_end);
    let mut f1 = f(t1);
    let mut prev_rising = f1 > f0;

    loop {
        let t2 = (t1 + step).min(t_end);
        if t2 <= t1 {
            break;
        }
        let f2 = f(t2);
        let now_rising = f2 > f1;

        if prev_rising && !now_rising {
            // Was rising, now falling → local maximum in [t0, t2]
            let (t_max, v_max) = maximize_tol(t0, t2, f, tol);
            result.push(Extremum {
                t: t_max,
                value: v_max,
                kind: ExtremumKind::Maximum,
            });
        } else if !prev_rising && now_rising {
            // Was falling, now rising → local minimum in [t0, t2]
            let (t_min, v_min) = minimize_tol(t0, t2, f, tol);
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
pub fn find_extrema_via_derivative<F: Fn(f64) -> f64>(
    t_start: f64,
    t_end: f64,
    step: f64,
    f: &F,
    fd_step: f64,
) -> Vec<Extremum> {
    let deriv = |t: f64| (f(t + fd_step) - f(t - fd_step)) / (2.0 * fd_step);

    let mut result = Vec::new();
    let mut t = t_start;
    let mut prev_d = deriv(t);

    while t < t_end {
        let next_t = (t + step).min(t_end);
        let next_d = deriv(next_t);

        if prev_d * next_d < 0.0 {
            if let Some(root) = root_finding::brent_with_values(t, next_t, prev_d, next_d, &deriv)
            {
                let val = f(root);
                let kind = if let Some(k) = classify(root, f) {
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
                    t: root,
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

    #[test]
    fn minimize_finds_parabola_vertex() {
        let (t, v) = minimize(-5.0, 5.0, &|t: f64| (t - 2.0).powi(2));
        assert!((t - 2.0).abs() < 1e-7, "t = {t}");
        assert!(v < 1e-12, "v = {v}");
    }

    #[test]
    fn maximize_finds_negative_parabola_peak() {
        let (t, v) = maximize(-5.0, 5.0, &|t: f64| -(t - 3.0).powi(2) + 10.0);
        assert!((t - 3.0).abs() < 1e-7, "t = {t}");
        assert!((v - 10.0).abs() < 1e-6, "v = {v}");
    }

    #[test]
    fn classify_detects_maximum() {
        let f = |t: f64| -(t * t);
        assert_eq!(classify(0.0, &f), Some(ExtremumKind::Maximum));
    }

    #[test]
    fn classify_detects_minimum() {
        let f = |t: f64| t * t;
        assert_eq!(classify(0.0, &f), Some(ExtremumKind::Minimum));
    }

    #[test]
    fn classify_returns_none_for_inflection() {
        let f = |t: f64| t * t * t;
        // At t=0 the cubic has an inflection point, not a max/min
        assert_eq!(classify(0.0, &f), None);
    }

    #[test]
    fn find_extrema_sine_wave() {
        // sin(2πt) over [0, 1] has max at t=0.25, min at t=0.75
        let f = |t: f64| (2.0 * std::f64::consts::PI * t).sin();
        let extrema = find_extrema(0.0, 1.0, 0.05, &f);

        assert_eq!(extrema.len(), 2, "expected 2 extrema, got {:?}", extrema);

        let max_ext = extrema.iter().find(|e| e.kind == ExtremumKind::Maximum).unwrap();
        let min_ext = extrema.iter().find(|e| e.kind == ExtremumKind::Minimum).unwrap();

        assert!((max_ext.t - 0.25).abs() < 1e-6, "max at {}", max_ext.t);
        assert!((max_ext.value - 1.0).abs() < 1e-6);
        assert!((min_ext.t - 0.75).abs() < 1e-6, "min at {}", min_ext.t);
        assert!((min_ext.value + 1.0).abs() < 1e-6);
    }

    #[test]
    fn find_extrema_via_derivative_sine() {
        // Shift slightly to avoid derivative sign-change at endpoints
        let f = |t: f64| (2.0 * std::f64::consts::PI * t).sin();
        // Use a step that doesn't align with extrema at t=0.25, 0.75
        let extrema = find_extrema_via_derivative(0.01, 0.99, 0.035, &f, 1e-5);

        assert!(extrema.len() >= 2, "expected ≥2 extrema, got {:?}", extrema);

        let max_ext = extrema.iter().find(|e| e.kind == ExtremumKind::Maximum).unwrap();
        let min_ext = extrema.iter().find(|e| e.kind == ExtremumKind::Minimum).unwrap();

        assert!((max_ext.t - 0.25).abs() < 1e-3, "max at {}", max_ext.t);
        assert!((min_ext.t - 0.75).abs() < 1e-3, "min at {}", min_ext.t);
    }

    #[test]
    fn find_extrema_constant_returns_empty() {
        let extrema = find_extrema(0.0, 10.0, 1.0, &|_| 42.0);
        assert!(extrema.is_empty());
    }

    #[test]
    fn find_extrema_monotone_returns_empty() {
        let extrema = find_extrema(0.0, 10.0, 0.5, &|t| t);
        assert!(extrema.is_empty());
    }

    #[test]
    fn find_extrema_multiple_oscillations() {
        // sin(6πt) over [0,1] → 3 maxima, 2–3 minima
        let f = |t: f64| (6.0 * std::f64::consts::PI * t).sin();
        let extrema = find_extrema(0.0, 1.0, 0.02, &f);

        let n_max = extrema.iter().filter(|e| e.kind == ExtremumKind::Maximum).count();
        let n_min = extrema.iter().filter(|e| e.kind == ExtremumKind::Minimum).count();

        assert!(n_max >= 2, "expected ≥2 maxima, got {n_max}");
        assert!(n_min >= 2, "expected ≥2 minima, got {n_min}");
    }
}
