// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Scan+Brent fallback baseline for threshold crossing discovery.

use crate::event::altitude::search::CROSSING_DEDUPE_EPS;
use crate::event::search::intervals::LabeledCrossing;
use crate::qtty::{Day, Quantity};
use crate::time::{Interval, ModifiedJulianDate};

use super::crossings::{eval_signal, precise_residual_days, SearchDiagnostics, POLY_ZERO_TOL};

type Days = Quantity<Day>;
type Mjd = ModifiedJulianDate;

const MIN_SEGMENT_DAYS: f64 = 1e-6;

/// Scan a window uniformly and refine threshold crossings with Brent's method.
pub(crate) fn scan_labelled_crossings<F>(
    period: Interval<Mjd>,
    step: Days,
    signal: &F,
    threshold_sin: f64,
    time_tolerance: Days,
    max_residual: f64,
    diagnostics: &mut SearchDiagnostics,
) -> Vec<LabeledCrossing>
where
    F: Fn(Mjd) -> f64,
{
    if period.end <= period.start {
        return Vec::new();
    }

    let step_days = step.value().max(MIN_SEGMENT_DAYS);
    let time_tol = time_tolerance.value().max(f64::EPSILON);
    let residual_tol = max_residual.max(POLY_ZERO_TOL);
    let t_start = period.start;
    let t_end = period.end;

    let mut labeled = Vec::new();
    let mut t = t_start;
    let mut prev = eval_signal(signal, t, diagnostics) - threshold_sin;

    while t < t_end {
        let next_t = {
            let proposed = mjd_from_days(mjd_days(t) + step_days);
            if proposed <= t_end {
                proposed
            } else {
                t_end
            }
        };
        let next_v = eval_signal(signal, next_t, diagnostics) - threshold_sin;

        if prev.signum() * next_v.signum() < 0.0 {
            let direction = if prev < 0.0 { 1 } else { -1 };
            if let Some(root_days) = brent_f64(
                mjd_days(t),
                mjd_days(next_t),
                prev,
                next_v,
                |days| precise_residual_days(signal, days, threshold_sin, diagnostics),
                time_tol,
                residual_tol,
            ) {
                let root = mjd_from_days(root_days);
                if root >= t_start && root <= t_end {
                    labeled.push(LabeledCrossing { t: root, direction });
                }
            }
        }

        t = next_t;
        prev = next_v;
    }

    sort_dedup_crossings(&mut labeled);
    labeled
}

pub(crate) fn brent_f64<F>(
    lo: f64,
    hi: f64,
    f_lo: f64,
    f_hi: f64,
    mut f: F,
    time_tol: f64,
    residual_tol: f64,
) -> Option<f64>
where
    F: FnMut(f64) -> f64,
{
    if !lo.is_finite() || !hi.is_finite() || hi < lo {
        return None;
    }
    if f_lo.abs() <= residual_tol {
        return Some(lo);
    }
    if f_hi.abs() <= residual_tol {
        return Some(hi);
    }
    if f_lo.signum() * f_hi.signum() > 0.0 {
        return None;
    }

    let mut a = lo;
    let mut b = hi;
    let mut fa = f_lo;
    let mut fb = f_hi;
    let mut c = a;
    let mut fc = fa;
    let mut d = b - a;
    let mut e = d;

    for _ in 0..100 {
        if fb.signum() * fc.signum() > 0.0 {
            c = a;
            fc = fa;
            d = b - a;
            e = d;
        }
        if fc.abs() < fb.abs() {
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
        }

        let tol1 = 2.0 * f64::EPSILON * b.abs() + time_tol * 0.5;
        let xm = 0.5 * (c - b);
        if xm.abs() <= tol1 {
            if brent_converged(fb, (c - b).abs(), time_tol, residual_tol) {
                return Some(b);
            }
        } else if fb.abs() <= residual_tol {
            return Some(b);
        }

        if e.abs() >= tol1 && fa.abs() > fb.abs() {
            let s = fb / fa;
            let (mut p, mut q) = if (a - c).abs() <= f64::EPSILON {
                (2.0 * xm * s, 1.0 - s)
            } else {
                let q = fa / fc;
                let r = fb / fc;
                (
                    s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0)),
                    (q - 1.0) * (r - 1.0) * (s - 1.0),
                )
            };
            if p > 0.0 {
                q = -q;
            }
            p = p.abs();

            let min1 = 3.0 * xm * q - (tol1 * q).abs();
            let min2 = (e * q).abs();
            if 2.0 * p < min1.min(min2) {
                e = d;
                d = p / q;
            } else {
                d = xm;
                e = d;
            }
        } else {
            d = xm;
            e = d;
        }

        a = b;
        fa = fb;
        if d.abs() > tol1 {
            b += d;
        } else {
            b += tol1.copysign(xm);
        }
        fb = f(b);
        if !fb.is_finite() {
            return None;
        }
    }

    if brent_converged(fb, (c - b).abs(), time_tol, residual_tol) {
        Some(b)
    } else {
        None
    }
}

#[inline]
fn brent_converged(fb: f64, bracket_width: f64, time_tol: f64, residual_tol: f64) -> bool {
    if !fb.is_finite() {
        return false;
    }
    if fb.abs() <= residual_tol {
        return true;
    }
    bracket_width <= time_tol * 4.0
}

#[inline]
fn mjd_days(t: Mjd) -> f64 {
    t.raw().value()
}

#[inline]
fn mjd_from_days(days: f64) -> Mjd {
    ModifiedJulianDate::new(days)
}

fn sort_dedup_crossings(crossings: &mut Vec<LabeledCrossing>) {
    crossings.sort_by(|a, b| a.t.partial_cmp(&b.t).unwrap_or(std::cmp::Ordering::Equal));
    crossings.dedup_by(|a, b| (a.t.raw() - b.t.raw()).abs() < CROSSING_DEDUPE_EPS);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn brent_uses_current_bracket_width_for_convergence() {
        // Flat function: original bracket is wide, current bracket is tiny, residual still bad.
        let root = brent_f64(0.0, 100.0, 1.0, 1.0, |_| 1.0, 1e-12, 1e-12);
        assert!(root.is_none());
    }

    #[test]
    fn brent_finds_linear_root_with_tolerances() {
        let root = brent_f64(0.0, 1.0, -0.5, 0.5, |x| x - 0.25, 1e-12, 1e-12);
        assert!(root.is_some());
        assert!((root.unwrap() - 0.25).abs() < 1e-9);
    }
}
