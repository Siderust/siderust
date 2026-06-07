// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Chebyshev-based crossing discovery for smooth altitude signals.
//!
//! This module fits `sin(altitude(t)) - sin(threshold)` on short time
//! segments, solves roots of the fitted polynomial on `[-1, 1]`, validates the
//! candidates against the precise signal, and falls back per segment to the
//! legacy scan+Brent path when the polynomial is not trustworthy.

use crate::event::altitude::search::{CrossingAlgorithm, SearchOptsV2};
use crate::event::search::intervals::LabeledCrossing;
use crate::qtty::{Day, Quantity};
use crate::time::{Interval, ModifiedJulianDate};

type Days = Quantity<Day>;
type Mjd = ModifiedJulianDate;

const UNIT_ROOT_TOL: f64 = 1e-13;
const POLY_ZERO_TOL: f64 = 1e-12;
const DEDUPE_X_EPS: f64 = 1e-10;
const DEDUPE_T_EPS: Days = Days::new(1e-8);
const MIN_SEGMENT_DAYS: f64 = 1e-6;

/// Runtime counters for Chebyshev crossing searches.
#[derive(Debug, Default, Clone, Copy)]
pub(crate) struct SearchDiagnostics {
    /// Chebyshev segments attempted.
    pub(crate) segments: usize,
    /// Candidate polynomial roots produced before validation.
    pub(crate) polynomial_roots: usize,
    /// Segments that used scan+Brent fallback.
    pub(crate) fallback_segments: usize,
    /// Precise signal evaluations performed by this search layer.
    pub(crate) precise_evaluations: usize,
}

#[derive(Debug, Clone, Copy)]
struct Segment {
    start: Mjd,
    end: Mjd,
}

impl Segment {
    fn from_interval(period: Interval<Mjd>) -> Self {
        Self {
            start: period.start,
            end: period.end,
        }
    }

    fn as_interval(self) -> Interval<Mjd> {
        Interval::new(self.start, self.end)
    }

    fn span_days(self) -> f64 {
        (self.end.raw() - self.start.raw()).value()
    }

    fn contains_positive_time(self) -> bool {
        self.span_days() > 0.0
    }

    fn split(self) -> (Self, Self) {
        let mid = mjd_from_days(mjd_days(self.start) + self.span_days() * 0.5);
        (
            Self {
                start: self.start,
                end: mid,
            },
            Self {
                start: mid,
                end: self.end,
            },
        )
    }

    fn time_from_unit(self, x: f64) -> Mjd {
        let u = ((x.clamp(-1.0, 1.0)) + 1.0) * 0.5;
        mjd_from_days(mjd_days(self.start) + self.span_days() * u)
    }
}

/// Chebyshev polynomial `sum c_k T_k(x)` on `[-1, 1]`.
#[derive(Debug, Clone)]
pub(crate) struct ChebPoly {
    coeffs: Vec<f64>,
}

impl ChebPoly {
    /// Fit a Chebyshev polynomial of `degree` on `[-1, 1]`.
    pub(crate) fn fit<F>(degree: usize, mut f: F) -> Self
    where
        F: FnMut(f64) -> f64,
    {
        let n = degree.saturating_add(1).max(2);
        let mut values = Vec::with_capacity(n);
        for j in 0..n {
            let theta = std::f64::consts::PI * (j as f64 + 0.5) / n as f64;
            values.push(f(theta.cos()));
        }

        let mut coeffs = vec![0.0; n];
        for (k, coeff) in coeffs.iter_mut().enumerate() {
            let mut sum = 0.0;
            for (j, value) in values.iter().enumerate() {
                let theta = std::f64::consts::PI * (j as f64 + 0.5) / n as f64;
                sum += value * (k as f64 * theta).cos();
            }
            *coeff = 2.0 * sum / n as f64;
        }
        coeffs[0] *= 0.5;

        Self { coeffs }
    }

    /// Evaluate the polynomial using Clenshaw recurrence.
    #[inline]
    pub(crate) fn eval(&self, x: f64) -> f64 {
        if self.coeffs.is_empty() {
            return 0.0;
        }

        let mut b1 = 0.0;
        let mut b2 = 0.0;
        for &coeff in self.coeffs.iter().skip(1).rev() {
            let b0 = 2.0 * x * b1 - b2 + coeff;
            b2 = b1;
            b1 = b0;
        }
        x * b1 - b2 + self.coeffs[0]
    }

    /// Return the derivative polynomial.
    pub(crate) fn derivative(&self) -> Self {
        let n = self.coeffs.len().saturating_sub(1);
        if n == 0 {
            return Self { coeffs: vec![0.0] };
        }
        if n == 1 {
            return Self {
                coeffs: vec![self.coeffs[1]],
            };
        }

        let mut deriv = vec![0.0; n];
        deriv[n - 1] = 2.0 * n as f64 * self.coeffs[n];
        if n >= 2 {
            for k in (1..=(n - 2)).rev() {
                let next = if k + 2 < deriv.len() {
                    deriv[k + 2]
                } else {
                    0.0
                };
                deriv[k] = next + 2.0 * (k as f64 + 1.0) * self.coeffs[k + 1];
            }
            deriv[0] = deriv.get(2).copied().unwrap_or(0.0) * 0.5 + self.coeffs[1];
        }

        Self { coeffs: deriv }
    }

    /// Sum of the last Chebyshev coefficients, used as a cheap fit-quality
    /// signal.
    pub(crate) fn tail_norm(&self) -> f64 {
        let tail = self.coeffs.len().min(4);
        self.coeffs.iter().rev().take(tail).map(|c| c.abs()).sum()
    }

    /// Roots in `[-1, 1]`, segmented by derivative roots and refined on the
    /// polynomial.
    pub(crate) fn roots(&self) -> Vec<f64> {
        let mut roots = self.roots_recursive(0);
        sort_dedup_f64(&mut roots, DEDUPE_X_EPS);
        roots
    }

    fn roots_recursive(&self, depth: usize) -> Vec<f64> {
        if self.coeffs.len() <= 1 {
            return Vec::new();
        }
        if self.coeffs.len() == 2 {
            let a = self.coeffs[1];
            if a.abs() <= POLY_ZERO_TOL {
                return Vec::new();
            }
            let x = -self.coeffs[0] / a;
            if (-1.0 - UNIT_ROOT_TOL..=1.0 + UNIT_ROOT_TOL).contains(&x) {
                return vec![x.clamp(-1.0, 1.0)];
            }
            return Vec::new();
        }

        let mut points = vec![-1.0, 1.0];
        if depth < self.coeffs.len() {
            points.extend(self.derivative().roots_recursive(depth + 1));
        }
        sort_dedup_f64(&mut points, DEDUPE_X_EPS);

        let mut roots = Vec::new();
        for &x in &points {
            if self.eval(x).abs() <= POLY_ZERO_TOL {
                roots.push(x);
            }
        }

        for pair in points.windows(2) {
            let a = pair[0];
            let b = pair[1];
            if b - a <= UNIT_ROOT_TOL {
                continue;
            }
            let fa = self.eval(a);
            let fb = self.eval(b);
            if fa.abs() <= POLY_ZERO_TOL || fb.abs() <= POLY_ZERO_TOL {
                continue;
            }
            if fa.signum() * fb.signum() < 0.0 {
                if let Some(root) = brent_f64(a, b, fa, fb, |x| self.eval(x), UNIT_ROOT_TOL) {
                    roots.push(root.clamp(-1.0, 1.0));
                }
            }
        }

        roots
    }
}

/// Find labelled crossings of a precise `sin_altitude` signal.
pub(crate) fn find_directed_crossings<F>(
    period: Interval<Mjd>,
    fallback_step: Days,
    signal: &F,
    threshold_sin: f64,
    opts: SearchOptsV2,
) -> (Vec<LabeledCrossing>, bool, SearchDiagnostics)
where
    F: Fn(Mjd) -> f64,
{
    let mut diagnostics = SearchDiagnostics::default();
    if period.end <= period.start {
        return (Vec::new(), false, diagnostics);
    }

    let start_above = eval_signal(signal, period.start, &mut diagnostics) > threshold_sin;
    let fallback_step = opts.scan_step_days.unwrap_or(fallback_step);

    let segment_days = opts.chebyshev.segment_length.value();
    let window_days = (period.end.raw() - period.start.raw()).value();
    let short_window = window_days <= segment_days.max(fallback_step.value()) * 0.5;
    let debug_long_window =
        cfg!(debug_assertions) && opts.algorithm == CrossingAlgorithm::Auto && window_days > 60.0;
    let use_scan = opts.algorithm == CrossingAlgorithm::ScanBrent
        || opts.scan_step_days.is_some()
        || (opts.algorithm == CrossingAlgorithm::Auto && short_window)
        || debug_long_window;

    if use_scan {
        diagnostics.fallback_segments = 1;
        let crossings = scan_directed_crossings(
            period,
            fallback_step,
            signal,
            threshold_sin,
            &mut diagnostics,
        );
        return (crossings, start_above, diagnostics);
    }

    let mut crossings = Vec::new();
    let mut t = period.start;
    while t < period.end {
        let next = {
            let proposed = mjd_from_days((t.raw() + opts.chebyshev.segment_length).value());
            if proposed <= period.end {
                proposed
            } else {
                period.end
            }
        };
        let segment = Segment {
            start: t,
            end: next,
        };
        let mut search = SegmentSearch {
            period,
            fallback_step,
            signal,
            threshold_sin,
            opts,
            out: &mut crossings,
            diagnostics: &mut diagnostics,
        };
        append_segment_crossings(&mut search, segment, 0);
        t = next;
    }

    sort_dedup_crossings(&mut crossings);
    (crossings, start_above, diagnostics)
}

struct SegmentSearch<'a, F>
where
    F: Fn(Mjd) -> f64,
{
    period: Interval<Mjd>,
    fallback_step: Days,
    signal: &'a F,
    threshold_sin: f64,
    opts: SearchOptsV2,
    out: &'a mut Vec<LabeledCrossing>,
    diagnostics: &'a mut SearchDiagnostics,
}

fn append_segment_crossings<F>(search: &mut SegmentSearch<'_, F>, segment: Segment, depth: usize)
where
    F: Fn(Mjd) -> f64,
{
    if !segment.contains_positive_time() {
        return;
    }

    search.diagnostics.segments += 1;
    if segment.span_days() < MIN_SEGMENT_DAYS {
        search.diagnostics.fallback_segments += 1;
        search.out.extend(scan_directed_crossings(
            segment.as_interval(),
            search.fallback_step,
            search.signal,
            search.threshold_sin,
            search.diagnostics,
        ));
        return;
    }

    let degree = search.opts.chebyshev.degree.clamp(4, 64);
    let poly = ChebPoly::fit(degree, |x| {
        eval_signal(search.signal, segment.time_from_unit(x), search.diagnostics)
            - search.threshold_sin
    });

    let tail_norm = poly.tail_norm();
    if !tail_norm.is_finite() || tail_norm > search.opts.chebyshev.max_tail_norm {
        if search.opts.chebyshev.adaptive_split && depth < search.opts.chebyshev.max_split_depth {
            let (left, right) = segment.split();
            append_segment_crossings(search, left, depth + 1);
            append_segment_crossings(search, right, depth + 1);
            return;
        }

        search.diagnostics.fallback_segments += 1;
        search.out.extend(scan_directed_crossings(
            segment.as_interval(),
            search.fallback_step,
            search.signal,
            search.threshold_sin,
            search.diagnostics,
        ));
        return;
    }

    let derivative = poly.derivative();
    let roots = poly.roots();
    search.diagnostics.polynomial_roots += roots.len();

    let mut segment_crossings = Vec::new();
    for root_x in roots {
        let slope = derivative.eval(root_x).abs();
        if slope < search.opts.chebyshev.min_slope {
            search.diagnostics.fallback_segments += 1;
            search.out.extend(scan_directed_crossings(
                segment.as_interval(),
                search.fallback_step,
                search.signal,
                search.threshold_sin,
                search.diagnostics,
            ));
            return;
        }

        let Some(crossing) = validate_candidate(
            search.period,
            segment,
            root_x,
            search.signal,
            search.threshold_sin,
            search.opts,
            search.diagnostics,
        ) else {
            search.diagnostics.fallback_segments += 1;
            search.out.extend(scan_directed_crossings(
                segment.as_interval(),
                search.fallback_step,
                search.signal,
                search.threshold_sin,
                search.diagnostics,
            ));
            return;
        };
        segment_crossings.push(crossing);
    }

    search.out.extend(segment_crossings);
}

fn validate_candidate<F>(
    period: Interval<Mjd>,
    segment: Segment,
    root_x: f64,
    signal: &F,
    threshold_sin: f64,
    opts: SearchOptsV2,
    diagnostics: &mut SearchDiagnostics,
) -> Option<LabeledCrossing>
where
    F: Fn(Mjd) -> f64,
{
    let initial = segment.time_from_unit(root_x);
    let refined = if opts.chebyshev.refine {
        refine_precise_root(segment, initial, signal, threshold_sin, opts, diagnostics)?
    } else {
        initial
    };

    let residual = (eval_signal(signal, refined, diagnostics) - threshold_sin).abs();
    if !residual.is_finite() || residual > opts.chebyshev.max_residual {
        return None;
    }

    let direction = classify_direction(period, refined, signal, threshold_sin, opts, diagnostics)?;
    Some(LabeledCrossing {
        t: refined,
        direction,
    })
}

fn refine_precise_root<F>(
    segment: Segment,
    initial: Mjd,
    signal: &F,
    threshold_sin: f64,
    opts: SearchOptsV2,
    diagnostics: &mut SearchDiagnostics,
) -> Option<Mjd>
where
    F: Fn(Mjd) -> f64,
{
    let segment_start = mjd_days(segment.start);
    let segment_end = mjd_days(segment.end);
    let margin = opts.chebyshev.refine_margin.value().max(1e-5);
    let mut t = mjd_days(initial).clamp(segment_start, segment_end);
    let mut best_t = t;
    let mut best_residual = precise_residual_days(signal, t, threshold_sin, diagnostics).abs();

    for _ in 0..6 {
        if best_residual <= opts.chebyshev.max_residual {
            return Some(mjd_from_days(best_t));
        }

        let h = margin.min((segment_end - segment_start) * 0.05).max(1e-7);
        let lo = (t - h).max(segment_start);
        let hi = (t + h).min(segment_end);
        if hi <= lo {
            break;
        }

        let f_lo = precise_residual_days(signal, lo, threshold_sin, diagnostics);
        let f_hi = precise_residual_days(signal, hi, threshold_sin, diagnostics);
        let slope = (f_hi - f_lo) / (hi - lo);
        if !slope.is_finite() || slope.abs() < 1e-12 {
            break;
        }

        let f_t = precise_residual_days(signal, t, threshold_sin, diagnostics);
        let next = (t - f_t / slope).clamp(segment_start, segment_end);
        if (next - t).abs() < opts.time_tolerance.value().max(1e-12)
            || (next - mjd_days(initial)).abs() <= margin * 4.0
        {
            t = next;
        } else {
            break;
        }

        let residual = precise_residual_days(signal, t, threshold_sin, diagnostics).abs();
        if residual < best_residual {
            best_residual = residual;
            best_t = t;
        }
    }

    if best_residual <= opts.chebyshev.max_residual {
        return Some(mjd_from_days(best_t));
    }

    let mut radius = margin;
    for _ in 0..8 {
        let lo = (mjd_days(initial) - radius).max(segment_start);
        let hi = (mjd_days(initial) + radius).min(segment_end);
        if hi > lo {
            let f_lo = precise_residual_days(signal, lo, threshold_sin, diagnostics);
            let f_hi = precise_residual_days(signal, hi, threshold_sin, diagnostics);
            if f_lo.abs() <= opts.chebyshev.max_residual {
                return Some(mjd_from_days(lo));
            }
            if f_hi.abs() <= opts.chebyshev.max_residual {
                return Some(mjd_from_days(hi));
            }
            if f_lo.signum() * f_hi.signum() < 0.0 {
                let tol = opts.time_tolerance.value().max(1e-12);
                return brent_f64(
                    lo,
                    hi,
                    f_lo,
                    f_hi,
                    |days| precise_residual_days(signal, days, threshold_sin, diagnostics),
                    tol,
                )
                .map(mjd_from_days);
            }
        }
        radius *= 2.0;
        if looser_radius_exhausts_segment(mjd_days(initial), radius, segment_start, segment_end) {
            break;
        }
    }

    None
}

fn classify_direction<F>(
    period: Interval<Mjd>,
    t: Mjd,
    signal: &F,
    threshold_sin: f64,
    opts: SearchOptsV2,
    diagnostics: &mut SearchDiagnostics,
) -> Option<i32>
where
    F: Fn(Mjd) -> f64,
{
    let margin = opts.chebyshev.refine_margin.value();
    let probe = margin.clamp(1e-6, 1e-3);
    let t_days = mjd_days(t);
    let start = mjd_days(period.start);
    let end = mjd_days(period.end);

    let left_days = (t_days - probe).max(start);
    let right_days = (t_days + probe).min(end);
    if right_days <= left_days {
        return None;
    }

    let left = precise_residual_days(signal, left_days, threshold_sin, diagnostics);
    let right = precise_residual_days(signal, right_days, threshold_sin, diagnostics);
    let eps = opts.chebyshev.max_residual.max(1e-12);

    if left <= eps && right > eps {
        Some(1)
    } else if left > eps && right <= eps {
        Some(-1)
    } else if left < -eps && right >= -eps {
        Some(1)
    } else if left >= -eps && right < -eps {
        Some(-1)
    } else {
        None
    }
}

fn scan_directed_crossings<F>(
    period: Interval<Mjd>,
    step: Days,
    signal: &F,
    threshold_sin: f64,
    diagnostics: &mut SearchDiagnostics,
) -> Vec<LabeledCrossing>
where
    F: Fn(Mjd) -> f64,
{
    if period.end <= period.start {
        return Vec::new();
    }

    let step_days = step.value().max(MIN_SEGMENT_DAYS);
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
            let tol = 1e-9;
            if let Some(root_days) = brent_f64(
                mjd_days(t),
                mjd_days(next_t),
                prev,
                next_v,
                |days| precise_residual_days(signal, days, threshold_sin, diagnostics),
                tol,
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

fn brent_f64<F>(lo: f64, hi: f64, f_lo: f64, f_hi: f64, mut f: F, tol: f64) -> Option<f64>
where
    F: FnMut(f64) -> f64,
{
    if !lo.is_finite() || !hi.is_finite() || hi < lo {
        return None;
    }
    if f_lo.abs() <= POLY_ZERO_TOL {
        return Some(lo);
    }
    if f_hi.abs() <= POLY_ZERO_TOL {
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

        let tol1 = 2.0 * f64::EPSILON * b.abs() + tol * 0.5;
        let xm = 0.5 * (c - b);
        if xm.abs() <= tol1 || fb.abs() <= POLY_ZERO_TOL {
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

    Some(b)
}

fn eval_signal<F>(signal: &F, t: Mjd, diagnostics: &mut SearchDiagnostics) -> f64
where
    F: Fn(Mjd) -> f64,
{
    diagnostics.precise_evaluations += 1;
    signal(t)
}

fn precise_residual_days<F>(
    signal: &F,
    days: f64,
    threshold_sin: f64,
    diagnostics: &mut SearchDiagnostics,
) -> f64
where
    F: Fn(Mjd) -> f64,
{
    eval_signal(signal, mjd_from_days(days), diagnostics) - threshold_sin
}

#[inline]
fn mjd_days(t: Mjd) -> f64 {
    t.raw().value()
}

#[inline]
fn mjd_from_days(days: f64) -> Mjd {
    ModifiedJulianDate::new(days)
}

fn looser_radius_exhausts_segment(center: f64, radius: f64, start: f64, end: f64) -> bool {
    center - radius <= start && center + radius >= end
}

fn sort_dedup_f64(values: &mut Vec<f64>, eps: f64) {
    values.retain(|v| v.is_finite());
    values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    values.dedup_by(|a, b| (*a - *b).abs() <= eps);
}

fn sort_dedup_crossings(crossings: &mut Vec<LabeledCrossing>) {
    crossings.sort_by(|a, b| a.t.partial_cmp(&b.t).unwrap_or(std::cmp::Ordering::Equal));
    crossings.dedup_by(|a, b| (a.t.raw() - b.t.raw()).abs() < DEDUPE_T_EPS);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::event::altitude::search::ChebyshevOptions;

    fn period(a: f64, b: f64) -> Interval<Mjd> {
        Interval::new(ModifiedJulianDate::new(a), ModifiedJulianDate::new(b))
    }

    #[test]
    fn cheb_poly_evaluates_fit() {
        let poly = ChebPoly::fit(12, |x| x * x - 0.25);
        assert!(poly.eval(0.5).abs() < 1e-12);
        assert!((poly.eval(0.0) + 0.25).abs() < 1e-12);
    }

    #[test]
    fn cheb_poly_derivative_is_consistent() {
        let poly = ChebPoly::fit(8, |x| 2.0 * x * x - 1.0);
        let deriv = poly.derivative();
        assert!((deriv.eval(0.25) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn polynomial_roots_find_crossings() {
        let poly = ChebPoly::fit(12, |x| x * x - 0.25);
        let roots = poly.roots();
        assert_eq!(roots.len(), 2, "{roots:?}");
        assert!((roots[0] + 0.5).abs() < 1e-10, "{roots:?}");
        assert!((roots[1] - 0.5).abs() < 1e-10, "{roots:?}");
    }

    #[test]
    fn directed_crossings_match_sine_wave() {
        let opts = SearchOptsV2 {
            chebyshev: ChebyshevOptions {
                segment_length: Days::new(0.5),
                degree: 14,
                max_tail_norm: 1e-10,
                max_residual: 1e-10,
                ..ChebyshevOptions::default()
            },
            ..SearchOptsV2::default()
        };
        let signal = |t: Mjd| (2.0 * std::f64::consts::PI * (t.raw().value() + 0.05)).sin();
        let (crossings, start_above, diagnostics) =
            find_directed_crossings(period(0.0, 1.0), Days::new(0.02), &signal, 0.0, opts);

        assert!(start_above);
        assert_eq!(crossings.len(), 2, "{crossings:?} {diagnostics:?}");
        assert_eq!(crossings[0].direction, -1);
        assert_eq!(crossings[1].direction, 1);
        assert!((crossings[0].t.raw().value() - 0.45).abs() < 1e-8);
        assert!((crossings[1].t.raw().value() - 0.95).abs() < 1e-8);
        assert_eq!(diagnostics.fallback_segments, 0);
    }

    #[test]
    fn explicit_scan_uses_fallback_path() {
        let opts = SearchOptsV2 {
            algorithm: CrossingAlgorithm::ScanBrent,
            ..SearchOptsV2::default()
        };
        let signal = |t: Mjd| t.raw().value() - 0.5;
        let (crossings, _, diagnostics) =
            find_directed_crossings(period(0.0, 1.0), Days::new(0.1), &signal, 0.0, opts);
        assert_eq!(crossings.len(), 1);
        assert_eq!(diagnostics.fallback_segments, 1);
    }
}
