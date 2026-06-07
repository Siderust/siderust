// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Chebyshev-first labelled crossing discovery for smooth altitude signals.
//!
//! Fits `sin(altitude(t)) - sin(threshold)` on short time segments, solves roots
//! of the fitted polynomial on `[-1, 1]`, validates candidates against the
//! precise signal, and falls back per segment to scan+Brent when the polynomial
//! is not trustworthy.

use cheby::fit_dyn_from_fn;

use crate::event::altitude::search::{CrossingAlgorithm, SearchOptsV2};
use crate::event::search::intervals::LabeledCrossing;
use crate::event::search::scan_fallback;
use crate::qtty::{Day, Quantity};
use crate::time::{Interval, ModifiedJulianDate};

type Days = Quantity<Day>;
type Mjd = ModifiedJulianDate;

pub(crate) const POLY_ZERO_TOL: f64 = 1e-12;
const DEDUPE_T_EPS: Days = Days::new(1e-8);
const MIN_SEGMENT_DAYS: f64 = 1e-6;
const TAIL_NORM_COEFFS: usize = 4;

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

/// Internal primitive: find labelled crossings of a precise `sin_altitude` signal.
pub(crate) fn find_labelled_crossings<F>(
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

    if opts.algorithm == CrossingAlgorithm::ScanBrent || opts.scan_step_days.is_some() {
        diagnostics.fallback_segments = 1;
        let crossings = scan_fallback::scan_labelled_crossings(
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
        fallback_segment(search, segment);
        return;
    }

    let degree = search.opts.chebyshev.degree.clamp(4, 64);
    let Ok(poly) = fit_dyn_from_fn(degree, |x| {
        eval_signal(search.signal, segment.time_from_unit(x), search.diagnostics)
            - search.threshold_sin
    }) else {
        fallback_segment(search, segment);
        return;
    };

    let tail_norm = poly.tail_norm(TAIL_NORM_COEFFS);
    if !tail_norm.is_finite() || tail_norm > search.opts.chebyshev.max_tail_norm {
        if search.opts.chebyshev.adaptive_split && depth < search.opts.chebyshev.max_split_depth {
            let (left, right) = segment.split();
            append_segment_crossings(search, left, depth + 1);
            append_segment_crossings(search, right, depth + 1);
            return;
        }
        fallback_segment(search, segment);
        return;
    }

    let derivative = poly.derivative();
    let roots = poly.roots();
    search.diagnostics.polynomial_roots += roots.len();

    let mut segment_crossings = Vec::new();
    for root_x in roots {
        let slope = derivative.evaluate(root_x).abs();
        if slope < search.opts.chebyshev.min_slope {
            fallback_segment(search, segment);
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
            fallback_segment(search, segment);
            return;
        };
        segment_crossings.push(crossing);
    }

    search.out.extend(segment_crossings);
}

fn fallback_segment<F>(search: &mut SegmentSearch<'_, F>, segment: Segment)
where
    F: Fn(Mjd) -> f64,
{
    search.diagnostics.fallback_segments += 1;
    search.out.extend(scan_fallback::scan_labelled_crossings(
        segment.as_interval(),
        search.fallback_step,
        search.signal,
        search.threshold_sin,
        search.diagnostics,
    ));
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
                return scan_fallback::brent_f64(
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

pub(crate) fn eval_signal<F>(signal: &F, t: Mjd, diagnostics: &mut SearchDiagnostics) -> f64
where
    F: Fn(Mjd) -> f64,
{
    diagnostics.precise_evaluations += 1;
    signal(t)
}

pub(crate) fn precise_residual_days<F>(
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
    fn labelled_crossings_match_sine_wave() {
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
            find_labelled_crossings(period(0.0, 1.0), Days::new(0.02), &signal, 0.0, opts);

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
            find_labelled_crossings(period(0.0, 1.0), Days::new(0.1), &signal, 0.0, opts);
        assert_eq!(crossings.len(), 1);
        assert_eq!(diagnostics.fallback_segments, 1);
    }
}
