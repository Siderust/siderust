// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Interval Assembly — "In‑Range" Period Detection
//!
//! Assembles time periods where a scalar function satisfies
//! `h_min ≤ f(t) ≤ h_max`, handling:
//!
//! * Multiple crossings per window (fast movers, satellites)
//! * Tangencies (touching a threshold without crossing)
//! * Always‑above / always‑below cases
//! * Single‑threshold (above only) as a special case
//!
//! All routines operate on `Period<MJD>` time windows and
//! closures `Fn(ModifiedJulianDate) → Quantity<V>`.

use crate::time::{complement_within, intersect_periods, ModifiedJulianDate, Period, MJD};
use qtty::{Day, Quantity, Unit};

use super::root_finding;

type Mjd = ModifiedJulianDate;
type Days = Quantity<Day>;

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// Tiny probe offset for classifying crossing direction.
const PROBE_DT: Days = Days::new(1e-7);

/// Deduplication epsilon for crossings.
const DEDUPE_EPS: Days = Days::new(1e-8);

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

/// A labelled crossing: time and direction.
///  +1 = entering the region,  −1 = exiting the region.
#[derive(Debug, Clone, Copy)]
pub struct LabeledCrossing {
    pub t: ModifiedJulianDate,
    pub direction: i32,
}

// ---------------------------------------------------------------------------
// Core: find crossings of f(t) = threshold via scan + Brent
// ---------------------------------------------------------------------------

/// Scan `period` at `step` intervals, find all roots of
/// `f(t) − threshold` using Brent's method.  Returns unsorted crossing times.
pub fn find_crossings<V, F>(
    period: Period<MJD>,
    step: Days,
    f: &F,
    threshold: Quantity<V>,
) -> Vec<ModifiedJulianDate>
where
    V: Unit,
    F: Fn(ModifiedJulianDate) -> Quantity<V>,
{
    let g = |t: Mjd| -> Quantity<V> { f(t) - threshold };

    let step_v = step;
    let t_start_v = period.start;
    let t_end_v = period.end;
    let zero = Quantity::<V>::new(0.0);

    let mut crossings = Vec::new();
    let mut t = t_start_v;
    let mut prev = g(t);
    while t < t_end_v {
        let next_t = (t + step_v).min(t_end_v);
        let next_v = g(next_t);

        if (prev < zero && next_v > zero) || (prev > zero && next_v < zero) {
            if let Some(root_t) =
                root_finding::brent_with_values(Period::new(t, next_t), prev, next_v, g)
            {
                if root_t >= t_start_v && root_t <= t_end_v {
                    crossings.push(root_t);
                }
            }
        }

        t = next_t;
        prev = next_v;
    }

    crossings
}

/// Find crossings within pre‑computed key‑time segments (e.g. between
/// successive culminations).  At most one crossing per segment is expected.
pub fn find_crossings_in_segments<V, F>(
    key_times: &[ModifiedJulianDate],
    f: &F,
    threshold: Quantity<V>,
    period: Period<MJD>,
) -> Vec<ModifiedJulianDate>
where
    V: Unit,
    F: Fn(ModifiedJulianDate) -> Quantity<V>,
{
    let g = |t: Mjd| -> Quantity<V> { f(t) - threshold };
    let t_start_v = period.start;
    let t_end_v = period.end;
    let zero = Quantity::<V>::new(0.0);

    let mut crossings = Vec::new();

    for window in key_times.windows(2) {
        let a = window[0];
        let b = window[1];
        if a >= b {
            continue;
        }

        let fa = g(a);
        let fb = g(b);

        const ROOT_EPS: f64 = 1e-12;
        if fa.abs() < ROOT_EPS {
            crossings.push(a);
            continue;
        }
        if fb.abs() < ROOT_EPS {
            crossings.push(b);
            continue;
        }

        if (fa < zero && fb > zero) || (fa > zero && fb < zero) {
            if let Some(rv) = root_finding::brent_with_values(Period::new(a, b), fa, fb, g) {
                if rv >= t_start_v && rv <= t_end_v {
                    crossings.push(rv);
                }
            }
        }
    }

    crossings
}

// ---------------------------------------------------------------------------
// Core: label crossings (entering / exiting)
// ---------------------------------------------------------------------------

/// Label each crossing as entering (+1) or exiting (−1) the above‑threshold
/// region, by probing f(t±ε).  Tangencies (same side before and after) are
/// skipped.
pub fn label_crossings<V, F>(
    crossings: &mut Vec<ModifiedJulianDate>,
    f: &F,
    threshold: Quantity<V>,
) -> Vec<LabeledCrossing>
where
    V: Unit,
    F: Fn(ModifiedJulianDate) -> Quantity<V>,
{
    crossings.sort_by(|a, b| a.partial_cmp(b).unwrap());
    crossings.dedup_by(|a, b| (*a - *b).abs() < DEDUPE_EPS);

    let is_above = |v: Quantity<V>| v > threshold;

    crossings
        .iter()
        .filter_map(|&t| {
            let tv = t;
            let before = is_above(f(tv - PROBE_DT));
            let after = is_above(f(tv + PROBE_DT));
            if !before && after {
                Some(LabeledCrossing { t, direction: 1 })
            } else if before && !after {
                Some(LabeledCrossing { t, direction: -1 })
            } else {
                None // tangency
            }
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Core: assemble "above threshold" periods from labelled crossings
// ---------------------------------------------------------------------------

/// Build time periods where `f(t) > threshold` from pre‑labelled crossings.
///
/// `labeled` must be sorted chronologically.  The caller supplies
/// `start_above` (whether f(t_start) > threshold) to handle the leading
/// edge correctly.
///
/// A midpoint validation check is performed for each candidate period.
pub fn build_above_periods<V, F>(
    labeled: &[LabeledCrossing],
    period: Period<MJD>,
    start_above: bool,
    f: &F,
    threshold: Quantity<V>,
) -> Vec<Period<MJD>>
where
    V: Unit,
    F: Fn(ModifiedJulianDate) -> Quantity<V>,
{
    let t_start = period.start;
    let t_end = period.end;
    let is_above = |v: Quantity<V>| v > threshold;
    let mut periods = Vec::new();

    if labeled.is_empty() {
        if start_above {
            return vec![Period::new(t_start, t_end)];
        }
        return Vec::new();
    }

    let mut i = 0;

    // Leading partial period: we start above and first crossing exits
    if start_above && labeled[0].direction == -1 {
        let exit_t = labeled[0].t;
        let mid_v = t_start + (exit_t - t_start) * 0.5;
        if is_above(f(mid_v)) {
            periods.push(Period::new(period.start, exit_t));
        }
        i = 1;
    }

    // Process remaining enter/exit pairs
    while i < labeled.len() {
        if labeled[i].direction == 1 {
            let enter_t = labeled[i].t;
            let exit_t = if i + 1 < labeled.len() && labeled[i + 1].direction == -1 {
                let t = labeled[i + 1].t;
                i += 2;
                t
            } else {
                i += 1;
                t_end
            };

            let mid_v = enter_t + (exit_t - enter_t) * 0.5;
            if mid_v >= t_start && mid_v <= t_end && is_above(f(mid_v)) {
                periods.push(Period::new(enter_t, exit_t));
            }
        } else {
            i += 1;
        }
    }

    periods
}

// ---------------------------------------------------------------------------
// High‑level: above‑threshold periods (scan approach)
// ---------------------------------------------------------------------------

/// Find all periods where `f(t) > threshold` in `period`,
/// using a coarse scan at `step` followed by Brent refinement and
/// crossing classification.
pub fn above_threshold_periods<V, F>(
    period: Period<MJD>,
    step: Days,
    f: &F,
    threshold: Quantity<V>,
) -> Vec<Period<MJD>>
where
    V: Unit,
    F: Fn(ModifiedJulianDate) -> Quantity<V>,
{
    let mut crossings = find_crossings(period, step, f, threshold);
    let labeled = label_crossings(&mut crossings, f, threshold);
    let start_above = f(period.start) > threshold;
    build_above_periods(&labeled, period, start_above, f, threshold)
}

/// Find all periods where `f(t) > threshold` using pre‑computed key‑time
/// segments (culmination‑based).
pub fn above_threshold_periods_segmented<V, F>(
    key_times: &[ModifiedJulianDate],
    period: Period<MJD>,
    f: &F,
    threshold: Quantity<V>,
) -> Vec<Period<MJD>>
where
    V: Unit,
    F: Fn(ModifiedJulianDate) -> Quantity<V>,
{
    let mut crossings = find_crossings_in_segments(key_times, f, threshold, period);
    let labeled = label_crossings(&mut crossings, f, threshold);
    let start_above = f(period.start) > threshold;
    build_above_periods(&labeled, period, start_above, f, threshold)
}

// ---------------------------------------------------------------------------
// High‑level: range [h_min, h_max] periods
// ---------------------------------------------------------------------------

/// Find all periods where `h_min < f(t) < h_max` in `period`.
///
/// Computed as `above(h_min) ∩ complement(above(h_max))`.
pub fn in_range_periods<V, F>(
    period: Period<MJD>,
    step: Days,
    f: &F,
    h_min: Quantity<V>,
    h_max: Quantity<V>,
) -> Vec<Period<MJD>>
where
    V: Unit,
    F: Fn(ModifiedJulianDate) -> Quantity<V>,
{
    let above_min = above_threshold_periods(period, step, f, h_min);
    let above_max = above_threshold_periods(period, step, f, h_max);
    let below_max = complement(period, &above_max);
    intersect(&above_min, &below_max)
}

/// Like [`in_range_periods`] but using key‑time segments.
pub fn in_range_periods_segmented<V, F>(
    key_times: &[ModifiedJulianDate],
    period: Period<MJD>,
    f: &F,
    h_min: Quantity<V>,
    h_max: Quantity<V>,
) -> Vec<Period<MJD>>
where
    V: Unit,
    F: Fn(ModifiedJulianDate) -> Quantity<V>,
{
    let above_min = above_threshold_periods_segmented(key_times, period, f, h_min);
    let above_max = above_threshold_periods_segmented(key_times, period, f, h_max);
    let below_max = complement(period, &above_max);
    intersect(&above_min, &below_max)
}

// ---------------------------------------------------------------------------
// Period algebra
// ---------------------------------------------------------------------------

/// Complement of `periods` within `within`.
pub fn complement(within: Period<MJD>, periods: &[Period<MJD>]) -> Vec<Period<MJD>> {
    complement_within(within, periods)
}

/// Intersection of two sorted, non‑overlapping period lists.
pub fn intersect(a: &[Period<MJD>], b: &[Period<MJD>]) -> Vec<Period<MJD>> {
    intersect_periods(a, b)
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
    fn period(a: f64, b: f64) -> Period<MJD> {
        Period::new(mjd(a), mjd(b))
    }

    fn mjd_scalar(t: Mjd) -> f64 {
        t / Days::new(1.0)
    }

    #[test]
    fn above_threshold_sine_wave() {
        // sin(2π(t+0.05)) > 0 roughly on (−0.05, 0.45) → within [0,1]: (0, ~0.45) and (~0.95, 1)
        // Use a shift to avoid exact zero at t=0
        let f = |t: Mjd| Radians::new((2.0 * std::f64::consts::PI * (mjd_scalar(t) + 0.05)).sin());
        let periods =
            above_threshold_periods(period(0.0, 1.0), Days::new(0.01), &f, Radians::new(0.0));

        // Should find at least one above-threshold period
        assert!(!periods.is_empty(), "got {:?}", periods);
        // Total above duration should be roughly 0.5
        let total = periods
            .iter()
            .fold(Days::new(0.0), |sum, p| sum + p.duration_days());
        assert!(
            (total - Days::new(0.5)).abs() < Days::new(0.05),
            "total = {total}"
        );
    }

    #[test]
    fn above_threshold_always_above() {
        let periods = above_threshold_periods(
            period(0.0, 10.0),
            Days::new(1.0),
            &|_: Mjd| Radians::new(5.0),
            Radians::new(0.0),
        );
        assert_eq!(periods.len(), 1);
        assert_eq!(periods[0], period(0.0, 10.0));
    }

    #[test]
    fn above_threshold_always_below() {
        let periods = above_threshold_periods(
            period(0.0, 10.0),
            Days::new(1.0),
            &|_: Mjd| Radians::new(-5.0),
            Radians::new(0.0),
        );
        assert!(periods.is_empty());
    }

    #[test]
    fn in_range_periods_band() {
        // sin(2π(t+0.05)) in range [-0.5, 0.5]
        // The band where |sin| < 0.5 occupies 1/3 of each cycle.
        let f = |t: Mjd| Radians::new((2.0 * std::f64::consts::PI * (mjd_scalar(t) + 0.05)).sin());
        let periods = in_range_periods(
            period(0.0, 1.0),
            Days::new(0.01),
            &f,
            Radians::new(-0.5),
            Radians::new(0.5),
        );

        assert!(!periods.is_empty(), "should find some in-range periods");

        let total = periods
            .iter()
            .fold(Days::new(0.0), |sum, p| sum + p.duration_days());
        // Analytically, sin(x) ∈ (-0.5, 0.5) for 1/3 of each full cycle ≈ 0.333
        assert!(
            total > Days::new(0.25) && total < Days::new(0.45),
            "total = {total}"
        );
    }

    #[test]
    fn complement_basic() {
        let intervals = vec![period(2.0, 4.0), period(6.0, 8.0)];
        let gaps = complement(period(0.0, 10.0), &intervals);

        assert_eq!(gaps.len(), 3);
        assert_eq!(gaps[0], period(0.0, 2.0));
        assert_eq!(gaps[1], period(4.0, 6.0));
        assert_eq!(gaps[2], period(8.0, 10.0));
    }

    #[test]
    fn complement_empty_input() {
        let gaps = complement(period(0.0, 10.0), &[]);
        assert_eq!(gaps.len(), 1);
        assert!((gaps[0].duration_days() - Days::new(10.0)).abs() < Days::new(1e-10));
    }

    #[test]
    fn complement_full_coverage() {
        let gaps = complement(period(0.0, 10.0), &[period(0.0, 10.0)]);
        assert!(gaps.is_empty());
    }

    #[test]
    fn intersect_basic() {
        let a = vec![period(0.0, 5.0), period(7.0, 10.0)];
        let b = vec![period(3.0, 8.0)];
        let result = intersect(&a, &b);

        assert_eq!(result.len(), 2);
        assert_eq!(result[0], period(3.0, 5.0));
        assert_eq!(result[1], period(7.0, 8.0));
    }

    #[test]
    fn intersect_no_overlap() {
        let a = vec![period(0.0, 2.0)];
        let b = vec![period(3.0, 5.0)];
        assert!(intersect(&a, &b).is_empty());
    }

    #[test]
    fn multiple_crossings_per_window() {
        // Fast oscillation: sin(20πt) > 0  → 10 positive intervals in [0,1]
        let f = |t: Mjd| Radians::new((20.0 * std::f64::consts::PI * mjd_scalar(t)).sin());
        let periods =
            above_threshold_periods(period(0.0, 1.0), Days::new(0.005), &f, Radians::new(0.0));
        assert!(
            periods.len() >= 8,
            "expected ≥8 intervals for fast oscillation, got {}",
            periods.len()
        );
    }

    #[test]
    fn segmented_matches_scan() {
        // Use shifted sine to avoid exact zeros at grid points
        let f = |t: Mjd| Radians::new((2.0 * std::f64::consts::PI * (mjd_scalar(t) + 0.05)).sin());

        // Build key times at 0.1 intervals
        let key_times: Vec<Mjd> = (0..=10).map(|i| mjd(i as f64 * 0.1)).collect();

        let p = period(0.0, 1.0);

        let scan = above_threshold_periods(p, Days::new(0.01), &f, Radians::new(0.0));
        let segmented = above_threshold_periods_segmented(&key_times, p, &f, Radians::new(0.0));

        assert_eq!(
            scan.len(),
            segmented.len(),
            "scan={scan:?} seg={segmented:?}"
        );
        for (s, g) in scan.iter().zip(segmented.iter()) {
            assert!((s.start - g.start).abs() < Days::new(0.05));
            assert!((s.end - g.end).abs() < Days::new(0.05));
        }
    }

    #[test]
    fn find_crossings_in_segments_basic() {
        let f = |t: Mjd| Radians::new(mjd_scalar(t) - 5.0);
        let keys: Vec<Mjd> = vec![mjd(0.0), mjd(3.0), mjd(7.0), mjd(10.0)];
        let crossings = find_crossings_in_segments(&keys, &f, Radians::new(0.0), period(0.0, 10.0));
        assert_eq!(crossings.len(), 1);
        assert!((crossings[0] - mjd(5.0)).abs() < Days::new(1e-8));
    }

    #[test]
    fn label_crossings_skips_tangency() {
        // f(t) = t^2 touches 0 at t=0 but doesn't cross
        let f = |t: Mjd| {
            let t = mjd_scalar(t);
            Radians::new(t * t)
        };
        let mut crossings = vec![mjd(0.0)];
        let labeled = label_crossings(&mut crossings, &f, Radians::new(0.0));
        // t=0 is a tangency (above on both sides), should be skipped
        assert!(labeled.is_empty(), "tangency should be skipped");
    }

    #[test]
    fn build_above_periods_from_labeled() {
        let labeled = vec![
            LabeledCrossing {
                t: mjd(2.0),
                direction: 1,
            },
            LabeledCrossing {
                t: mjd(5.0),
                direction: -1,
            },
            LabeledCrossing {
                t: mjd(7.0),
                direction: 1,
            },
            LabeledCrossing {
                t: mjd(9.0),
                direction: -1,
            },
        ];
        let f = |t: Mjd| -> Radians {
            let tv = mjd_scalar(t);
            if (tv > 2.0 && tv < 5.0) || (tv > 7.0 && tv < 9.0) {
                Radians::new(1.0)
            } else {
                Radians::new(-1.0)
            }
        };
        let periods =
            build_above_periods(&labeled, period(0.0, 10.0), false, &f, Radians::new(0.0));
        assert_eq!(periods.len(), 2);
        assert_eq!(periods[0], period(2.0, 5.0));
        assert_eq!(periods[1], period(7.0, 9.0));
    }
}
