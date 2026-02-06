// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Interval Assembly — "In‑Range" Period Detection
//!
//! Assembles time intervals where a scalar function satisfies
//! `h_min ≤ f(t) ≤ h_max`, handling:
//!
//! * Multiple crossings per window (fast movers, satellites)
//! * Tangencies (touching a threshold without crossing)
//! * Always‑above / always‑below cases
//! * Single‑threshold (above only) as a special case
//!
//! All routines are astronomy‑agnostic — they operate on plain `f64`.

use super::root_finding;

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// Tiny probe offset for classifying crossing direction.
const PROBE_DT: f64 = 1e-7;

/// Deduplication epsilon for crossings.
const DEDUPE_EPS: f64 = 1e-8;

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

/// A time interval `[start, end]`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Interval {
    pub start: f64,
    pub end: f64,
}

impl Interval {
    pub fn new(start: f64, end: f64) -> Self {
        Self { start, end }
    }

    pub fn duration(&self) -> f64 {
        self.end - self.start
    }
}

/// A labelled crossing: time and direction.
///  +1 = entering the region,  −1 = exiting the region.
#[derive(Debug, Clone, Copy)]
pub struct LabeledCrossing {
    pub t: f64,
    pub direction: i32,
}

// ---------------------------------------------------------------------------
// Core: find crossings of f(t) = threshold via scan + Brent
// ---------------------------------------------------------------------------

/// Scan `[t_start, t_end]` at `step` intervals, find all roots of
/// `f(t) − threshold` using Brent's method.  Returns unsorted crossing times.
pub fn find_crossings<F: Fn(f64) -> f64>(
    t_start: f64,
    t_end: f64,
    step: f64,
    f: &F,
    threshold: f64,
) -> Vec<f64> {
    let g = |t: f64| f(t) - threshold;
    let mut crossings = Vec::new();
    let mut t = t_start;
    let mut prev = g(t);

    while t < t_end {
        let next_t = (t + step).min(t_end);
        let next_v = g(next_t);

        if prev * next_v < 0.0 {
            if let Some(root) = root_finding::brent_with_values(t, next_t, prev, next_v, &g) {
                if root >= t_start && root <= t_end {
                    crossings.push(root);
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
pub fn find_crossings_in_segments<F: Fn(f64) -> f64>(
    key_times: &[f64],
    f: &F,
    threshold: f64,
    t_start: f64,
    t_end: f64,
) -> Vec<f64> {
    let g = |t: f64| f(t) - threshold;
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

        if fa * fb < 0.0 {
            if let Some(root) = root_finding::brent_with_values(a, b, fa, fb, &g) {
                if root >= t_start && root <= t_end {
                    crossings.push(root);
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
pub fn label_crossings<F: Fn(f64) -> f64>(
    crossings: &mut Vec<f64>,
    f: &F,
    threshold: f64,
) -> Vec<LabeledCrossing> {
    crossings.sort_by(|a, b| a.partial_cmp(b).unwrap());
    crossings.dedup_by(|a, b| (*a - *b).abs() < DEDUPE_EPS);

    let is_above = |v: f64| v > threshold;

    crossings
        .iter()
        .filter_map(|&t| {
            let before = is_above(f(t - PROBE_DT));
            let after = is_above(f(t + PROBE_DT));
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

/// Build time intervals where `f(t) > threshold` from pre‑labelled crossings.
///
/// `labeled` must be sorted chronologically.  The caller supplies
/// `start_above` (whether f(t_start) > threshold) to handle the leading
/// edge correctly.
///
/// A midpoint validation check is performed for each candidate interval.
pub fn build_above_periods<F: Fn(f64) -> f64>(
    labeled: &[LabeledCrossing],
    t_start: f64,
    t_end: f64,
    start_above: bool,
    f: &F,
    threshold: f64,
) -> Vec<Interval> {
    let is_above = |v: f64| v > threshold;
    let mut periods = Vec::new();

    if labeled.is_empty() {
        if start_above {
            return vec![Interval::new(t_start, t_end)];
        }
        return Vec::new();
    }

    let mut i = 0;

    // Leading partial interval: we start above and first crossing exits
    if start_above && labeled[0].direction == -1 {
        let exit_t = labeled[0].t;
        let mid = 0.5 * (t_start + exit_t);
        if is_above(f(mid)) {
            periods.push(Interval::new(t_start, exit_t));
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

            let mid = 0.5 * (enter_t + exit_t);
            if mid >= t_start && mid <= t_end && is_above(f(mid)) {
                periods.push(Interval::new(enter_t, exit_t));
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

/// Find all intervals where `f(t) > threshold` in `[t_start, t_end]`,
/// using a coarse scan at `step` followed by Brent refinement and
/// crossing classification.
pub fn above_threshold_periods<F: Fn(f64) -> f64>(
    t_start: f64,
    t_end: f64,
    step: f64,
    f: &F,
    threshold: f64,
) -> Vec<Interval> {
    let mut crossings = find_crossings(t_start, t_end, step, f, threshold);
    let labeled = label_crossings(&mut crossings, f, threshold);
    let start_above = f(t_start) > threshold;
    build_above_periods(&labeled, t_start, t_end, start_above, f, threshold)
}

/// Find all intervals where `f(t) > threshold` using pre‑computed key‑time
/// segments (culmination‑based).
pub fn above_threshold_periods_segmented<F: Fn(f64) -> f64>(
    key_times: &[f64],
    t_start: f64,
    t_end: f64,
    f: &F,
    threshold: f64,
) -> Vec<Interval> {
    let mut crossings = find_crossings_in_segments(key_times, f, threshold, t_start, t_end);
    let labeled = label_crossings(&mut crossings, f, threshold);
    let start_above = f(t_start) > threshold;
    build_above_periods(&labeled, t_start, t_end, start_above, f, threshold)
}

// ---------------------------------------------------------------------------
// High‑level: range [h_min, h_max] periods
// ---------------------------------------------------------------------------

/// Find all intervals where `h_min < f(t) < h_max` in `[t_start, t_end]`.
///
/// Computed as `above(h_min) ∩ complement(above(h_max))`.
pub fn in_range_periods<F: Fn(f64) -> f64>(
    t_start: f64,
    t_end: f64,
    step: f64,
    f: &F,
    h_min: f64,
    h_max: f64,
) -> Vec<Interval> {
    let above_min = above_threshold_periods(t_start, t_end, step, f, h_min);
    let above_max = above_threshold_periods(t_start, t_end, step, f, h_max);
    let below_max = complement(t_start, t_end, &above_max);
    intersect(&above_min, &below_max)
}

/// Like [`in_range_periods`] but using key‑time segments.
pub fn in_range_periods_segmented<F: Fn(f64) -> f64>(
    key_times: &[f64],
    t_start: f64,
    t_end: f64,
    f: &F,
    h_min: f64,
    h_max: f64,
) -> Vec<Interval> {
    let above_min = above_threshold_periods_segmented(key_times, t_start, t_end, f, h_min);
    let above_max = above_threshold_periods_segmented(key_times, t_start, t_end, f, h_max);
    let below_max = complement(t_start, t_end, &above_max);
    intersect(&above_min, &below_max)
}

// ---------------------------------------------------------------------------
// Interval algebra
// ---------------------------------------------------------------------------

/// Complement of `intervals` within `[t_start, t_end]`.
pub fn complement(t_start: f64, t_end: f64, intervals: &[Interval]) -> Vec<Interval> {
    let mut gaps = Vec::new();
    let mut cursor = t_start;
    for iv in intervals {
        if iv.start > cursor {
            gaps.push(Interval::new(cursor, iv.start));
        }
        if iv.end > cursor {
            cursor = iv.end;
        }
    }
    if cursor < t_end {
        gaps.push(Interval::new(cursor, t_end));
    }
    gaps
}

/// Intersection of two sorted, non‑overlapping interval lists.
pub fn intersect(a: &[Interval], b: &[Interval]) -> Vec<Interval> {
    let mut result = Vec::new();
    let (mut i, mut j) = (0, 0);
    while i < a.len() && j < b.len() {
        let start = a[i].start.max(b[j].start);
        let end = a[i].end.min(b[j].end);
        if start < end {
            result.push(Interval::new(start, end));
        }
        if a[i].end <= b[j].end {
            i += 1;
        } else {
            j += 1;
        }
    }
    result
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn above_threshold_sine_wave() {
        // sin(2π(t+0.05)) > 0 roughly on (−0.05, 0.45) → within [0,1]: (0, ~0.45) and (~0.95, 1)
        // Use a shift to avoid exact zero at t=0
        let f = |t: f64| (2.0 * std::f64::consts::PI * (t + 0.05)).sin();
        let periods = above_threshold_periods(0.0, 1.0, 0.01, &f, 0.0);

        // Should find at least one above-threshold period
        assert!(!periods.is_empty(), "got {:?}", periods);
        // Total above duration should be roughly 0.5
        let total: f64 = periods.iter().map(|p| p.duration()).sum();
        assert!((total - 0.5).abs() < 0.05, "total = {total}");
    }

    #[test]
    fn above_threshold_always_above() {
        let periods = above_threshold_periods(0.0, 10.0, 1.0, &|_| 5.0, 0.0);
        assert_eq!(periods.len(), 1);
        assert!((periods[0].start - 0.0).abs() < 1e-10);
        assert!((periods[0].end - 10.0).abs() < 1e-10);
    }

    #[test]
    fn above_threshold_always_below() {
        let periods = above_threshold_periods(0.0, 10.0, 1.0, &|_| -5.0, 0.0);
        assert!(periods.is_empty());
    }

    #[test]
    fn in_range_periods_band() {
        // sin(2π(t+0.05)) in range [-0.5, 0.5]
        // The band where |sin| < 0.5 occupies 1/3 of each cycle.
        let f = |t: f64| (2.0 * std::f64::consts::PI * (t + 0.05)).sin();
        let periods = in_range_periods(0.0, 1.0, 0.01, &f, -0.5, 0.5);

        assert!(!periods.is_empty(), "should find some in-range periods");

        let total: f64 = periods.iter().map(|p| p.duration()).sum();
        // Analytically, sin(x) ∈ (-0.5, 0.5) for 1/3 of each full cycle ≈ 0.333
        assert!(total > 0.25 && total < 0.45, "total = {total}");
    }

    #[test]
    fn complement_basic() {
        let intervals = vec![Interval::new(2.0, 4.0), Interval::new(6.0, 8.0)];
        let gaps = complement(0.0, 10.0, &intervals);

        assert_eq!(gaps.len(), 3);
        assert!((gaps[0].start - 0.0).abs() < 1e-10);
        assert!((gaps[0].end - 2.0).abs() < 1e-10);
        assert!((gaps[1].start - 4.0).abs() < 1e-10);
        assert!((gaps[1].end - 6.0).abs() < 1e-10);
        assert!((gaps[2].start - 8.0).abs() < 1e-10);
        assert!((gaps[2].end - 10.0).abs() < 1e-10);
    }

    #[test]
    fn complement_empty_input() {
        let gaps = complement(0.0, 10.0, &[]);
        assert_eq!(gaps.len(), 1);
        assert!((gaps[0].duration() - 10.0).abs() < 1e-10);
    }

    #[test]
    fn complement_full_coverage() {
        let gaps = complement(0.0, 10.0, &[Interval::new(0.0, 10.0)]);
        assert!(gaps.is_empty());
    }

    #[test]
    fn intersect_basic() {
        let a = vec![Interval::new(0.0, 5.0), Interval::new(7.0, 10.0)];
        let b = vec![Interval::new(3.0, 8.0)];
        let result = intersect(&a, &b);

        assert_eq!(result.len(), 2);
        assert!((result[0].start - 3.0).abs() < 1e-10);
        assert!((result[0].end - 5.0).abs() < 1e-10);
        assert!((result[1].start - 7.0).abs() < 1e-10);
        assert!((result[1].end - 8.0).abs() < 1e-10);
    }

    #[test]
    fn intersect_no_overlap() {
        let a = vec![Interval::new(0.0, 2.0)];
        let b = vec![Interval::new(3.0, 5.0)];
        assert!(intersect(&a, &b).is_empty());
    }

    #[test]
    fn multiple_crossings_per_window() {
        // Fast oscillation: sin(20πt) > 0  → 10 positive intervals in [0,1]
        let f = |t: f64| (20.0 * std::f64::consts::PI * t).sin();
        let periods = above_threshold_periods(0.0, 1.0, 0.005, &f, 0.0);
        assert!(
            periods.len() >= 8,
            "expected ≥8 intervals for fast oscillation, got {}",
            periods.len()
        );
    }

    #[test]
    fn segmented_matches_scan() {
        // Use shifted sine to avoid exact zeros at grid points
        let f = |t: f64| (2.0 * std::f64::consts::PI * (t + 0.05)).sin();

        // Build key times at 0.1 intervals
        let key_times: Vec<f64> = (0..=10).map(|i| i as f64 * 0.1).collect();

        let scan = above_threshold_periods(0.0, 1.0, 0.01, &f, 0.0);
        let segmented = above_threshold_periods_segmented(&key_times, 0.0, 1.0, &f, 0.0);

        assert_eq!(scan.len(), segmented.len(), "scan={scan:?} seg={segmented:?}");
        for (s, g) in scan.iter().zip(segmented.iter()) {
            assert!((s.start - g.start).abs() < 0.05);
            assert!((s.end - g.end).abs() < 0.05);
        }
    }

    #[test]
    fn find_crossings_in_segments_basic() {
        let f = |t: f64| t - 5.0;
        let keys = vec![0.0, 3.0, 7.0, 10.0];
        let crossings = find_crossings_in_segments(&keys, &f, 0.0, 0.0, 10.0);
        assert_eq!(crossings.len(), 1);
        assert!((crossings[0] - 5.0).abs() < 1e-8);
    }

    #[test]
    fn label_crossings_skips_tangency() {
        // f(t) = t^2 touches 0 at t=0 but doesn't cross
        let f = |t: f64| t * t;
        let mut crossings = vec![0.0];
        let labeled = label_crossings(&mut crossings, &f, 0.0);
        // t=0 is a tangency (above on both sides), should be skipped
        assert!(labeled.is_empty(), "tangency should be skipped");
    }

    #[test]
    fn build_above_periods_from_labeled() {
        let labeled = vec![
            LabeledCrossing { t: 2.0, direction: 1 },
            LabeledCrossing { t: 5.0, direction: -1 },
            LabeledCrossing { t: 7.0, direction: 1 },
            LabeledCrossing { t: 9.0, direction: -1 },
        ];
        let f = |t: f64| {
            if (t > 2.0 && t < 5.0) || (t > 7.0 && t < 9.0) {
                1.0
            } else {
                -1.0
            }
        };
        let periods = build_above_periods(&labeled, 0.0, 10.0, false, &f, 0.0);
        assert_eq!(periods.len(), 2);
        assert!((periods[0].start - 2.0).abs() < 1e-10);
        assert!((periods[0].end - 5.0).abs() < 1e-10);
        assert!((periods[1].start - 7.0).abs() < 1e-10);
        assert!((periods[1].end - 9.0).abs() < 1e-10);
    }
}
