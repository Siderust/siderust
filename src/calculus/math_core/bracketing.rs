// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Bracketing / Seeding Policies
//!
//! Strategies for producing candidate brackets where a scalar function may
//! cross a threshold or reach an extremum.
//!
//! All routines operate on `Period<ModifiedJulianDate>` time windows and
//! closures `Fn(ModifiedJulianDate) → Quantity<V>`.
//!
//! ## Provided policies
//!
//! | Policy | Description |
//! |--------|-------------|
//! | [`fixed_step_brackets`] | Uniform step across the window |
//! | [`adaptive_step_brackets`] | Narrows step near suspected events |
//! | [`extrema_based_brackets`] | Find extrema first, bracket crossings around each |

use crate::time::{ModifiedJulianDate, Period};
use qtty::*;

use super::extrema::{self, ExtremumKind};

type MJD = ModifiedJulianDate;

// ---------------------------------------------------------------------------
// Fixed‑step seeding
// ---------------------------------------------------------------------------

/// Generate uniform brackets at a fixed step size.
///
/// Returns sign‑change brackets for `f(t) − threshold`.
pub fn fixed_step_brackets<V, F>(
    period: Period<MJD>,
    step: Days,
    f: &F,
    threshold: Quantity<V>,
) -> Vec<Period<MJD>>
where
    V: Unit,
    F: Fn(ModifiedJulianDate) -> Quantity<V>,
{
    let g = |t: MJD| -> Quantity<V> { f(t) - threshold };

    let mut brackets = Vec::new();
    let mut t = period.start;
    let mut prev = g(t);

    while t < period.end {
        let next_t = (t + step).min(period.end);
        let next_v = g(next_t);

        if prev.value() * next_v.value() < 0.0 {
            brackets.push(Period::new(t, next_t));
        }

        t = next_t;
        prev = next_v;
    }

    brackets
}

// ---------------------------------------------------------------------------
// Adaptive‑step seeding
// ---------------------------------------------------------------------------

/// Generate brackets with adaptive step: starts at `initial_step`, subdivides
/// where the function changes rapidly (derivative exceeds `rate_threshold`).
///
/// `min_step` prevents infinite subdivision.
pub fn adaptive_step_brackets<V, F>(
    period: Period<MJD>,
    initial_step: Days,
    min_step: Days,
    f: &F,
    threshold: Quantity<V>,
) -> Vec<Period<MJD>>
where
    V: Unit,
    F: Fn(ModifiedJulianDate) -> Quantity<V>,
{
    let g = |t: MJD| -> Quantity<V> { f(t) - threshold };

    let mut brackets = Vec::new();

    struct Frame<V: Unit> {
        period: Period<MJD>,
        g_lo: Quantity<V>,
        g_hi: Quantity<V>,
    }

    let mut stack: Vec<Frame<V>> = Vec::new();

    // Initial pass at coarse step
    let mut t = period.start;
    let mut prev = g(t);
    while t < period.end {
        let next_t = (t + initial_step).min(period.end);
        let next_v = g(next_t);
        stack.push(Frame {
            period: Period::new(t, next_t),
            g_lo: prev,
            g_hi: next_v,
        });
        t = next_t;
        prev = next_v;
    }

    // Process stack: subdivide large steps where sign change detected
    while let Some(frame) = stack.pop() {
        if frame.g_lo.value() * frame.g_hi.value() < 0.0 {
            let width = frame.period.end - frame.period.start;
            if width <= min_step + min_step {
                brackets.push(frame.period);
            } else {
                // Subdivide to find tighter bracket
                let mid = frame.period.start + width * 0.5;
                let g_mid = g(mid);
                // Push both halves (will be processed)
                stack.push(Frame {
                    period: Period::new(frame.period.start, mid),
                    g_lo: frame.g_lo,
                    g_hi: g_mid,
                });
                stack.push(Frame {
                    period: Period::new(mid, frame.period.end),
                    g_lo: g_mid,
                    g_hi: frame.g_hi,
                });
            }
        }
    }

    brackets.sort_by(|a, b| a.start.partial_cmp(&b.start).unwrap());
    brackets
}

// ---------------------------------------------------------------------------
// Extrema‑based seeding (for satellite passes / fast movers)
// ---------------------------------------------------------------------------

/// Find extrema of `f` first, then generate crossing brackets around each
/// extremum where `f(extremum) > threshold` (or `< threshold` for minima that
/// dip below).
///
/// This is ideal for satellite pass detection: find the altitude peak of each
/// pass, then bracket the rise/set crossings on either side.
pub fn extrema_based_brackets<V, F>(
    period: Period<MJD>,
    extrema_step: Days,
    f: &F,
    threshold: Quantity<V>,
) -> Vec<Period<MJD>>
where
    V: Unit,
    F: Fn(ModifiedJulianDate) -> Quantity<V>,
{
    let extrema = extrema::find_extrema(period, extrema_step, f);

    let mut brackets = Vec::new();

    for ext in &extrema {
        if ext.kind == ExtremumKind::Maximum && ext.value.value() > threshold.value() {
            // This maximum is above threshold → there must be a rising crossing
            // before it and a setting crossing after it.
            // Search backward from the extremum for the rising crossing
            let rise_bracket = search_crossing_backward(
                Period::new(period.start, ext.t),
                f,
                threshold,
            );
            if let Some(br) = rise_bracket {
                brackets.push(br);
            }

            // Search forward from the extremum for the setting crossing
            let set_bracket = search_crossing_forward(
                Period::new(ext.t, period.end),
                f,
                threshold,
            );
            if let Some(br) = set_bracket {
                brackets.push(br);
            }
        }
    }

    // Deduplicate overlapping brackets
    brackets.sort_by(|a, b| a.start.partial_cmp(&b.start).unwrap());
    brackets.dedup_by(|a, b| {
        const TOL: Days = Days::new(1e-8); // 0.00000001 days ~ 0.86 seconds
        let dt_start = (a.start - b.start).abs();
        let dt_end = (a.end - b.end).abs();
        dt_start < TOL && dt_end < TOL
    });
    brackets
}

/// Search backward from `search_period.end` to `search_period.start` for a sign change in `f(t) − threshold`.
fn search_crossing_backward<V, F>(
    search_period: Period<MJD>,
    f: &F,
    threshold: Quantity<V>,
) -> Option<Period<MJD>>
where
    V: Unit,
    F: Fn(ModifiedJulianDate) -> Quantity<V>,
{
    let g = |t: MJD| f(t).value() - threshold.value();
    let range = search_period.duration();

    // Expanding search: start near the extremum and step backward
    let mut bracket = Period::new(search_period.end, search_period.end);
    let mut g_hi = g(bracket.end);
    let mut step = Days::new(range.value() * 0.1);
    if step.value() < 1e-10 {
        step = Days::new(range.value() * 0.5);
    }

    bracket.start = (bracket.end - step).max(search_period.start);
    let mut g_lo = g(bracket.start);

    while bracket.start > search_period.start {
        if g_lo * g_hi < 0.0 {
            return Some(bracket);
        }
        bracket.end = bracket.start;
        g_hi = g_lo;
        bracket.start = (bracket.start - step).max(search_period.start);
        g_lo = g(bracket.start);
    }

    // Check the final segment
    if g_lo * g_hi < 0.0 {
        Some(bracket)
    } else {
        None
    }
}

/// Search forward from `search_period.start` to `search_period.end` for a sign change.
fn search_crossing_forward<V, F>(
    search_period: Period<MJD>,
    f: &F,
    threshold: Quantity<V>,
) -> Option<Period<MJD>>
where
    V: Unit,
    F: Fn(ModifiedJulianDate) -> Quantity<V>,
{
    let g = |t: MJD| f(t).value() - threshold.value();
    let range = search_period.duration();

    let mut bracket = Period::new(search_period.start, search_period.start);
    let mut g_lo = g(bracket.start);
    let mut step = Days::new(range.value() * 0.1);
    if step.value() < 1e-10 {
        step = Days::new(range.value() * 0.5);
    }

    bracket.end = (bracket.start + step).min(search_period.end);
    let mut g_hi = g(bracket.end);

    while bracket.end < search_period.end {
        if g_lo * g_hi < 0.0 {
            return Some(bracket);
        }
        bracket.start = bracket.end;
        g_lo = g_hi;
        bracket.end = (bracket.end + step).min(search_period.end);
        g_hi = g(bracket.end);
    }

    if g_lo * g_hi < 0.0 {
        Some(bracket)
    } else {
        None
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use qtty::Radian;

    type Radians = Quantity<Radian>;

    fn mjd(v: f64) -> MJD {
        MJD::new(v)
    }
    fn period(a: f64, b: f64) -> Period<MJD> {
        Period::new(mjd(a), mjd(b))
    }

    #[test]
    fn fixed_step_finds_sine_crossings() {
        // Shift to avoid exact zeros at grid points
        let f = |t: MJD| Radians::new((2.0 * std::f64::consts::PI * (t.value() + 0.03)).sin());
        let brackets = fixed_step_brackets(
            period(0.0, 1.0),
            Days::new(0.05),
            &f,
            Radians::new(0.0),
        );
        // Shifted sin crosses 0 twice in [0,1] (near 0.47 and 0.97)
        assert!(
            brackets.len() >= 2,
            "got {} brackets: {:?}",
            brackets.len(),
            brackets
        );
    }

    #[test]
    fn adaptive_step_finds_crossings() {
        let f = |t: MJD| Radians::new((2.0 * std::f64::consts::PI * (t.value() + 0.03)).sin());
        let brackets = adaptive_step_brackets(
            period(0.0, 1.0),
            Days::new(0.2),
            Days::new(0.01),
            &f,
            Radians::new(0.0),
        );
        assert!(brackets.len() >= 2, "got {} brackets", brackets.len());
    }

    #[test]
    fn extrema_based_finds_satellite_pass() {
        // Simulate a satellite pass: altitude peaks at t=5 with max=30°
        let f = |t: MJD| Radians::new(-2.0 * (t.value() - 5.0).powi(2) + 30.0);
        let brackets = extrema_based_brackets(
            period(0.0, 10.0),
            Days::new(0.5),
            &f,
            Radians::new(10.0),
        );

        // Should find rise and set brackets around the peak
        assert!(brackets.len() >= 2, "got {} brackets", brackets.len());

        // Verify brackets contain actual crossings
        for br in &brackets {
            let g_lo = f(br.start).value() - 10.0;
            let g_hi = f(br.end).value() - 10.0;
            assert!(
                g_lo * g_hi <= 0.0,
                "bracket [{}, {}] doesn't contain crossing: g_lo={}, g_hi={}",
                br.start.value(),
                br.end.value(),
                g_lo,
                g_hi
            );
        }
    }

    #[test]
    fn extrema_based_no_pass_above_threshold() {
        // Satellite never reaches threshold
        let f = |t: MJD| Radians::new(-2.0 * (t.value() - 5.0).powi(2) + 5.0);
        let brackets = extrema_based_brackets(
            period(0.0, 10.0),
            Days::new(0.5),
            &f,
            Radians::new(10.0),
        );
        assert!(brackets.is_empty(), "no passes above threshold");
    }

    #[test]
    fn fixed_step_no_crossings() {
        let brackets = fixed_step_brackets(
            period(0.0, 10.0),
            Days::new(1.0),
            &|_: MJD| Radians::new(5.0),
            Radians::new(0.0),
        );
        assert!(brackets.is_empty());
    }
}
