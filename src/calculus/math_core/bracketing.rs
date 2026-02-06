// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Bracketing / Seeding Policies
//!
//! Strategies for producing candidate brackets where a scalar function may
//! cross a threshold or reach an extremum.  These are astronomy‑agnostic —
//! they only know about scalar functions of `f64`.
//!
//! ## Provided policies
//!
//! | Policy | Description |
//! |--------|-------------|
//! | [`FixedStep`] | Uniform step across the window |
//! | [`AdaptiveStep`] | Narrows step near suspected events |
//! | [`ExtremaBasedSeeds`] | Find extrema first, bracket crossings around each |

use super::extrema::{self, ExtremumKind};

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

/// A bracket: a pair `(lo, hi)` that may contain a root or extremum.
#[derive(Debug, Clone, Copy)]
pub struct Bracket {
    pub lo: f64,
    pub hi: f64,
}

// ---------------------------------------------------------------------------
// Fixed‑step seeding
// ---------------------------------------------------------------------------

/// Generate uniform brackets at a fixed step size.
///
/// Returns sign‑change brackets for `f(t) − threshold`.
pub fn fixed_step_brackets<F: Fn(f64) -> f64>(
    t_start: f64,
    t_end: f64,
    step: f64,
    f: &F,
    threshold: f64,
) -> Vec<Bracket> {
    let g = |t: f64| f(t) - threshold;
    let mut brackets = Vec::new();
    let mut t = t_start;
    let mut prev = g(t);

    while t < t_end {
        let next_t = (t + step).min(t_end);
        let next_v = g(next_t);

        if prev * next_v < 0.0 {
            brackets.push(Bracket { lo: t, hi: next_t });
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
pub fn adaptive_step_brackets<F: Fn(f64) -> f64>(
    t_start: f64,
    t_end: f64,
    initial_step: f64,
    min_step: f64,
    f: &F,
    threshold: f64,
) -> Vec<Bracket> {
    let g = |t: f64| f(t) - threshold;
    let mut brackets = Vec::new();

    struct Frame {
        lo: f64,
        hi: f64,
        g_lo: f64,
        g_hi: f64,
    }

    let mut stack: Vec<Frame> = Vec::new();

    // Initial pass at coarse step
    let mut t = t_start;
    let mut prev = g(t);
    while t < t_end {
        let next_t = (t + initial_step).min(t_end);
        let next_v = g(next_t);
        stack.push(Frame {
            lo: t,
            hi: next_t,
            g_lo: prev,
            g_hi: next_v,
        });
        t = next_t;
        prev = next_v;
    }

    // Process stack: subdivide large steps where sign change detected
    while let Some(frame) = stack.pop() {
        if frame.g_lo * frame.g_hi < 0.0 {
            let width = frame.hi - frame.lo;
            if width <= min_step * 2.0 {
                brackets.push(Bracket {
                    lo: frame.lo,
                    hi: frame.hi,
                });
            } else {
                // Subdivide to find tighter bracket
                let mid = 0.5 * (frame.lo + frame.hi);
                let g_mid = g(mid);
                // Push both halves (will be processed)
                stack.push(Frame {
                    lo: frame.lo,
                    hi: mid,
                    g_lo: frame.g_lo,
                    g_hi: g_mid,
                });
                stack.push(Frame {
                    lo: mid,
                    hi: frame.hi,
                    g_lo: g_mid,
                    g_hi: frame.g_hi,
                });
            }
        }
    }

    brackets.sort_by(|a, b| a.lo.partial_cmp(&b.lo).unwrap());
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
pub fn extrema_based_brackets<F: Fn(f64) -> f64>(
    t_start: f64,
    t_end: f64,
    extrema_step: f64,
    f: &F,
    threshold: f64,
) -> Vec<Bracket> {
    let extrema = extrema::find_extrema(t_start, t_end, extrema_step, f);

    let mut brackets = Vec::new();

    for ext in &extrema {
        if ext.kind == ExtremumKind::Maximum && ext.value > threshold {
            // This maximum is above threshold → there must be a rising crossing
            // before it and a setting crossing after it.
            // Search backward from the extremum for the rising crossing
            let rise_bracket = search_crossing_backward(t_start, ext.t, f, threshold);
            if let Some(br) = rise_bracket {
                brackets.push(br);
            }

            // Search forward from the extremum for the setting crossing
            let set_bracket = search_crossing_forward(ext.t, t_end, f, threshold);
            if let Some(br) = set_bracket {
                brackets.push(br);
            }
        }
    }

    // Deduplicate overlapping brackets
    brackets.sort_by(|a, b| a.lo.partial_cmp(&b.lo).unwrap());
    brackets.dedup_by(|a, b| (a.lo - b.lo).abs() < 1e-8 && (a.hi - b.hi).abs() < 1e-8);
    brackets
}

/// Search backward from `t_max` to `t_start` for a sign change in `f(t) − threshold`.
fn search_crossing_backward<F: Fn(f64) -> f64>(
    t_start: f64,
    t_max: f64,
    f: &F,
    threshold: f64,
) -> Option<Bracket> {
    let g = |t: f64| f(t) - threshold;

    // Expanding search: start near the extremum and step backward
    let mut hi = t_max;
    let mut g_hi = g(hi);
    let mut step = (t_max - t_start) * 0.1;
    if step < 1e-10 {
        step = (t_max - t_start) * 0.5;
    }

    let mut lo = (hi - step).max(t_start);
    let mut g_lo = g(lo);

    while lo > t_start {
        if g_lo * g_hi < 0.0 {
            return Some(Bracket { lo, hi });
        }
        hi = lo;
        g_hi = g_lo;
        lo = (lo - step).max(t_start);
        g_lo = g(lo);
    }

    // Check the final segment
    if g_lo * g_hi < 0.0 {
        Some(Bracket { lo, hi })
    } else {
        None
    }
}

/// Search forward from `t_min` to `t_end` for a sign change.
fn search_crossing_forward<F: Fn(f64) -> f64>(
    t_min: f64,
    t_end: f64,
    f: &F,
    threshold: f64,
) -> Option<Bracket> {
    let g = |t: f64| f(t) - threshold;

    let mut lo = t_min;
    let mut g_lo = g(lo);
    let mut step = (t_end - t_min) * 0.1;
    if step < 1e-10 {
        step = (t_end - t_min) * 0.5;
    }

    let mut hi = (lo + step).min(t_end);
    let mut g_hi = g(hi);

    while hi < t_end {
        if g_lo * g_hi < 0.0 {
            return Some(Bracket { lo, hi });
        }
        lo = hi;
        g_lo = g_hi;
        hi = (hi + step).min(t_end);
        g_hi = g(hi);
    }

    if g_lo * g_hi < 0.0 {
        Some(Bracket { lo, hi })
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

    #[test]
    fn fixed_step_finds_sine_crossings() {
        // Shift to avoid exact zeros at grid points
        let f = |t: f64| (2.0 * std::f64::consts::PI * (t + 0.03)).sin();
        let brackets = fixed_step_brackets(0.0, 1.0, 0.05, &f, 0.0);
        // Shifted sin crosses 0 twice in [0,1] (near 0.47 and 0.97)
        assert!(brackets.len() >= 2, "got {} brackets: {:?}", brackets.len(), brackets);
    }

    #[test]
    fn adaptive_step_finds_crossings() {
        let f = |t: f64| (2.0 * std::f64::consts::PI * (t + 0.03)).sin();
        let brackets = adaptive_step_brackets(0.0, 1.0, 0.2, 0.01, &f, 0.0);
        assert!(brackets.len() >= 2, "got {} brackets", brackets.len());
    }

    #[test]
    fn extrema_based_finds_satellite_pass() {
        // Simulate a satellite pass: altitude peaks at t=5 with max=30°
        let f = |t: f64| -2.0 * (t - 5.0).powi(2) + 30.0;
        let brackets = extrema_based_brackets(0.0, 10.0, 0.5, &f, 10.0);

        // Should find rise and set brackets around the peak
        assert!(brackets.len() >= 2, "got {} brackets", brackets.len());

        // Verify brackets contain actual crossings
        for br in &brackets {
            let g_lo = f(br.lo) - 10.0;
            let g_hi = f(br.hi) - 10.0;
            assert!(
                g_lo * g_hi <= 0.0,
                "bracket [{}, {}] doesn't contain crossing: g_lo={}, g_hi={}",
                br.lo, br.hi, g_lo, g_hi
            );
        }
    }

    #[test]
    fn extrema_based_no_pass_above_threshold() {
        // Satellite never reaches threshold
        let f = |t: f64| -2.0 * (t - 5.0).powi(2) + 5.0;
        let brackets = extrema_based_brackets(0.0, 10.0, 0.5, &f, 10.0);
        assert!(brackets.is_empty(), "no passes above threshold");
    }

    #[test]
    fn fixed_step_no_crossings() {
        let brackets = fixed_step_brackets(0.0, 10.0, 1.0, &|_| 5.0, 0.0);
        assert!(brackets.is_empty());
    }
}
