// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Untyped numerical kernels for sampled spectra
//!
//! ## Scientific scope
//!
//! Synthetic photometry, atmospheric extinction, and night-sky
//! brightness pipelines reduce, at their innermost loop, to two
//! operations on a tabulated spectral function: piecewise-linear
//! interpolation at a query wavelength, and trapezoidal integration of
//! a (possibly weighted) product over a wavelength interval. The exact
//! semantics of those operations — boundary handling, monotonicity
//! guarantees, weighting — must match the reference pipelines
//! (`numpy.interp`, the SVO Filter Profile Service convention,
//! Bessell-style synthetic magnitudes) bit-for-bit, otherwise published
//! magnitudes drift by tens of millimagnitudes between toolchains.
//!
//! These kernels are the canonical implementation of those operations
//! in the crate. Both the typed [`crate::spectra::SampledSpectrum`]
//! surface and the legacy NSB code path call into them, so the two
//! routes return identical numerical results.
//!
//! ## Technical scope
//!
//! All entry points operate on plain `&[f64]` slices and return either
//! `f64` or `Result<f64, SpectrumError>`:
//!
//! - [`interp_linear`] — piecewise-linear interpolation honouring the
//!   given [`OutOfRange`] policy (matches `numpy.interp` exactly when
//!   `oor == ClampToEndpoints`).
//! - [`interp_nearest`] — nearest-neighbour interpolation with deterministic
//!   half-way tie-breaking toward the *lower* index (matches
//!   `scipy.interpolate.interp1d(kind="nearest")`).
//! - [`interp_step_left`] / [`interp_step_right`] — left- and
//!   right-continuous piecewise-constant ("step") interpolation
//!   (matches `scipy.interpolate.interp1d(kind="previous" | "next")`).
//! - [`CubicSplineCoeffs::natural`] / [`interp_cubic_spline`] — natural
//!   cubic-spline interpolation. The second-derivative coefficients are
//!   computed once in `O(n)` via the tridiagonal Thomas algorithm and then
//!   reused for every query.
//! - [`trapz`] — trapezoidal integral over the entire sampled domain.
//! - [`trapz_range`] — trapezoidal integral restricted to a sub-interval.
//! - [`trapz_weighted`] — trapezoidal integral of `f(x) · w(x)` against
//!   a separately tabulated weight, with policy-controlled boundary
//!   handling on each input.
//! - [`validate`] — pre-condition check (length match, ≥ 2 samples,
//!   strict monotonic increase).
//!
//! Pre-conditions are debug-asserted in release builds; the typed
//! constructor in [`crate::spectra::sampled`] is the safe entry point.
//!
//! ## References
//!
//! - Atkinson, K. E. (1989). *An Introduction to Numerical Analysis*,
//!   2nd ed., §5.2 (composite trapezoidal rule). John Wiley & Sons.
//!   ISBN 978-0-471-62489-9.
//! - NumPy developers. *numpy.interp* documentation
//!   (boundary-handling semantics).
//! - Press, W. H., Teukolsky, S. A., Vetterling, W. T., Flannery, B. P.
//!   (1992). *Numerical Recipes in C*, 2nd ed., §3.1, §4.1.
//!   Cambridge University Press.

use super::interp::{Interpolation, OutOfRange};
use super::SpectrumError;

/// Linear-interpolated value with the given out-of-range policy.
///
/// Matches `numpy.interp` exactly when `oor == OutOfRange::ClampToEndpoints`.
#[inline]
pub fn interp_linear(
    xs: &[f64],
    ys: &[f64],
    x: f64,
    oor: OutOfRange,
) -> Result<f64, SpectrumError> {
    debug_assert_eq!(xs.len(), ys.len());
    debug_assert!(xs.len() >= 2);
    let lo = xs[0];
    let hi = *xs.last().unwrap();
    if x <= lo {
        return match oor {
            OutOfRange::ClampToEndpoints => Ok(ys[0]),
            OutOfRange::Zero => Ok(0.0),
            OutOfRange::Error if x < lo => Err(SpectrumError::OutOfRange { x, lo, hi }),
            OutOfRange::Error => Ok(ys[0]),
        };
    }
    if x >= hi {
        return match oor {
            OutOfRange::ClampToEndpoints => Ok(*ys.last().unwrap()),
            OutOfRange::Zero => Ok(0.0),
            OutOfRange::Error if x > hi => Err(SpectrumError::OutOfRange { x, lo, hi }),
            OutOfRange::Error => Ok(*ys.last().unwrap()),
        };
    }
    let i = xs.partition_point(|&xi| xi <= x);
    let (x0, x1) = (xs[i - 1], xs[i]);
    let (y0, y1) = (ys[i - 1], ys[i]);
    let t = (x - x0) / (x1 - x0);
    Ok(y0 + t * (y1 - y0))
}

/// Nearest-neighbour interpolated value with the given out-of-range policy.
///
/// Ties (exact midpoint between adjacent samples) resolve to the *lower*
/// index, matching `scipy.interpolate.interp1d(kind="nearest")`.
#[inline]
pub fn interp_nearest(
    xs: &[f64],
    ys: &[f64],
    x: f64,
    oor: OutOfRange,
) -> Result<f64, SpectrumError> {
    debug_assert_eq!(xs.len(), ys.len());
    debug_assert!(xs.len() >= 2);
    let lo = xs[0];
    let hi = *xs.last().unwrap();
    if x <= lo {
        if x == lo {
            return Ok(ys[0]);
        }
        return out_of_range_low(ys, x, lo, hi, oor);
    }
    if x >= hi {
        if x == hi {
            return Ok(*ys.last().unwrap());
        }
        return out_of_range_high(ys, x, lo, hi, oor);
    }
    let i = xs.partition_point(|&xi| xi <= x);
    let (x0, x1) = (xs[i - 1], xs[i]);
    let mid = 0.5 * (x0 + x1);
    Ok(if x <= mid { ys[i - 1] } else { ys[i] })
}

/// Left-continuous piecewise-constant ("previous") interpolation:
/// `y(x) = ys[i-1]` for `x` in `[xs[i-1], xs[i])`.
#[inline]
pub fn interp_step_left(
    xs: &[f64],
    ys: &[f64],
    x: f64,
    oor: OutOfRange,
) -> Result<f64, SpectrumError> {
    debug_assert_eq!(xs.len(), ys.len());
    debug_assert!(xs.len() >= 2);
    let lo = xs[0];
    let hi = *xs.last().unwrap();
    if x < lo {
        return out_of_range_low(ys, x, lo, hi, oor);
    }
    if x > hi {
        return out_of_range_high(ys, x, lo, hi, oor);
    }
    if x == hi {
        return Ok(*ys.last().unwrap());
    }
    // partition_point finds the first index where xi > x; previous sample owns x.
    let i = xs.partition_point(|&xi| xi <= x);
    if i == 0 {
        Ok(ys[0])
    } else {
        Ok(ys[i - 1])
    }
}

/// Right-continuous piecewise-constant ("next") interpolation:
/// `y(x) = ys[i]` for `x` in `(xs[i-1], xs[i]]`.
#[inline]
pub fn interp_step_right(
    xs: &[f64],
    ys: &[f64],
    x: f64,
    oor: OutOfRange,
) -> Result<f64, SpectrumError> {
    debug_assert_eq!(xs.len(), ys.len());
    debug_assert!(xs.len() >= 2);
    let lo = xs[0];
    let hi = *xs.last().unwrap();
    if x < lo {
        return out_of_range_low(ys, x, lo, hi, oor);
    }
    if x > hi {
        return out_of_range_high(ys, x, lo, hi, oor);
    }
    if x == lo {
        return Ok(ys[0]);
    }
    let i = xs.partition_point(|&xi| xi < x);
    Ok(ys[i])
}

#[inline]
fn out_of_range_low(
    ys: &[f64],
    x: f64,
    lo: f64,
    hi: f64,
    oor: OutOfRange,
) -> Result<f64, SpectrumError> {
    match oor {
        OutOfRange::ClampToEndpoints => Ok(ys[0]),
        OutOfRange::Zero => Ok(0.0),
        OutOfRange::Error => Err(SpectrumError::OutOfRange { x, lo, hi }),
    }
}

#[inline]
fn out_of_range_high(
    ys: &[f64],
    x: f64,
    lo: f64,
    hi: f64,
    oor: OutOfRange,
) -> Result<f64, SpectrumError> {
    match oor {
        OutOfRange::ClampToEndpoints => Ok(*ys.last().unwrap()),
        OutOfRange::Zero => Ok(0.0),
        OutOfRange::Error => Err(SpectrumError::OutOfRange { x, lo, hi }),
    }
}

/// Precomputed natural-cubic-spline second derivatives.
///
/// Constructed once in `O(n)` time via the tridiagonal Thomas algorithm
/// and reused for every query through [`interp_cubic_spline`].
#[derive(Debug, Clone)]
pub struct CubicSplineCoeffs {
    /// Second derivatives `y''(xs[i])`, co-indexed with the source samples.
    pub y2: Vec<f64>,
}

impl CubicSplineCoeffs {
    /// Build the natural cubic spline coefficients for the given samples.
    ///
    /// The "natural" spline imposes `y''(x_0) = y''(x_{n-1}) = 0`.
    ///
    /// Returns [`SpectrumError::TooFewSamples`] if `xs.len() < 2`.
    pub fn natural(xs: &[f64], ys: &[f64]) -> Result<Self, SpectrumError> {
        debug_assert_eq!(xs.len(), ys.len());
        let n = xs.len();
        if n < 2 {
            return Err(SpectrumError::TooFewSamples(n));
        }
        let mut y2 = vec![0.0; n];
        if n == 2 {
            // Degenerates to linear interpolation; natural BCs already zero.
            return Ok(Self { y2 });
        }
        // Thomas-algorithm scratch.
        let mut u = vec![0.0; n];
        for i in 1..n - 1 {
            let sig = (xs[i] - xs[i - 1]) / (xs[i + 1] - xs[i - 1]);
            let p = sig * y2[i - 1] + 2.0;
            y2[i] = (sig - 1.0) / p;
            let num = (ys[i + 1] - ys[i]) / (xs[i + 1] - xs[i])
                - (ys[i] - ys[i - 1]) / (xs[i] - xs[i - 1]);
            u[i] = (6.0 * num / (xs[i + 1] - xs[i - 1]) - sig * u[i - 1]) / p;
        }
        for k in (1..n - 1).rev() {
            y2[k] = y2[k] * y2[k + 1] + u[k];
        }
        Ok(Self { y2 })
    }
}

/// Evaluate a natural cubic spline at `x` using precomputed coefficients.
#[inline]
pub fn interp_cubic_spline(
    xs: &[f64],
    ys: &[f64],
    coeffs: &CubicSplineCoeffs,
    x: f64,
    oor: OutOfRange,
) -> Result<f64, SpectrumError> {
    debug_assert_eq!(xs.len(), ys.len());
    debug_assert_eq!(xs.len(), coeffs.y2.len());
    debug_assert!(xs.len() >= 2);
    let lo = xs[0];
    let hi = *xs.last().unwrap();
    if x <= lo {
        if x == lo {
            return Ok(ys[0]);
        }
        return out_of_range_low(ys, x, lo, hi, oor);
    }
    if x >= hi {
        if x == hi {
            return Ok(*ys.last().unwrap());
        }
        return out_of_range_high(ys, x, lo, hi, oor);
    }
    let i = xs.partition_point(|&xi| xi <= x);
    let (x0, x1) = (xs[i - 1], xs[i]);
    let (y0, y1) = (ys[i - 1], ys[i]);
    let (y2_0, y2_1) = (coeffs.y2[i - 1], coeffs.y2[i]);
    let h = x1 - x0;
    let a = (x1 - x) / h;
    let b = (x - x0) / h;
    Ok(a * y0 + b * y1 + ((a * a * a - a) * y2_0 + (b * b * b - b) * y2_1) * (h * h) / 6.0)
}

/// Dispatches an interpolation policy to its implementation.
///
/// For [`Interpolation::CubicSpline`] this builds the coefficients on every
/// call. Callers that evaluate many points on the same table should hold a
/// [`CubicSplineCoeffs`] and invoke [`interp_cubic_spline`] directly; the
/// typed [`crate::spectra::SampledSpectrum`] does this automatically.
#[inline]
pub fn interp(
    xs: &[f64],
    ys: &[f64],
    x: f64,
    interp: Interpolation,
    oor: OutOfRange,
) -> Result<f64, SpectrumError> {
    match interp {
        Interpolation::Linear => interp_linear(xs, ys, x, oor),
        Interpolation::Nearest => interp_nearest(xs, ys, x, oor),
        Interpolation::PiecewiseConstantLeft => interp_step_left(xs, ys, x, oor),
        Interpolation::PiecewiseConstantRight => interp_step_right(xs, ys, x, oor),
        Interpolation::CubicSpline => {
            let coeffs = CubicSplineCoeffs::natural(xs, ys)?;
            interp_cubic_spline(xs, ys, &coeffs, x, oor)
        }
    }
}

/// Trapezoidal integral over the full sampled domain.
#[inline]
pub fn trapz(xs: &[f64], ys: &[f64]) -> f64 {
    debug_assert_eq!(xs.len(), ys.len());
    let mut s = 0.0;
    for i in 1..xs.len() {
        s += 0.5 * (ys[i] + ys[i - 1]) * (xs[i] - xs[i - 1]);
    }
    s
}

/// Trapezoidal integral over `[lo, hi]` using linear interpolation at the
/// (possibly partial) endpoints.
///
/// This matches the behaviour of NSB's `Spectrum::integrate_range`: each
/// internal segment that overlaps `[lo, hi]` is clipped, the endpoint values
/// are recomputed via linear interpolation against `OutOfRange::ClampToEndpoints`,
/// and the trapezoidal area of the clipped segment is added to the sum.
#[inline]
pub fn trapz_range(xs: &[f64], ys: &[f64], lo: f64, hi: f64) -> f64 {
    debug_assert_eq!(xs.len(), ys.len());
    let mut s = 0.0;
    for i in 1..xs.len() {
        let (a, b) = (xs[i - 1], xs[i]);
        if b < lo || a > hi {
            continue;
        }
        let lo_clip = a.max(lo);
        let hi_clip = b.min(hi);
        let ya = interp_linear(xs, ys, lo_clip, OutOfRange::ClampToEndpoints).unwrap();
        let yb = interp_linear(xs, ys, hi_clip, OutOfRange::ClampToEndpoints).unwrap();
        s += 0.5 * (ya + yb) * (hi_clip - lo_clip);
    }
    s
}

/// Trapezoidal integral of `source(x) * weight(x)` over `weight`'s grid.
///
/// `source` is sampled (with linear interpolation and endpoint clamping) at
/// each `weight` x-value. This is the standard "filter integral" used in
/// photometric and atmospheric pipelines.
#[inline]
pub fn trapz_weighted(
    source_xs: &[f64],
    source_ys: &[f64],
    weight_xs: &[f64],
    weight_ys: &[f64],
) -> f64 {
    debug_assert_eq!(source_xs.len(), source_ys.len());
    debug_assert_eq!(weight_xs.len(), weight_ys.len());
    let mut sum = 0.0;
    for i in 1..weight_xs.len() {
        let a = weight_xs[i - 1];
        let b = weight_xs[i];
        let fa = interp_linear(source_xs, source_ys, a, OutOfRange::ClampToEndpoints).unwrap()
            * weight_ys[i - 1];
        let fb = interp_linear(source_xs, source_ys, b, OutOfRange::ClampToEndpoints).unwrap()
            * weight_ys[i];
        sum += 0.5 * (fa + fb) * (b - a);
    }
    sum
}

/// Validate that `xs` and `ys` form a well-formed sampled spectrum.
pub fn validate(xs: &[f64], ys: &[f64]) -> Result<(), SpectrumError> {
    if xs.len() != ys.len() {
        return Err(SpectrumError::LengthMismatch {
            xs: xs.len(),
            ys: ys.len(),
        });
    }
    if xs.len() < 2 {
        return Err(SpectrumError::TooFewSamples(xs.len()));
    }
    for i in 1..xs.len() {
        if xs[i].partial_cmp(&xs[i - 1]) != Some(std::cmp::Ordering::Greater) {
            return Err(SpectrumError::NotMonotonic { index: i });
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    fn linear() -> (Vec<f64>, Vec<f64>) {
        // y(x) = 2x + 1 sampled on [0, 1, 2, 3, 4]
        let xs: Vec<f64> = (0..5).map(|i| i as f64).collect();
        let ys: Vec<f64> = xs.iter().map(|x| 2.0 * x + 1.0).collect();
        (xs, ys)
    }

    #[test]
    fn linear_interp_matches_analytic() {
        let (xs, ys) = linear();
        for x in [0.0, 0.5, 1.5, 2.25, 4.0] {
            let want = 2.0 * x + 1.0;
            let got = interp_linear(&xs, &ys, x, OutOfRange::ClampToEndpoints).unwrap();
            assert_abs_diff_eq!(got, want, epsilon = 1e-12);
        }
    }

    #[test]
    fn clamp_below_returns_endpoint() {
        let (xs, ys) = linear();
        let v = interp_linear(&xs, &ys, -10.0, OutOfRange::ClampToEndpoints).unwrap();
        assert_eq!(v, ys[0]);
    }

    #[test]
    fn zero_below_returns_zero() {
        let (xs, ys) = linear();
        let v = interp_linear(&xs, &ys, -10.0, OutOfRange::Zero).unwrap();
        assert_eq!(v, 0.0);
    }

    #[test]
    fn error_below_errors() {
        let (xs, ys) = linear();
        let err = interp_linear(&xs, &ys, -10.0, OutOfRange::Error).unwrap_err();
        assert!(matches!(err, SpectrumError::OutOfRange { .. }));
    }

    #[test]
    fn trapz_constant_function() {
        let xs: Vec<f64> = (0..11).map(|i| i as f64).collect();
        let ys = vec![3.0; xs.len()];
        // ∫_0^10 3 dx = 30
        assert_abs_diff_eq!(trapz(&xs, &ys), 30.0, epsilon = 1e-12);
    }

    #[test]
    fn trapz_linear_function_exact() {
        // Trapezoidal rule is exact for linear integrands.
        let (xs, ys) = linear();
        // ∫_0^4 (2x+1) dx = x²+x = 16 + 4 = 20
        assert_abs_diff_eq!(trapz(&xs, &ys), 20.0, epsilon = 1e-12);
    }

    #[test]
    fn trapz_range_clips_correctly() {
        let (xs, ys) = linear();
        // ∫_{0.5}^{3.5} (2x+1) dx = x²+x evaluated 3.5 - 0.5
        //                        = (12.25 + 3.5) - (0.25 + 0.5) = 15.0
        assert_abs_diff_eq!(trapz_range(&xs, &ys, 0.5, 3.5), 15.0, epsilon = 1e-12);
    }

    #[test]
    fn trapz_weighted_against_analytic_rectangle() {
        // source(x) = x on [0, 4], weight is unity rectangle on [1, 3].
        let sx: Vec<f64> = (0..=4).map(|i| i as f64).collect();
        let sy: Vec<f64> = sx.clone();
        let wx = vec![1.0, 3.0];
        let wy = vec![1.0, 1.0];
        // ∫_1^3 x dx = (9-1)/2 = 4
        assert_abs_diff_eq!(trapz_weighted(&sx, &sy, &wx, &wy), 4.0, epsilon = 1e-12);
    }

    #[test]
    fn validate_rejects_non_monotonic() {
        let xs = vec![0.0, 1.0, 0.5, 2.0];
        let ys = vec![0.0; 4];
        assert!(matches!(
            validate(&xs, &ys),
            Err(SpectrumError::NotMonotonic { index: 2 })
        ));
    }

    #[test]
    fn validate_rejects_length_mismatch() {
        let xs = vec![0.0, 1.0];
        let ys = vec![0.0];
        assert!(matches!(
            validate(&xs, &ys),
            Err(SpectrumError::LengthMismatch { xs: 2, ys: 1 })
        ));
    }

    // -------- Nearest-neighbour --------

    #[test]
    fn nearest_exact_samples() {
        let (xs, ys) = linear();
        for (i, x) in xs.iter().enumerate() {
            let v = interp_nearest(&xs, &ys, *x, OutOfRange::ClampToEndpoints).unwrap();
            assert_eq!(v, ys[i]);
        }
    }

    #[test]
    fn nearest_ties_resolve_to_lower_index() {
        let xs = vec![0.0_f64, 1.0, 2.0];
        let ys = vec![10.0_f64, 20.0, 30.0];
        // Exact midpoint between xs[0] and xs[1] is 0.5; tie -> lower index (ys[0]).
        let v = interp_nearest(&xs, &ys, 0.5, OutOfRange::ClampToEndpoints).unwrap();
        assert_eq!(v, 10.0);
    }

    #[test]
    fn nearest_off_midpoint() {
        let xs = vec![0.0_f64, 1.0, 2.0];
        let ys = vec![10.0_f64, 20.0, 30.0];
        let v = interp_nearest(&xs, &ys, 0.51, OutOfRange::ClampToEndpoints).unwrap();
        assert_eq!(v, 20.0);
        let v = interp_nearest(&xs, &ys, 0.49, OutOfRange::ClampToEndpoints).unwrap();
        assert_eq!(v, 10.0);
    }

    #[test]
    fn nearest_out_of_range_zero() {
        let (xs, ys) = linear();
        assert_eq!(
            interp_nearest(&xs, &ys, -5.0, OutOfRange::Zero).unwrap(),
            0.0
        );
        assert_eq!(
            interp_nearest(&xs, &ys, 50.0, OutOfRange::Zero).unwrap(),
            0.0
        );
    }

    // -------- Step (left / right) --------

    #[test]
    fn step_left_discontinuity() {
        let xs = vec![0.0_f64, 1.0, 2.0];
        let ys = vec![10.0_f64, 20.0, 30.0];
        // Left-continuous: y(0.999) = 10, y(1.0) = 20.
        assert_eq!(
            interp_step_left(&xs, &ys, 0.999, OutOfRange::ClampToEndpoints).unwrap(),
            10.0
        );
        assert_eq!(
            interp_step_left(&xs, &ys, 1.0, OutOfRange::ClampToEndpoints).unwrap(),
            20.0
        );
        assert_eq!(
            interp_step_left(&xs, &ys, 2.0, OutOfRange::ClampToEndpoints).unwrap(),
            30.0
        );
    }

    #[test]
    fn step_right_discontinuity() {
        let xs = vec![0.0_f64, 1.0, 2.0];
        let ys = vec![10.0_f64, 20.0, 30.0];
        // Right-continuous: y(0.0) = 10, y(0.001) = 20, y(1.0) = 20, y(1.001) = 30.
        assert_eq!(
            interp_step_right(&xs, &ys, 0.0, OutOfRange::ClampToEndpoints).unwrap(),
            10.0
        );
        assert_eq!(
            interp_step_right(&xs, &ys, 0.001, OutOfRange::ClampToEndpoints).unwrap(),
            20.0
        );
        assert_eq!(
            interp_step_right(&xs, &ys, 1.0, OutOfRange::ClampToEndpoints).unwrap(),
            20.0
        );
        assert_eq!(
            interp_step_right(&xs, &ys, 1.001, OutOfRange::ClampToEndpoints).unwrap(),
            30.0
        );
    }

    #[test]
    fn step_left_out_of_range_error() {
        let xs = vec![0.0_f64, 1.0];
        let ys = vec![0.0_f64, 1.0];
        assert!(matches!(
            interp_step_left(&xs, &ys, -0.1, OutOfRange::Error),
            Err(SpectrumError::OutOfRange { .. })
        ));
        assert!(matches!(
            interp_step_right(&xs, &ys, 1.5, OutOfRange::Error),
            Err(SpectrumError::OutOfRange { .. })
        ));
    }

    // -------- Cubic spline --------

    #[test]
    fn cubic_spline_reproduces_linear_function() {
        // Trapezoidal cubic spline of a strictly linear function should match
        // the analytic line within roundoff: y''(x) = 0 everywhere, so the
        // natural BCs are exact for any segment count.
        let xs: Vec<f64> = (0..6).map(|i| i as f64).collect();
        let ys: Vec<f64> = xs.iter().map(|x| 2.0 * x + 1.0).collect();
        let c = CubicSplineCoeffs::natural(&xs, &ys).unwrap();
        for x in [0.0, 0.25, 1.5, 3.7, 5.0] {
            let v = interp_cubic_spline(&xs, &ys, &c, x, OutOfRange::ClampToEndpoints).unwrap();
            assert_abs_diff_eq!(v, 2.0 * x + 1.0, epsilon = 1e-10);
        }
    }

    #[test]
    fn cubic_spline_passes_through_samples() {
        let xs = vec![0.0_f64, 1.0, 2.5, 4.0, 6.0];
        let ys = vec![1.0_f64, 0.0, -2.0, 0.5, 3.0];
        let c = CubicSplineCoeffs::natural(&xs, &ys).unwrap();
        for (i, x) in xs.iter().enumerate() {
            let v = interp_cubic_spline(&xs, &ys, &c, *x, OutOfRange::ClampToEndpoints).unwrap();
            assert_abs_diff_eq!(v, ys[i], epsilon = 1e-10);
        }
    }

    #[test]
    fn cubic_spline_natural_endpoints_have_zero_second_derivative() {
        let xs = vec![0.0_f64, 1.0, 2.5, 4.0, 6.0];
        let ys = vec![1.0_f64, 0.0, -2.0, 0.5, 3.0];
        let c = CubicSplineCoeffs::natural(&xs, &ys).unwrap();
        assert_eq!(c.y2[0], 0.0);
        assert_eq!(*c.y2.last().unwrap(), 0.0);
    }

    #[test]
    fn cubic_spline_two_samples_is_linear() {
        let xs = vec![0.0_f64, 4.0];
        let ys = vec![1.0_f64, 9.0];
        let c = CubicSplineCoeffs::natural(&xs, &ys).unwrap();
        for x in [0.0, 1.0, 2.0, 3.0, 4.0] {
            let v = interp_cubic_spline(&xs, &ys, &c, x, OutOfRange::ClampToEndpoints).unwrap();
            assert_abs_diff_eq!(v, 1.0 + 2.0 * x, epsilon = 1e-12);
        }
    }

    #[test]
    fn dispatch_routes_each_mode() {
        let (xs, ys) = linear();
        for mode in [
            Interpolation::Linear,
            Interpolation::Nearest,
            Interpolation::PiecewiseConstantLeft,
            Interpolation::PiecewiseConstantRight,
            Interpolation::CubicSpline,
        ] {
            let v = interp(&xs, &ys, 2.0, mode, OutOfRange::ClampToEndpoints).unwrap();
            // Every mode must reproduce the value exactly at a sample point.
            assert_abs_diff_eq!(v, 5.0, epsilon = 1e-10);
        }
    }
}
