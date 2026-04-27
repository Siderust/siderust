// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Untyped numerical kernels for sampled spectra.
//!
//! These functions operate on plain `&[f64]` slices and are the reference
//! implementation that the typed [`crate::spectra::SampledSpectrum`] API
//! delegates to. They are exposed because legacy code paths (notably NSB) and
//! external loaders may already hold raw vectors; keeping a single
//! implementation guarantees bit-for-bit parity between the typed and untyped
//! surfaces.
//!
//! Pre-conditions (debug-asserted):
//! - `xs.len() == ys.len()`
//! - `xs.len() >= 2`
//! - `xs` is strictly monotonically increasing.

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

/// Dispatches an interpolation policy to its implementation.
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
        Interpolation::Nearest
        | Interpolation::PiecewiseConstantLeft
        | Interpolation::PiecewiseConstantRight
        | Interpolation::CubicSpline => Err(SpectrumError::Parse(format!(
            "interpolation {interp:?} is not implemented yet",
        ))),
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
        return Err(SpectrumError::LengthMismatch { xs: xs.len(), ys: ys.len() });
    }
    if xs.len() < 2 {
        return Err(SpectrumError::TooFewSamples(xs.len()));
    }
    for i in 1..xs.len() {
        if !(xs[i] > xs[i - 1]) {
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
}
