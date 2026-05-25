// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Interpolation policy types for sampled spectra
//!
//! ## Scientific scope
//!
//! Selects how a [`crate::spectra::SampledSpectrum`] evaluates `y(x)`
//! between adjacent samples. Synthetic photometry and atmospheric
//! transmission codes converge on piecewise-linear interpolation as the
//! conservative default — it preserves monotonic regions, has no
//! ringing, and matches the reference behaviour of `numpy.interp` and
//! `scipy.interpolate.interp1d(kind="linear")`. Nearest-neighbour and
//! piecewise-constant (left/right) variants are appropriate for
//! categorical or tabulated step data, and natural cubic splines provide
//! C² smoothness when the sampled function is known to be smooth (e.g.
//! laboratory-calibrated detector responses).
//!
//! ## Technical scope
//!
//! - [`Interpolation`] enum with variants `Linear`, `Nearest`,
//!   `PiecewiseConstantLeft`, `PiecewiseConstantRight`, and
//!   `CubicSpline` — all implemented.
//! - [`OutOfRange`] re-exported from [`crate::interp`] so the two
//!   policies can be referenced from a single import path.
//!
//! ## Tradeoffs
//!
//! | Mode                     | Continuity | Best for                                                       |
//! |--------------------------|------------|----------------------------------------------------------------|
//! | `Linear`                 | C⁰         | Default; conservative, no ringing.                             |
//! | `Nearest`                | step       | Categorical / quantised data; ties resolve to the lower index. |
//! | `PiecewiseConstantLeft`  | step       | Left-continuous histograms or holding-period quantities.       |
//! | `PiecewiseConstantRight` | step       | Right-continuous histograms or sample-on-demand semantics.     |
//! | `CubicSpline`            | C²         | Smoothly tabulated functions; may overshoot near sharp edges.  |
//!
//! ## References
//!
//! - NumPy developers. *numpy.interp* documentation
//!   (linear interpolation reference).
//! - SciPy developers. *scipy.interpolate.interp1d* documentation.
//! - Press, W. H. et al. (1992). *Numerical Recipes in C*, 2nd ed., §3.3
//!   (natural cubic spline).

pub use crate::interp::OutOfRange;

/// How to evaluate `y(x)` between samples.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Interpolation {
    /// Piecewise linear interpolation between adjacent samples. Matches
    /// `numpy.interp` and `scipy.interpolate.interp1d(kind="linear")`.
    #[default]
    Linear,
    /// Nearest-neighbour interpolation. Ties resolve to the lower index,
    /// matching `scipy.interpolate.interp1d(kind="nearest")`.
    Nearest,
    /// Left-continuous step interpolation: `y(x) = ys[i-1]` for
    /// `x in [xs[i-1], xs[i])`. Equivalent to
    /// `scipy.interpolate.interp1d(kind="previous")`.
    PiecewiseConstantLeft,
    /// Right-continuous step interpolation: `y(x) = ys[i]` for
    /// `x in (xs[i-1], xs[i]]`. Equivalent to
    /// `scipy.interpolate.interp1d(kind="next")`.
    PiecewiseConstantRight,
    /// Natural cubic spline (`y''` at both endpoints clamped to zero).
    /// Coefficients are precomputed once at construction time when
    /// installed on a [`crate::spectra::SampledSpectrum`].
    CubicSpline,
}
