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
//! `scipy.interpolate.interp1d(kind="linear")`. The other variants in
//! the enum are reserved for future implementations and currently
//! reject at construction time.
//!
//! ## Technical scope
//!
//! - [`Interpolation`] enum with variants `Linear` (implemented),
//!   `Nearest`, `PiecewiseConstantLeft`, `PiecewiseConstantRight`,
//!   `CubicSpline` (reserved).
//! - [`OutOfRange`] re-exported from [`crate::interp`] so the two
//!   policies can be referenced from a single import path.
//!
//! ## References
//!
//! - NumPy developers. *numpy.interp* documentation
//!   (linear interpolation reference).
//! - SciPy developers. *scipy.interpolate.interp1d* documentation.

pub use crate::interp::OutOfRange;

/// How to evaluate `y(x)` between samples.
///
/// Only [`Interpolation::Linear`] is implemented today; the other variants are
/// reserved for future work and will be hard errors at construction time when
/// passed.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Interpolation {
    /// Piecewise linear interpolation between adjacent samples. This matches
    /// `numpy.interp` and `scipy.interpolate.interp1d(kind="linear")`.
    #[default]
    Linear,
    /// Reserved for a future nearest-neighbour implementation.
    Nearest,
    /// Reserved for a future left-continuous step interpolation.
    PiecewiseConstantLeft,
    /// Reserved for a future right-continuous step interpolation.
    PiecewiseConstantRight,
    /// Reserved for a future natural cubic spline implementation.
    CubicSpline,
}
