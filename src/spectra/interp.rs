// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Interpolation policies for [`crate::spectra::SampledSpectrum`].
//!
//! [`OutOfRange`] lives in [`crate::interp`] so it can be shared with
//! [`crate::tables`]; it is re-exported here for backwards compatibility.

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
