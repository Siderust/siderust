// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Interpolation and out-of-range policies for [`crate::spectra::SampledSpectrum`].

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

/// How to evaluate `y(x)` outside the sampled domain `[xs[0], xs[-1]]`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum OutOfRange {
    /// Hold the nearest endpoint value (matches `numpy.interp` default).
    #[default]
    ClampToEndpoints,
    /// Return zero outside the domain.
    Zero,
    /// Surface a [`crate::spectra::SpectrumError::OutOfRange`].
    Error,
}
