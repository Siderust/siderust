// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Errors raised by [`crate::spectra`].

use core::fmt;

/// Reasons a [`crate::spectra::SampledSpectrum`] operation can fail.
#[derive(Debug, Clone, PartialEq)]
pub enum SpectrumError {
    /// `xs` and `ys` had different lengths.
    LengthMismatch { xs: usize, ys: usize },
    /// Fewer than two samples were provided.
    TooFewSamples(usize),
    /// `xs` was not strictly monotonically increasing.
    NotMonotonic { index: usize },
    /// A query was outside the spectrum's domain and the configured
    /// [`crate::spectra::OutOfRange`] policy was [`crate::spectra::OutOfRange::Error`].
    OutOfRange { x: f64, lo: f64, hi: f64 },
    /// A loader failed to parse an input source.
    Parse(String),
}

impl fmt::Display for SpectrumError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SpectrumError::LengthMismatch { xs, ys } => {
                write!(f, "spectrum length mismatch: xs.len()={xs}, ys.len()={ys}")
            }
            SpectrumError::TooFewSamples(n) => {
                write!(f, "spectrum needs at least 2 samples (got {n})")
            }
            SpectrumError::NotMonotonic { index } => {
                write!(f, "spectrum xs not strictly increasing at index {index}")
            }
            SpectrumError::OutOfRange { x, lo, hi } => {
                write!(f, "spectrum query x={x} outside domain [{lo}, {hi}]")
            }
            SpectrumError::Parse(msg) => write!(f, "spectrum parse error: {msg}"),
        }
    }
}

impl std::error::Error for SpectrumError {}
