// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Errors surfaced by the gridded-table API.

use core::fmt;

/// Errors returned by [`Grid1D`](super::Grid1D), [`Grid2D`](super::Grid2D),
/// and the [`algo`](super::algo) kernels.
#[derive(Debug, Clone, PartialEq)]
pub enum TableError {
    /// Axis and table dimensions disagree.
    ShapeMismatch {
        expected_x: usize,
        expected_y: usize,
        actual_rows: usize,
        actual_cols: usize,
    },
    /// One of the axis arrays is shorter than two samples.
    TooFewSamples { axis: &'static str, len: usize },
    /// An axis is not strictly increasing.
    NotMonotonic { axis: &'static str, at_index: usize },
    /// Evaluation point falls outside the sampled domain and the
    /// active [`OutOfRange::Error`](super::OutOfRange::Error) policy
    /// surfaces the failure.
    OutOfRange {
        axis: &'static str,
        value: f64,
        lo: f64,
        hi: f64,
    },
}

impl fmt::Display for TableError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            TableError::ShapeMismatch {
                expected_x,
                expected_y,
                actual_rows,
                actual_cols,
            } => write!(
                f,
                "table shape mismatch: expected {expected_y}×{expected_x}, got {actual_rows}×{actual_cols}"
            ),
            TableError::TooFewSamples { axis, len } => {
                write!(f, "{axis} axis has too few samples ({len}); need ≥ 2")
            }
            TableError::NotMonotonic { axis, at_index } => write!(
                f,
                "{axis} axis is not strictly increasing at index {at_index}"
            ),
            TableError::OutOfRange { axis, value, lo, hi } => write!(
                f,
                "{axis} value {value} is outside the sampled domain [{lo}, {hi}]"
            ),
        }
    }
}

impl std::error::Error for TableError {}
