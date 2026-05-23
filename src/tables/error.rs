// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Error taxonomy for gridded tables
//!
//! ## Scientific scope
//!
//! Errors raised when constructing or evaluating a typed gridded table.
//! They cover invariants the typed containers enforce — axis/value
//! shape agreement, minimum sample count per axis, strict monotonicity
//! of each axis (either ascending or descending) — and out-of-domain
//! queries that the configured boundary policy turns into surfaced
//! errors. These are infrastructure errors, not physical failure modes.
//!
//! ## Technical scope
//!
//! [`TableError`] enumerates: `ShapeMismatch`, `TooFewSamples`,
//! `NotMonotonic`, `OutOfRange`. Implements [`core::fmt::Display`].
//!
//! ## References
//!
//! None — error infrastructure.

use core::fmt;

/// Errors returned by [`Grid1D`](super::Grid1D), [`Grid2D`](super::Grid2D),
/// and the [`algo`](super::algo) kernels.
#[derive(Debug, Clone, PartialEq)]
pub enum TableError {
    /// Axis and table dimensions disagree.
    ShapeMismatch {
        /// Expected number of samples along the x axis.
        expected_x: usize,
        /// Expected number of samples along the y axis.
        expected_y: usize,
        /// Actual number of rows in the provided value matrix.
        actual_rows: usize,
        /// Actual number of columns in the provided value matrix.
        actual_cols: usize,
    },
    /// One of the axis arrays is shorter than two samples.
    TooFewSamples {
        /// Name of the undersampled axis.
        axis: &'static str,
        /// Number of samples present on that axis.
        len: usize,
    },
    /// An axis is not strictly monotonic (equal consecutive values or a
    /// direction change after the first step).
    NotMonotonic {
        /// Name of the axis that violates monotonicity.
        axis: &'static str,
        /// Index of the first offending sample.
        at_index: usize,
    },
    /// Evaluation point falls outside the sampled domain and the
    /// active [`OutOfRange::Error`](super::OutOfRange::Error) policy
    /// surfaces the failure.
    OutOfRange {
        /// Name of the axis whose query is out of range.
        axis: &'static str,
        /// Queried coordinate value.
        value: f64,
        /// Lower bound of the sampled interval.
        lo: f64,
        /// Upper bound of the sampled interval.
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
                "{axis} axis is not strictly monotonic at index {at_index}"
            ),
            TableError::OutOfRange { axis, value, lo, hi } => write!(
                f,
                "{axis} value {value} is outside the sampled domain [{lo}, {hi}]"
            ),
        }
    }
}

impl std::error::Error for TableError {}
