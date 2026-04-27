// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Shared interpolation policies.
//!
//! This module defines policy enums that are reused across the
//! [`crate::spectra`] and [`crate::tables`] feature modules so callers can
//! mix typed sampled spectra and typed gridded tables without adapter glue.
//!
//! Currently provides:
//!
//! - [`OutOfRange`] — what to do when an evaluation point falls outside the
//!   sampled domain.

/// How to evaluate `y(x)` (or `f(x, y)`) outside the sampled domain.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum OutOfRange {
    /// Hold the nearest endpoint value (matches `numpy.interp` default).
    #[default]
    ClampToEndpoints,
    /// Return zero outside the domain.
    Zero,
    /// Surface an out-of-range error to the caller.
    Error,
}
