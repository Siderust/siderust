// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Shared interpolation policy types
//!
//! ## Scientific scope
//!
//! Tabulated astronomical functions — filter throughputs, atmospheric
//! transmission, sky-brightness models, ephemeris-derived quantities —
//! must specify what they do at *query points outside the tabulated
//! domain*. The conventional choices are: clamp to the nearest sampled
//! endpoint (the `numpy.interp` default and the most common assumption
//! in observational pipelines), return zero (appropriate for compactly
//! supported quantities like passband transmissions outside their
//! design range), or surface an explicit out-of-range error to the
//! caller (appropriate for safety-critical numerics where extrapolation
//! is never acceptable). Mixing these policies silently across modules
//! produces inconsistent answers.
//!
//! Centralising the policy enum here lets the `spectra` and `tables`
//! feature modules share a single typed
//! contract for out-of-range handling, so callers can mix sampled
//! spectra and gridded tables without writing adapter glue.
//!
//! ## Technical scope
//!
//! Defines [`OutOfRange`] with three variants:
//!
//! - [`OutOfRange::ClampToEndpoints`] — default; matches `numpy.interp`.
//! - [`OutOfRange::Zero`] — return `0` outside the domain.
//! - [`OutOfRange::Error`] — surface an out-of-range error from the
//!   typed evaluation API.
//!
//! No interpolation algorithm lives here; the kernels live in the
//! feature-gated `spectra::algo` and `tables::algo` modules.
//!
//! ## References
//!
//! - NumPy developers. *numpy.interp* documentation.
//!   <https://numpy.org/doc/stable/reference/generated/numpy.interp.html>
//!   (default behaviour: clamp to endpoints).

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
