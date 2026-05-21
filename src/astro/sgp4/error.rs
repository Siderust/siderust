// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Unified error type for the SGP4/SDP4 propagator.

use thiserror::Error;

/// Errors produced while constructing or evaluating an [`Sgp4Propagator`](crate::sgp4::Sgp4Propagator).
///
/// The variants split cleanly into three groups:
///
/// * **Initialisation** ([`InvalidElements`](Sgp4Error::InvalidElements),
///   [`InvalidEpoch`](Sgp4Error::InvalidEpoch)) — the TLE record is structurally
///   well-formed but its mean elements cannot be turned into SGP4 secular /
///   periodic constants (e.g. eccentricity outside `[0, 1)`, sub-orbital mean
///   motion, malformed UTC epoch).
/// * **Propagation** ([`Propagation`](Sgp4Error::Propagation)) — the underlying
///   SGP4/SDP4 polynomial diverged at the requested epoch (typical for very
///   large `|Δt|` on decayed objects, or when Lyddane choice fails on a
///   geostationary regime — see Vallado 2006).
/// * **Time conversion** ([`TimeConversion`](Sgp4Error::TimeConversion)) — the
///   target epoch supplied as a [`tempoch::JulianDate<tempoch::UTC>`] cannot be
///   represented as a calendar instant, e.g. it falls outside the leap-second
///   table covered by the active [`tempoch::TimeContext`].
///
/// # Examples
///
/// ```
/// use siderust::astro::sgp4::Sgp4Error;
/// let e = Sgp4Error::InvalidEpoch("year out of range");
/// assert!(matches!(e, Sgp4Error::InvalidEpoch(_)));
/// assert!(e.to_string().contains("epoch"));
/// ```
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum Sgp4Error {
    /// The TLE mean elements failed SGP4 initialisation. The `details` payload
    /// is the human-readable diagnostic produced by the underlying propagator.
    #[error("invalid TLE elements for SGP4: {details}")]
    InvalidElements {
        /// Free-form diagnostic from the SGP4 initialiser.
        details: String,
    },

    /// The TLE epoch could not be expressed as a calendar instant.
    #[error("invalid TLE epoch: {0}")]
    InvalidEpoch(&'static str),

    /// The SGP4/SDP4 polynomial returned an error at the requested epoch.
    #[error("SGP4 propagation failed: {details}")]
    Propagation {
        /// Free-form diagnostic from the SGP4 evaluator.
        details: String,
    },

    /// Conversion between [`tempoch::JulianDate`] and a calendar instant failed.
    #[error("UTC time conversion failed: {0}")]
    TimeConversion(String),
}
