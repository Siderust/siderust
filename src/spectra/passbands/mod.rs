// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Photometric filter passband datasets.
//!
//! Each sub-module provides lazily-initialised, statically-cached
//! [`SampledSpectrum`](crate::spectra::SampledSpectrum) constants for a
//! particular photometric system, following the same `OnceLock` pattern as
//! [`crate::atmosphere::ozone`].
//!
//! ## Units
//!
//! Filter throughputs are dimensionless values in \[0, 1\] represented by the
//! [`Throughput`] unit marker defined in this module.  It is analogous to
//! [`crate::atmosphere::ozone::Transmittance`] but scoped to filter passbands
//! so that the two concepts remain distinguishable in type signatures.

pub mod bessell1990;

use crate::ext_qtty::Unit;

/// Dimensionless unit marker for filter throughput values ∈ [0, 1].
///
/// `RATIO = 1.0`; canonical SI dimension is dimensionless.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Throughput;

impl Unit for Throughput {
    const RATIO: f64 = 1.0;
    type Dim = crate::ext_qtty::Dimensionless;
    const SYMBOL: &'static str = "";
}
