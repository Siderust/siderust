// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Free-function integration helpers.
//!
//! These re-export the underlying [`crate::spectra::algo`] primitives under
//! shorter names for callers that prefer free functions over the
//! [`crate::spectra::SampledSpectrum`] inherent methods.

pub use super::algo::{trapz as integrate, trapz_range as integrate_range, trapz_weighted as integrate_weighted};
