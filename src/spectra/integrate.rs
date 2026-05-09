// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Free-function integration helpers
//!
//! ## Scientific scope
//!
//! Re-exports the trapezoidal integration kernels from
//! [`crate::spectra::algo`] under shorter free-function names for
//! callers that prefer top-level functions over the inherent methods on
//! [`crate::spectra::SampledSpectrum`]. The underlying numerical method
//! (composite trapezoidal rule on a strictly-monotonic axis) is the same
//! used by every reference synthetic-photometry pipeline.
//!
//! ## Technical scope
//!
//! - [`integrate`] = [`crate::spectra::algo::trapz`].
//! - [`integrate_range`] = [`crate::spectra::algo::trapz_range`].
//! - [`integrate_weighted`] = [`crate::spectra::algo::trapz_weighted`].
//!
//! All operate on raw `&[f64]` slices.
//!
//! ## References
//!
//! - See [`crate::spectra::algo`] for the underlying trapezoidal-rule
//!   citations (Atkinson 1989; Press et al. 1992).

pub use super::algo::{
    trapz as integrate, trapz_range as integrate_range, trapz_weighted as integrate_weighted,
};
