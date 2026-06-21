// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

use crate::healpix::HealpixError;

/// Result alias for stellar surface-brightness map operations.
pub type Result<T> = std::result::Result<T, StellarMapError>;

/// Error type for stellar surface-brightness map construction and validation.
#[derive(Debug, thiserror::Error)]
pub enum StellarMapError {
    /// Wrapped HEALPix error.
    #[error(transparent)]
    Healpix(#[from] HealpixError),
    /// A magnitude value was invalid.
    #[error("apparent magnitude must be finite, got {0}")]
    InvalidMagnitude(f64),
    /// A catalogue weight was invalid.
    #[error("catalogue weight must be finite and non-negative, got {0}")]
    InvalidWeight(f64),
    /// A radiance scale was invalid.
    #[error("integrated_per_v_s10 must be finite and non-negative, got {0}")]
    InvalidRadianceScale(f64),
    /// Smoothing was requested but is not implemented in the v1 builder.
    #[error("stellar-map smoothing is not implemented in the v1 builder")]
    UnsupportedSmoothing,
    /// No catalogue records survived filtering.
    #[error("no stellar catalogue records survived filtering")]
    EmptyFilteredCatalogue,
    /// No record contributed usable B or V photometric flux.
    #[error("no stellar catalogue records contributed usable B or V photometric flux")]
    NoUsablePhotometry,
    /// A validation check failed.
    #[error("stellar map validation failed: {0}")]
    Validation(&'static str),
}
