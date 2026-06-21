//! Error types for starlight map construction and validation.

use crate::healpix::HealpixError;

/// Result alias for starlight map operations.
pub type Result<T> = std::result::Result<T, StellarMapError>;

/// Error type for stellar surface-brightness map construction and validation.
#[derive(Debug, thiserror::Error)]
pub enum StellarMapError {
    /// Wrapped HEALPix grid, indexing, or map-storage error.
    #[error(transparent)]
    Healpix(#[from] HealpixError),
    /// A magnitude value was not finite.
    #[error("apparent magnitude must be finite, got {0}")]
    InvalidMagnitude(f64),
    /// A catalogue weight was not finite or was negative.
    #[error("catalogue weight must be finite and nonnegative, got {0}")]
    InvalidWeight(f64),
    /// The radiance scale was not finite or was negative.
    #[error("integrated_per_v_s10 must be finite and nonnegative, got {0}")]
    InvalidRadianceScale(f64),
    /// No catalogue records survived the configured magnitude cuts.
    #[error("no stellar catalogue records survived filtering")]
    EmptyFilteredCatalogue,
    /// Records survived filtering but none contributed usable B- or V-band flux.
    #[error("no stellar catalogue records contributed usable B or V photometric flux")]
    NoUsablePhotometry,
    /// A map validation check failed.
    #[error("stellar map validation failed: {0}")]
    Validation(&'static str),
}
