use crate::healpix::HealpixError;

pub type Result<T> = std::result::Result<T, StellarMapError>;

#[derive(Debug, thiserror::Error)]
pub enum StellarMapError {
    #[error(transparent)]
    Healpix(#[from] HealpixError),
    #[error("apparent magnitude must be finite, got {0}")]
    InvalidMagnitude(f64),
    #[error("catalogue weight must be finite and nonnegative, got {0}")]
    InvalidWeight(f64),
    #[error("integrated_per_v_s10 must be finite and nonnegative, got {0}")]
    InvalidRadianceScale(f64),
    #[error("no stellar catalogue records survived filtering")]
    EmptyFilteredCatalogue,
    #[error("no stellar catalogue records contributed usable B or V photometric flux")]
    NoUsablePhotometry,
    #[error("stellar map validation failed: {0}")]
    Validation(&'static str),
}
