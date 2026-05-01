//! Generic typed sampled spectra — interpolation, integration, passbands.
//!
//! Re-exported from [`siderust`] as `siderust::spectra`. Downstream crates
//! that only need the spectra layer can depend on this crate directly.

pub mod spectra;
pub use spectra::*;

