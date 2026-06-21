// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Integrated stellar surface-brightness map construction.
//!
//! This module builds Galactic HEALPix maps from stellar catalogue records. The
//! v1 model bins catalogue-star flux into B- and V-band tenth-magnitude-star
//! (`S10`) surface-brightness diagnostics, then scales V-band S10 to an
//! integrated photon radiance using a caller-provided calibration factor.
//!
//! The implementation is intended for reproducible data-asset generation. It is
//! not a full spectral synthesis or passband integration model.
//!
//! # References
//!
//! - Leinert, C. et al. (1998), "The 1997 reference of diffuse night sky
//!   brightness", Astronomy and Astrophysics Supplement Series, 127, 1.
//! - Gorski, K. M. et al. (2005), "HEALPix: A Framework for High-Resolution
//!   Discretization and Fast Analysis of Data Distributed on the Sphere",
//!   Astrophysical Journal, 622, 759.

mod brightness;
mod builder;
pub mod csv;
mod error;
mod magnitude;
mod map;
mod photometry;
mod provenance;
mod record;
mod validation;

#[cfg(test)]
mod tests;

pub use brightness::StellarSurfaceBrightness;
pub use builder::StellarSurfaceBrightnessMapBuilder;
pub use error::{Result, StellarMapError};
pub use magnitude::ApparentMagnitude;
pub use map::StellarSurfaceBrightnessMap;
pub use photometry::flux_10mag_units;
pub use provenance::StellarMapProvenance;
pub use record::StellarCatalogueRecord;
pub use validation::{
    validate_flux_conservation, validate_no_longitude_wrap_artifact,
    validate_plane_pole_contrast, validate_stellar_map_values,
};
