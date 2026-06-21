//! HEALPix grids, pixel indices, and frame-typed map containers.
//!
//! This module implements reusable library primitives for the Hierarchical Equal
//! Area isoLatitude Pixelisation (HEALPix) of the unit sphere. Public
//! pixelization APIs operate on typed [`crate::coordinates::cartesian::Direction`]
//! values so the same grid abstraction can be used with Galactic, equatorial,
//! ecliptic, or other reference frames without exposing frame-specific angular
//! names such as right ascension, declination, longitude, or latitude.
//!
//! Internally the RING kernel uses the standard spherical HEALPix variables
//! `theta`, `phi`, and `z = cos(theta)`.
//!
//! # References
//!
//! - Gorski, K. M. et al. (2005), "HEALPix: A Framework for High-Resolution
//!   Discretization and Fast Analysis of Data Distributed on the Sphere",
//!   Astrophysical Journal, 622, 759.
//! - Calabretta, M. R. and Roukema, B. F. (2007), "Mapping on the HEALPix grid",
//!   Monthly Notices of the Royal Astronomical Society, 381, 865.

mod error;
mod grid;
mod index;
mod map;
mod nside;
mod ordering;
mod ring;
mod validation;

#[cfg(test)]
mod tests;

pub use error::{HealpixError, Result};
pub use grid::HealpixGrid;
pub use index::HealpixIndex;
pub use map::HealpixMap;
pub use nside::{Nside, MAX_NSIDE};
pub use ordering::HealpixOrdering;
pub(crate) use validation::validate_healpix_map_complete;
#[cfg(test)]
pub(crate) use ring::direction_from_theta_phi;
