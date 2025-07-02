//! Defines spherical coordinate types in the ecliptic frame, using latitude (B), longitude (L), and distance (R).
//!
//! These coordinates are commonly used in planetary and heliocentric astronomy for representing
//! object positions in the solar system relative to various centers (barycentric, heliocentric, etc.).
//!
//! # Coordinate Convention
//! The `Position<Center, Ecliptic>` type uses:
//!
//! - **Longitude (λ or L)**  → `azimuth`: angle from the ecliptic zero meridian, in degrees.
//! - **Latitude (β or B)**   → `polar`: angle from the ecliptic plane, in degrees.
//! - **Radial distance (R)** → distance from the reference center, in AU (astronomical units).
//!
//! Longitude is normalized to the [0°, 360°] range, and latitude to [-90°, 90°].
//!
//! # Example
//! ```rust
//! use siderust::coordinates::spherical::position::Ecliptic;
//! use siderust::units::Degrees;
//!
//! let coord = Ecliptic::new(
//!     Degrees::new(120.0), Degrees::new(5.0), 1.0
//! );
//! println!("lon = {}, lat = {}", coord.lon(), coord.lat());
//! ```

use super::*;
use crate::coordinates::{
    centers::*, frames,
    kinds::Kind,
};
use crate::units::{Degrees, Unit};

impl<C: ReferenceCenter, U: Unit> Position<C, frames::Ecliptic, U> {
    pub const fn new_const(lon: Degrees, lat: Degrees, distance: U) -> Self {
        Self::new_spherical_coord(lat, lon, Some(distance))
    }

    pub fn new(lon: Degrees, lat: Degrees, distance: U) -> Self {
        Self::new_const(
            lon.normalize(),
            lat.normalize_to_90_range(),
            distance)
    }
}

impl<C: ReferenceCenter> Direction<C, frames::Ecliptic> {
    /// Creates a new ecliptic direction with constant values.
    ///
    /// # Arguments
    /// - `lon`: Longitude (λ), in degrees.
    /// - `lat`: Latitude (β), in degrees.
    pub const fn new_const(lon: Degrees, lat: Degrees) -> Self {
        Self::new_spherical_coord(lat, lon, None)
    }

    /// Constructs a new ecliptic direction with normalized input angular.
    ///
    /// Longitude is normalized to [0°, 360°], latitude to [-90°, 90°].
    ///
    /// # Arguments
    /// - `lon`: Longitude (λ), in degrees.
    /// - `lat`: Latitude (β), in degrees.
    pub fn new(lon: Degrees, lat: Degrees) -> Self {
        Self::new_const(
            lon.normalize(),
            lat.normalize_to_90_range()
        )
    }
}

impl<C: ReferenceCenter, U: Unit, K: Kind> SphericalCoord<C, frames::Ecliptic, U, K> {
    /// Returns the Latitude (β) in degrees.
    pub fn lat(&self) -> Degrees { self.polar }

    /// Returns the Longitude (λ) in degrees.
    pub fn lon(&self) -> Degrees { self.azimuth }
}
