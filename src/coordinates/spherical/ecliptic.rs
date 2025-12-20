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
//! - **Radial distance (R)** → distance from the reference center, in AstronomicalUnits (astronomical units).
//!
//! Longitude is normalized to the [0°, 360°] range, and latitude to [-90°, 90°].
//!
//! # Example
//! ```rust
//! use siderust::coordinates::spherical::direction::Ecliptic;
//! use qtty::*;
//!
//! let coord = Ecliptic::new(
//!     120.0 * DEG, 5.0 * DEG
//! );
//! println!("lon = {}, lat = {}", coord.lon(), coord.lat());
//! ```

use super::*;
use crate::coordinates::{centers::*, frames};
use qtty::*;

impl<C: ReferenceCenter<Params = ()>, U: LengthUnit> Position<C, frames::Ecliptic, U> {
    pub const fn new_const(lon: Degrees, lat: Degrees, distance: Quantity<U>) -> Self {
        Self::new_raw(lat, lon, distance)
    }

    pub fn new<A, T>(lon: A, lat: A, distance: T) -> Self
    where
        T: Into<Quantity<U>>,
        A: Into<Degrees>,
    {
        Self::new_const(
            lon.into().normalize(),
            lat.into().wrap_quarter_fold(),
            distance.into(),
        )
    }
}

impl direction::Direction<frames::Ecliptic> {
    /// Creates a new ecliptic direction with constant values.
    ///
    /// # Arguments
    /// - `lon`: Longitude (λ), in degrees.
    /// - `lat`: Latitude (β), in degrees.
    pub fn new_ecliptic(lon: Degrees, lat: Degrees) -> Self {
        Self::new(lat.wrap_quarter_fold(), lon.normalize())
    }

    /// Returns the Latitude (β) in degrees.
    pub fn lat(&self) -> Degrees {
        self.polar
    }

    /// Returns the Longitude (λ) in degrees.
    pub fn lon(&self) -> Degrees {
        self.azimuth
    }
}

impl<C: ReferenceCenter, U: Unit> Position<C, frames::Ecliptic, U> {
    /// Returns the Latitude (β) in degrees.
    pub fn lat(&self) -> Degrees {
        self.polar
    }

    /// Returns the Longitude (λ) in degrees.
    pub fn lon(&self) -> Degrees {
        self.azimuth
    }
}
