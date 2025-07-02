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
//! # Provided Types
//! The following type aliases define common combinations of center and ecliptic frame:
//!
//! - `EclipticBarycentricPos` → Center: solar system barycenter.
//! - `EclipticHeliocentricPos` → Center: Sun (true heliocentric).
//! - `EclipticGeocentricPos` → Center: Earth (common for Moon or planets).
//! - `EclipticTopocentricPos` → Center: specific observer on Earth.
//!
//! # Example
//! ```rust
//! use siderust::coordinates::spherical::EclipticPos;
//! use siderust::units::Degrees;
//!
//! let coord = EclipticPos::new(
//!     Degrees::new(120.0), Degrees::new(5.0), 1.0
//! );
//! println!("lon = {}, lat = {}", coord.lon(), coord.lat());
//! ```

use super::*;
use crate::coordinates::{
    frames::*,
    centers::*,
    kinds::Kind,
};
use crate::units::{Degrees, Unit};


// Ecliptic Coordinate Types
// Polar   -> Latitude  (B) – the angle from the equator. [-90°, 90°]
// Azimuth -> Longitude (L) – the angle from a prime meridian. [0°, 360°]
// Radial  -> Distance  (R) – the distance between the source and the target.
pub type EclipticPos<U> = Position<Heliocentric, Ecliptic, U>; // L (l), B (b), R (r)
pub type EclipticDir = Direction<Heliocentric, Ecliptic>; // L (l), B (b), R (r)

impl<C: ReferenceCenter, U: Unit> Position<C, Ecliptic, U> {
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

// Direction type aliases for Ecliptic frame
pub type EclipticBarycentricDir  = Direction<Barycentric,  Ecliptic>;
pub type EclipticHeliocentricDir = Direction<Heliocentric, Ecliptic>;
pub type EclipticGeocentricDir   = Direction<Geocentric,   Ecliptic>;
pub type EclipticTopocentricDir  = Direction<Topocentric,  Ecliptic>;

impl<C: ReferenceCenter> Direction<C, Ecliptic> {
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

impl<C: ReferenceCenter, U: Unit, K: Kind> SphericalCoord<C, Ecliptic, U, K> {
    /// Returns the Latitude (β) in degrees.
    pub fn lat(&self) -> Degrees { self.polar }

    /// Returns the Longitude (λ) in degrees.
    pub fn lon(&self) -> Degrees { self.azimuth }
}
