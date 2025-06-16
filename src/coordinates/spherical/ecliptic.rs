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
//! use siderust::coordinates::spherical::EclipticHeliocentricPos;
//! use siderust::units::Degrees;
//!
//! let coord = EclipticHeliocentricPos::new(
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
use crate::units::Degrees;


// Ecliptic Coordinate Types
// Polar   -> Latitude  (B) – the angle from the equator. [-90°, 90°]
// Azimuth -> Longitude (L) – the angle from a prime meridian. [0°, 360°]
// Radial  -> Distance  (R) – the distance between the source and the target.
pub type EclipticBarycentricPos  = Position<Barycentric,  Ecliptic>;
pub type EclipticHeliocentricPos = Position<Heliocentric, Ecliptic>; // L (l), B (b), R (r)
pub type EclipticGeocentricPos   = Position<Geocentric,   Ecliptic>; // L (λ), B (β), R (Δ)
pub type EclipticTopocentricPos  = Position<Topocentric,  Ecliptic>;

impl<Center: ReferenceCenter> Position<Center, Ecliptic> {
    pub const fn new_const(lon: Degrees, lat: Degrees, distance: f64) -> Self {
        Position::new_spherical_coord(lat, lon, Some(distance))
    }

    pub fn new(lon: Degrees, lat: Degrees, distance: f64) -> Self {
        Position::<Center, Ecliptic>::new_const(
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

impl<Center: ReferenceCenter> Direction<Center, Ecliptic> {
    /// Creates a new ecliptic direction with constant values.
    ///
    /// # Arguments
    /// - `lon`: Longitude (λ), in degrees.
    /// - `lat`: Latitude (β), in degrees.
    pub const fn new_const(lon: Degrees, lat: Degrees) -> Self {
        Direction::new_spherical_coord(lat, lon, None)
    }

    /// Constructs a new ecliptic direction with normalized input angular.
    ///
    /// Longitude is normalized to [0°, 360°], latitude to [-90°, 90°].
    ///
    /// # Arguments
    /// - `lon`: Longitude (λ), in degrees.
    /// - `lat`: Latitude (β), in degrees.
    pub fn new(lon: Degrees, lat: Degrees) -> Self {
        Direction::<Center, Ecliptic>::new_const(
            lon.normalize(),
            lat.normalize_to_90_range()
        )
    }
}

impl<C: ReferenceCenter, K: Kind> SphericalCoord<C, Ecliptic, K> {
    /// Returns the Latitude (β) in degrees.
    pub fn lat(&self) -> Degrees { self.polar }

    /// Returns the Longitude (λ) in degrees.
    pub fn lon(&self) -> Degrees { self.azimuth }
}
