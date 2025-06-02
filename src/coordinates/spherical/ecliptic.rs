//! Defines spherical coordinate types in the ecliptic frame, using latitude (B), longitude (L), and distance (R).
//!
//! These coordinates are commonly used in planetary and heliocentric astronomy for representing
//! object positions in the solar system relative to various centers (barycentric, heliocentric, etc.).
//!
//! # Coordinate Convention
//! The `SphericalCoord<Center, Ecliptic>` type uses:
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
//! - `EclipticBarycentricSphericalCoord` → Center: solar system barycenter.
//! - `EclipticHeliocentricSphericalCoord` → Center: Sun (true heliocentric).
//! - `EclipticGeocentricSphericalCoord` → Center: Earth (common for Moon or planets).
//! - `EclipticTopocentricSphericalCoord` → Center: specific observer on Earth.
//!
//! # Example
//! ```rust
//! use siderust::coordinates::EclipticHeliocentricSphericalCoord;
//! use siderust::units::Degrees;
//!
//! let coord = EclipticHeliocentricSphericalCoord::new(
//!     Degrees::new(120.0), Degrees::new(5.0), 1.0
//! );
//! println!("lon = {}, lat = {}", coord.lon(), coord.lat());
//! ```

use super::SphericalCoord;
use crate::coordinates::{
    frames::*,
    centers::*,
};
use crate::units::Degrees;


// Ecliptic Coordinate Types
// Polar   -> Latitude  (B) – the angle from the equator. [-90°, 90°]
// Azimuth -> Longitude (L) – the angle from a prime meridian. [0°, 360°]
// Radial  -> Distance  (R) – the distance between the source and the target.
pub type EclipticBarycentricSphericalCoord  = SphericalCoord<Barycentric,  Ecliptic>;
pub type EclipticHeliocentricSphericalCoord = SphericalCoord<Heliocentric, Ecliptic>; // L (l), B (b), R (r)
pub type EclipticGeocentricSphericalCoord   = SphericalCoord<Geocentric,   Ecliptic>; // L (λ), B (β), R (Δ)
pub type EclipticTopocentricSphericalCoord  = SphericalCoord<Topocentric,  Ecliptic>;

impl<Center: ReferenceCenter> SphericalCoord<Center, Ecliptic> {
    pub const fn new_const(lon: Degrees, lat: Degrees, distance: f64) -> Self {
        SphericalCoord::new_spherical_coord(lat, lon, distance)
    }

    pub fn new(lon: Degrees, lat: Degrees, distance: f64) -> Self {
        SphericalCoord::<Center, Ecliptic>::new_const(
            lon.normalize(),
            lat.normalize_to_90_range(),
            distance)
    }

    pub fn lat(&self) -> Degrees { self.polar }
    pub fn lon(&self)  -> Degrees { self.azimuth }
}
