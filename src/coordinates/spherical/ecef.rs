//! Defines geographic coordinates using spherical representation
//! with respect to the Earth-Centered Earth-Fixed (ECEF) reference frame.
//!
//! This module provides:
//! - The `GeographicCoord` type, a spherical coordinate
//!   system using latitude, longitude, and radial height.
//! - Creation methods that normalize angular values to valid geodetic ranges.
//! - Traits and builders for structured instantiation.
//!
//! # Coordinate Meaning
//! The spherical coordinates follow geodetic conventions:
//!
//! - **Latitude (φ)**   → `polar`: angle from the equator in degrees, normalized to [-90°, 90°].
//! - **Longitude (λ)**  → `azimuth`: angle from the prime meridian in degrees, normalized to [-180°, 180°].
//! - **Radial distance (h)** → height above the reference ellipsoid (e.g., WGS84), in meters.
//!
//! # Type Aliases
//! - `GeographicCoord` = `SphericalCoord<Geocentric, ECEF>`
//!
//! # Example
//! ```rust
//! use siderust::coordinates::GeographicCoord;
//! use siderust::units::Degrees;
//!
//! let coord = GeographicCoord::new(Degrees::new(45.0), Degrees::new(7.0), 400.0);
//! println!("lat = {}, lon = {}", coord.lat(), coord.lon());
//! ```

use super::SphericalCoord;
use crate::coordinates::{
    frames::*,
    centers::*,
};
use crate::units::Degrees;

impl<Center: ReferenceCenter> SphericalCoord<Center, ECEF> {
    pub const fn new_const(lon: Degrees, lat: Degrees, distance: f64) -> Self {
        SphericalCoord::new_spherical_coord(lat, lon, distance)
    }

    /// Creates a new geographic coordinate with normalized latitude and longitude.
    ///
    /// # Arguments
    /// - `lat`: Latitude in degrees, will be normalized to [-90°, 90°].
    /// - `lon`: Longitude in degrees, will be normalized to [-180°, 180°].
    /// - `distance`: Height or distance above the ellipsoid, in meters.
    pub fn new(lon: Degrees, lat: Degrees, distance: f64) -> Self {
        Self::new_const(
            lat.normalize_to_90_range(),
            lon.normalize_to_180_range(),
            distance)
    }

    /// Returns the latitude (φ) in degrees.
    pub fn lat(&self) -> Degrees { self.polar }

    /// Returns the longitude (λ) in degrees.
    pub fn lon(&self) -> Degrees { self.azimuth }
}


// Define type alias for Geodetic/Geographic coordinates using SphericalCoord
// Polar   -> Latitude (φ) – the angle from the equator. [-90°, 90°]
// Azimuth -> Longitude (λ) – the angle from a prime meridian. [-180°, 180°]
// Radial  -> Height (h) – the elevation above the reference ellipsoid (such as WGS84).
pub type GeographicCoord = SphericalCoord<Geocentric, ECEF>;
