//! Defines geographic coordinates using spherical representation
//! with respect to the Earth-Centered Earth-Fixed (ECEF) reference frame.
//!
//! This module provides:
//! - The `Geographic` type, a spherical coordinate
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
//! # Example
//! ```rust
//! use siderust::coordinates::spherical::position::Geographic;
//! use siderust::units::{Degrees, KM};
//!
//! let coord = Geographic::new(Degrees::new(45.0), Degrees::new(7.0), 2.4*KM);
//! println!("lat = {}, lon = {}", coord.lat(), coord.lon());
//! ```

use super::*;
use crate::coordinates::{
    frames::*,
    centers::*,
    kinds::Kind,
};
use crate::units::{Degrees, Kilometers, Unit};
use crate::bodies::EARTH;

impl<C: ReferenceCenter,> Direction<C, ECEF> {
    pub const fn new_const(lon: Degrees, lat: Degrees) -> Self {
        Self::new_spherical_coord(
            lat,
            lon,
            None
        )
    }

    /// Creates a new geographic coordinate with normalized latitude and longitude.
    ///
    /// # Arguments
    /// - `lat`: Latitude in degrees, will be normalized to [-90°, 90°].
    /// - `lon`: Longitude in degrees, will be normalized to [-180°, 180°].
    pub fn new(lon: Degrees, lat: Degrees) -> Self {
        Self::new_const(
            lat.normalize_to_90_range(),
            lon.normalize_to_180_range()
        )
    }
}

impl<C: ReferenceCenter> Position<C, ECEF, Kilometers> {
    pub const fn new_const(lon: Degrees, lat: Degrees, alt: Kilometers) -> Self {
        Self::new_spherical_coord(
            lat,
            lon,
            Some(Kilometers::new(EARTH.radius.value() + alt.value()))
        )
    }

    /// Creates a new geographic coordinate with normalized latitude and longitude.
    ///
    /// # Arguments
    /// - `lat`: Latitude in degrees, will be normalized to [-90°, 90°].
    /// - `lon`: Longitude in degrees, will be normalized to [-180°, 180°].
    /// - `alt`: Altitude above the sea, in Kilometers.
    pub fn new(lon: Degrees, lat: Degrees, alt: Kilometers) -> Self {
        Self::new_const(
            lat.normalize_to_90_range(),
            lon.normalize_to_180_range(),
            alt)
    }
}

impl<C: ReferenceCenter, U: Unit, K: Kind> SphericalCoord<C, ECEF, U, K> {
    /// Returns the latitude (φ) in degrees.
    pub fn lat(&self) -> Degrees { self.polar }

    /// Returns the longitude (λ) in degrees.
    pub fn lon(&self) -> Degrees { self.azimuth }
}


