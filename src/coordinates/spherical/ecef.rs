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
};
use crate::units::*;
use crate::bodies::EARTH;

impl<C: ReferenceCenter,> Direction<C, ECEF> {
    pub const fn new_const(lon: Degrees, lat: Degrees) -> Self {
        Self::new_raw(
            lat,
            lon,
            Quantity::<Unitless>::new(1.0)
        )
    }

    /// Creates a new geographic coordinate with normalized latitude and longitude.
    ///
    /// # Arguments
    /// - `lat`: Latitude in degrees, will be normalized to [-90°, 90°].
    /// - `lon`: Longitude in degrees, will be normalized to [-180°, 180°].
    pub fn new<T>(lon: T, lat: T) -> Self
    where
        T: Into<Degrees>,
    {
        Self::new_const(
            lat.into().wrap_quarter_fold(),
            lon.into().wrap_signed_lo()
        )
    }
}

impl<C: ReferenceCenter> Position<C, ECEF, Kilometer> {
    pub const fn new_const(lon: Degrees, lat: Degrees, alt: Kilometers) -> Self {
        Self::new_raw(
            lat,
            lon,
            Kilometers::new(EARTH.radius.value() + alt.value())
        )
    }

    /// Creates a new geographic coordinate with normalized latitude and longitude.
    ///
    /// # Arguments
    /// - `lat`: Latitude in degrees, will be normalized to [-90°, 90°].
    /// - `lon`: Longitude in degrees, will be normalized to [-180°, 180°].
    /// - `alt`: Altitude above the sea, in KM.
    pub fn new<A, T>(lon: A, lat: A, alt: T) -> Self
    where
        T: Into<Kilometers>,
        A: Into<Degrees>,
    {
        Self::new_const(
            lat.into().wrap_quarter_fold(),
            lon.into().wrap_signed_lo(),
            alt.into())
    }
}

impl<C: ReferenceCenter, U: Unit> SphericalCoord<C, ECEF, U> {
    /// Returns the latitude (φ) in degrees.
    pub fn lat(&self) -> Degrees { self.polar }

    /// Returns the longitude (λ) in degrees.
    pub fn lon(&self) -> Degrees { self.azimuth }
}


