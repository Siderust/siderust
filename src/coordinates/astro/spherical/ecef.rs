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
//! use qtty::*;
//!
//! let coord = Geographic::new(Degrees::new(45.0), Degrees::new(7.0), 2.4*KM);
//! println!("lat = {}, lon = {}", coord.lat(), coord.lon());
//! ```

use super::*;
use crate::bodies::EARTH;
use crate::coordinates::algebra::{centers::*, frames::*};
use qtty::*;

impl direction::Direction<ECEF> {
    /// Creates a new geographic direction with normalized latitude and longitude.
    pub fn new_geographic(lon: Degrees, lat: Degrees) -> Self {
        Self::new(lat.wrap_quarter_fold(), lon.wrap_signed_lo())
    }

    /// Returns the latitude (φ) in degrees.
    pub fn lat(&self) -> Degrees {
        self.polar
    }

    /// Returns the longitude (λ) in degrees.
    pub fn lon(&self) -> Degrees {
        self.azimuth
    }
}

impl<C: ReferenceCenter<Params = ()>> Position<C, ECEF, Kilometer> {
    pub const fn new_const(lon: Degrees, lat: Degrees, alt: Kilometers) -> Self {
        Self::new_raw(
            lat,
            lon,
            Kilometers::new(EARTH.radius.value() + alt.value()),
        )
    }

    /// Creates a new geographic coordinate with normalized latitude and longitude.
    ///
    /// # Arguments
    /// - `lon`: Longitude in degrees, will be normalized to [-180°, 180°].
    /// - `lat`: Latitude in degrees, will be normalized to [-90°, 90°].
    /// - `alt`: Altitude above the sea, in KM.
    pub fn new<A, T>(lon: A, lat: A, alt: T) -> Self
    where
        T: Into<Kilometers>,
        A: Into<Degrees>,
    {
        Self::new_const(
            lon.into().wrap_signed_lo(),
            lat.into().wrap_quarter_fold(),
            alt.into(),
        )
    }
}

impl<C: ReferenceCenter, U: LengthUnit> Position<C, ECEF, U> {
    /// Returns the latitude (φ) in degrees.
    pub fn lat(&self) -> Degrees {
        self.polar
    }

    /// Returns the longitude (λ) in degrees.
    pub fn lon(&self) -> Degrees {
        self.azimuth
    }
}
