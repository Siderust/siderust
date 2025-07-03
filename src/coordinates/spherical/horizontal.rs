//! Defines spherical coordinate types in the horizontal frame, using altitude (Alt), azimuth (Az), and distance (d).
//!
//! These coordinates are commonly used in observational astronomy to describe the position of celestial objects
//! relative to an observer's local horizon.
//!
//! # Coordinate Convention
//! The `Position<Center, Horizontal>` type uses:
//!
//! - **Altitude (Alt)** → `polar`: angle above or below the horizon, in degrees.
//! - **Azimuth (Az)**   → `azimuth`: angle measured clockwise from the north, in degrees.
//! - **Radial distance (d)** → distance from the reference center, typically in astronomical units (AU).
//!
//! Altitude is normalized to the [-90°, 90°] range, and azimuth to the [0°, 360°] range.
//! Time is represented as Julian Date (JD).
//!
//! # Example
//! ```rust
//! use siderust::coordinates::spherical::position::Horizontal;
//! use siderust::units::Degrees;
//!
//! let coord = Horizontal::new(
//!     Degrees::new(45.0), Degrees::new(120.0), 1.0
//! );
//! println!("alt = {}, az = {}", coord.alt(), coord.az());
//! ```

use super::*;
use crate::coordinates::{
    frames::*,
    centers::*,
    kinds::Kind,
};
use crate::units::{Degrees, Distance};

impl<C: ReferenceCenter> Direction<C, Horizontal> {
    /// Creates a new horizontal spherical coordinate with constant values.
    ///
    /// # Arguments
    /// - `alt`: Altitude (α), in degrees.
    /// - `az`: Azimuth (θ), in degrees.
    ///
    /// # Returns
    /// A new `Position` in the horizontal frame.
    pub const fn new_const(alt: Degrees, az: Degrees) -> Self {
        Self::new_spherical_coord(alt, az, None)
    }

    /// Constructs a new horizontal spherical coordinate with normalized input angular.
    ///
    /// Altitude is normalized to the [-90°, 90°] range, and azimuth to the [0°, 360°] range.
    ///
    /// # Arguments
    /// - `alt`: Altitude (α), in degrees.
    /// - `az`: Azimuth (θ), in degrees.
    ///
    /// # Returns
    /// A `Direction` in the horizontal frame.
    pub fn new(alt: Degrees, az: Degrees) -> Self {
        Self::new_const(
            alt.normalize_to_90_range(),
            az.normalize()
        )
    }
}

impl<C: ReferenceCenter, U: Distance> Position<C, Horizontal, U> {
    /// Creates a new horizontal spherical coordinate with constant values.
    ///
    /// # Arguments
    /// - `alt`: Altitude (α), in degrees.
    /// - `az`: Azimuth (θ), in degrees.
    /// - `distance`: Distance to the object, typically in astronomical units (AU).
    ///
    /// # Returns
    /// A new `Position` in the horizontal frame.
    pub const fn new_const(alt: Degrees, az: Degrees, distance: U) -> Self {
        Self::new_spherical_coord(alt, az, Some(distance))
    }

    /// Constructs a new horizontal spherical coordinate with normalized input angular.
    ///
    /// Altitude is normalized to the [-90°, 90°] range, and azimuth to the [0°, 360°] range.
    ///
    /// # Arguments
    /// - `alt`: Altitude (α), in degrees.
    /// - `az`: Azimuth (θ), in degrees.
    /// - `distance`: Distance to the object, typically in astronomical units (AU).
    ///
    /// # Returns
    /// A `Position` in the horizontal frame.
    pub fn new(alt: Degrees, az: Degrees, distance: U) -> Self {
        Self::new_const(
            alt.normalize_to_90_range(),
            az.normalize(),
            distance)
    }
}

impl<C: ReferenceCenter, U: Distance, K: Kind> SphericalCoord<C, Horizontal, U, K> {
    /// Returns the Altitude (α) in degrees.
    pub fn alt(&self) -> Degrees { self.polar }

    /// Returns the Azimuth (θ) in degrees.
    pub fn az(&self) -> Degrees { self.azimuth }
}
