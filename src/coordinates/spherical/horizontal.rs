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
//! # Provided Types
//! The following type aliases define common combinations of center and horizontal frame:
//!
//! - `HorizontalBarycentricSphericalCoord` → Center: solar system barycenter.
//! - `HorizontalHeliocentricSphericalCoord` → Center: Sun.
//! - `HorizontalGeocentricSphericalCoord` → Center: Earth (common for ground-based observations).
//! - `HorizontalTopocentricSphericalCoord` → Center: specific observer on Earth.
//!
//! # Example
//! ```rust
//! use siderust::coordinates::spherical::HorizontalTopocentricSphericalCoord;
//! use siderust::units::Degrees;
//!
//! let coord = HorizontalTopocentricSphericalCoord::new(
//!     Degrees::new(45.0), Degrees::new(120.0), 1.0
//! );
//! println!("alt = {}, az = {}", coord.alt(), coord.az());
//! ```

use super::Position;
use crate::coordinates::{
    frames::*,
    centers::*,
};
use crate::units::Degrees;

// Horizontal Coordinate Types
// Polar   -> Alt (α) – the angle from the horizon. [-90°, 90°]
// Azimuth -> Az (θ) – the angle from a prime meridian. [0°, 360°]
// Radial  -> Distance (d) – the distance between the source and the target.
pub type HorizontalBarycentricSphericalCoord  = Position<Barycentric,  Horizontal>;
pub type HorizontalHeliocentricSphericalCoord = Position<Heliocentric, Horizontal>;
pub type HorizontalGeocentricSphericalCoord   = Position<Geocentric,   Horizontal>;
pub type HorizontalTopocentricSphericalCoord  = Position<Topocentric,  Horizontal>;

impl<Center: ReferenceCenter> Position<Center, Horizontal> {
    /// Creates a new horizontal spherical coordinate with constant values.
    ///
    /// # Arguments
    /// - `alt`: Altitude (α), in degrees.
    /// - `az`: Azimuth (θ), in degrees.
    /// - `distance`: Distance to the object, typically in astronomical units (AU).
    ///
    /// # Returns
    /// A new `Position` in the horizontal frame.
    pub const fn new_const(alt: Degrees, az: Degrees, distance: f64) -> Self {
        Position::new_spherical_coord(alt, az, Some(distance))
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
    pub fn new(alt: Degrees, az: Degrees, distance: f64) -> Self {
        Position::<Center, Horizontal>::new_const(
            alt.normalize_to_90_range(),
            az.normalize(),
            distance)
    }

    /// Returns the Altitude (α) in degrees.
    pub fn alt(&self) -> Degrees { self.polar }

    /// Returns the Azimuth (θ) in degrees.
    pub fn az(&self) -> Degrees { self.azimuth }
}
