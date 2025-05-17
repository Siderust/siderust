//! Defines spherical coordinate types in the horizontal frame, using altitude (Alt), azimuth (Az), and distance (d).
//!
//! These coordinates are commonly used in observational astronomy to describe the position of celestial objects
//! relative to an observer's local horizon.
//!
//! # Coordinate Convention
//! The `SphericalCoord<Center, Horizontal>` type uses:
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
//! use siderust::coordinates::HorizontalTopocentricSphericalCoord;
//! use siderust::units::Degrees;
//!
//! let coord = HorizontalTopocentricSphericalCoord::new(
//!     Degrees::new(45.0), Degrees::new(120.0), 1.0
//! );
//! println!("alt = {}, az = {}", coord.alt(), coord.az());
//! ```

use super::{SphericalCoord, SphericalBuilder};
use crate::coordinates::{
    frames::*,
    centers::*,
};
use crate::units::Degrees;

// Horizontal Coordinate Types
// Polar   -> Alt (α) – the angle from the horizon. [-90°, 90°]
// Azimuth -> Az (θ) – the angle from a prime meridian. [0°, 360°]
// Radial  -> Distance (d) – the distance between the source and the target.
pub type HorizontalBarycentricSphericalCoord  = SphericalCoord<Barycentric,  Horizontal>;
pub type HorizontalHeliocentricSphericalCoord = SphericalCoord<Heliocentric, Horizontal>;
pub type HorizontalGeocentricSphericalCoord   = SphericalCoord<Geocentric,   Horizontal>;
pub type HorizontalTopocentricSphericalCoord  = SphericalCoord<Topocentric,  Horizontal>;

impl<Center: ReferenceCenter> SphericalCoord<Center, Horizontal> {
    /// Creates a new horizontal spherical coordinate with constant values.
    ///
    /// # Arguments
    /// - `alt`: Altitude (α), in degrees.
    /// - `az`: Azimuth (θ), in degrees.
    /// - `radial_distance`: Distance to the object, typically in astronomical units (AU).
    ///
    /// # Returns
    /// A new `SphericalCoord` in the horizontal frame.
    pub const fn new_const(alt: Degrees, az: Degrees, radial_distance: f64) -> Self {
        SphericalCoord::new_spherical_coord(alt, az, radial_distance)
    }

    /// Constructs a new horizontal spherical coordinate with normalized input angular.
    ///
    /// Altitude is normalized to the [-90°, 90°] range, and azimuth to the [0°, 360°] range.
    ///
    /// # Arguments
    /// - `alt`: Altitude (α), in degrees.
    /// - `az`: Azimuth (θ), in degrees.
    /// - `radial_distance`: Distance to the object, typically in astronomical units (AU).
    ///
    /// # Returns
    /// A `SphericalCoord` in the horizontal frame.
    pub fn new(alt: Degrees, az: Degrees, radial_distance: f64) -> Self {
        SphericalCoord::<Center, Horizontal>::new_const(
            alt.normalize_to_90_range(),
            az.normalize(),
            radial_distance)
    }

    /// Returns the Altitude (α) in degrees.
    pub fn alt(&self) -> Degrees { self.polar }

    /// Returns the Azimuth (θ) in degrees.
    pub fn az(&self) -> Degrees { self.azimuth }
}

/// Builds a new horizontal spherical coordinate from Alt/Az, distance.
///
/// This is a convenience constructor for `SphericalBuilder` implementations,
/// allowing generic code to instantiate horizontal coordinates.
///
/// # Arguments
/// - `alt`: Altitude (α), in degrees.
/// - `az`: Azimuth (θ), in degrees.
/// - `r`: Radial distance, typically in astronomical units (AU).
///
/// # Returns
/// A new horizontal `SphericalCoord`.
impl<Center: ReferenceCenter> SphericalBuilder<Center, Horizontal>
    for SphericalCoord<Center, Horizontal>
{
    fn build(
        alt: Degrees,
        az: Degrees,
        r: f64
    ) -> SphericalCoord<Center, Horizontal> {
        SphericalCoord::<Center, Horizontal>::new(alt, az, r)
    }
}
