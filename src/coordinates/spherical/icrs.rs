//! Defines spherical coordinate types in the International Celestial Reference System (ICRS) frame,
//! using Right Ascension (RA), Declination (Dec), and distance.
//!
//! These coordinates are widely used in astrometry and celestial mechanics for describing the positions
//! of celestial objects relative to the ICRS, which is a standard celestial reference frame.
//!
//! # Coordinate Convention
//! The `SphericalCoord<Center, ICRS>` type uses:
//!
//! - **Right Ascension (RA or α)** → `azimuth`: angle measured eastward along the celestial equator from the vernal equinox, in degrees.
//! - **Declination (Dec or δ)**    → `polar`: angle above or below the celestial equator, in degrees.
//! - **Radial distance (d)**       → distance from the reference center, typically in astronomical units (AU).
//!
//! Right Ascension is normalized to the [0°, 360°] range, and Declination to the [-90°, 90°] range.
//!
//! # Provided Types
//! The following type aliases define common combinations of center and ICRS frame:
//!
//! - `ICRSBarycentricSphericalCoord` → Center: solar system barycenter.
//! - `ICRSHeliocentricSphericalCoord` → Center: Sun.
//! - `ICRSGeocentricSphericalCoord` → Center: Earth (common for observations).
//! - `ICRSTopocentricSphericalCoord` → Center: specific observer on Earth.
//!
//! # Example
//! ```rust
//! use siderust::coordinates::ICRSGeocentricSphericalCoord;
//! use siderust::units::Degrees;
//!
//! let coord = ICRSGeocentricSphericalCoord::new(
//!     Degrees::new(120.0), Degrees::new(45.0), 1.0
//! );
//! println!("RA = {}, Dec = {}", coord.ra(), coord.dec());
//! ```

use super::SphericalCoord;
use crate::coordinates::{
    frames::*,
    centers::*,
};
use crate::units::Degrees;

// ICRS Coordinate Types
// Polar   -> Dec (δ) – the angle from a prime meridian. [-90°, 90°]
// Azimuth -> RA (α) – the angle from the equator. [0°, 360°]
// Radial  -> Distance (d) – the distance between the source and the target.
pub type ICRSBarycentricSphericalCoord  = SphericalCoord<Barycentric,  ICRS>;
pub type ICRSHeliocentricSphericalCoord = SphericalCoord<Heliocentric, ICRS>;
pub type ICRSGeocentricSphericalCoord   = SphericalCoord<Geocentric,   ICRS>;
pub type ICRSTopocentricSphericalCoord  = SphericalCoord<Topocentric,  ICRS>;

impl<Center: ReferenceCenter> SphericalCoord<Center, ICRS> {
    /// Creates a new ICRS spherical coordinate with constant values.
    ///
    /// # Arguments
    /// - `ra`: Right Ascension (α), in degrees.
    /// - `dec`: Declination (δ), in degrees.
    /// - `distance`: Distance to the object, typically in astronomical units (AU).
    ///
    /// # Returns
    /// A new `SphericalCoord` in the ICRS frame.
    pub const fn new_const(ra: Degrees, dec: Degrees, distance: f64) -> Self {
        SphericalCoord::new_spherical_coord(dec, ra, distance)
    }

    /// Constructs a new ICRS spherical coordinate with normalized input angular.
    ///
    /// Right Ascension is normalized to the [0°, 360°] range, and Declination to the [-90°, 90°] range.
    ///
    /// # Arguments
    /// - `ra`: Right Ascension (α), in degrees.
    /// - `dec`: Declination (δ), in degrees.
    /// - `distance`: Distance to the object, typically in astronomical units (AU).
    ///
    /// # Returns
    /// A normalized `SphericalCoord` in the ICRS frame.
    pub fn new(ra: Degrees, dec: Degrees, distance: f64) -> Self {
        SphericalCoord::<Center, ICRS>::new_const(
            ra.normalize(),
            dec.normalize_to_90_range(),
            distance)
    }

    /// Returns the Declination (δ) in degrees.
    pub fn dec(&self) -> Degrees { self.polar }

    /// Returns the Right Ascension (α) in degrees.
    pub fn ra(&self) -> Degrees { self.azimuth }
}
