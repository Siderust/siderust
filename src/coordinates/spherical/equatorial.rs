//! Defines spherical coordinate types in the equatorial frame using Right Ascension (RA), Declination (Dec), and distance.
//!
//! These coordinates are widely used in observational and astrometric astronomy for locating stars,
//! galaxies, and solar system bodies relative to celestial equators and poles.
//!
//! # Coordinate Convention
//! The `Position<Center, Equatorial>` type uses:
//!
//! - **Right Ascension (α)** → `azimuth`: angle measured eastward along the celestial equator from the vernal equinox, in degrees.
//! - **Declination (δ)**     → `polar`: angle above or below the celestial equator, in degrees.
//! - **Radial distance (d)** → distance from the reference center, typically in astronomical units (AU).
//!
//! Right Ascension is normalized to the [0°, 360°] range, and Declination to [-90°, 90°].
//!
//! # Provided Types
//! The following type aliases define common equatorial coordinate systems with different centers:
//!
//! - `EquatorialBarycentricSphericalCoord` → For sources relative to the solar system barycenter.
//! - `EquatorialHeliocentricSphericalCoord` → For solar system objects relative to the Sun.
//! - `EquatorialGeocentricSphericalCoord` → For objects relative to the Earth (most common in observations).
//!
//! # Example
//! ```rust
//! use siderust::coordinates::EquatorialGeocentricSphericalCoord;
//! use siderust::units::Degrees;
//!
//! let coord = EquatorialGeocentricSphericalCoord::new(
//!     Degrees::new(120.0), Degrees::new(45.0), 1.0
//! );
//! println!("RA = {}, Dec = {}", coord.ra(), coord.dec());
//! ```

use super::Position;
use crate::coordinates::{
    frames::*,
    centers::*,
};
use crate::units::Degrees;

// Equatorial Coordinate Types
// Polar   -> Dec (δ) – the angle from a prime meridian.
// Azimuth -> RA (α) – the angle from the equator.
// Radial  -> Distance (d) – the distance between the source and the target.
pub type EquatorialBarycentricSphericalCoord  = Position<Barycentric,  Equatorial>;
pub type EquatorialHeliocentricSphericalCoord = Position<Heliocentric, Equatorial>;
pub type EquatorialGeocentricSphericalCoord   = Position<Geocentric,   Equatorial>;
pub type EquatorialTopocentricSphericalCoord  = Position<Topocentric,  Equatorial>;

impl<Center: ReferenceCenter> Position<Center, Equatorial> {
    pub const fn new_const(ra: Degrees, dec: Degrees, distance: f64) -> Self {
        Position::new_spherical_coord(dec, ra, Some(distance))
    }

    /// Constructs a new equatorial spherical coordinate with normalized input angular.
    ///
    /// Right Ascension is normalized to the [0°, 360°] range, and Declination to the [-90°, 90°] range.
    ///
    /// # Arguments
    /// - `ra`: Right Ascension (α), in degrees.
    /// - `dec`: Declination (δ), in degrees.
    /// - `distance`: Distance to the object, typically in astronomical units (AU).
    ///
    /// # Returns
    /// A normalized `Position` in the equatorial frame.
    pub fn new(ra: Degrees, dec: Degrees, distance: f64) -> Self {
        Position::<Center, Equatorial>::new_const(
            ra.normalize(),
            dec.normalize_to_90_range(),
            distance)
    }

    /// Returns the Declination (δ) in degrees.
    pub fn dec(&self) -> Degrees { self.polar }

    /// Returns the Right Ascension (α) in degrees.
    pub fn ra(&self)  -> Degrees { self.azimuth }

}
