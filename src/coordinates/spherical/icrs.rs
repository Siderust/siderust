//! Defines spherical coordinate types in the International Celestial Reference System (ICRS) frame,
//! using Right Ascension (RA), Declination (Dec), and distance.
//!
//! These coordinates are widely used in astrometry and celestial mechanics for describing the positions
//! of celestial objects relative to the ICRS, which is a standard celestial reference frame.
//!
//! # Coordinate Convention
//! The `Position<Center, ICRS>` type uses:
//!
//! - **Right Ascension (RA or α)** → `azimuth`: angle measured eastward along the celestial equator from the vernal equinox, in degrees.
//! - **Declination (Dec or δ)**    → `polar`: angle above or below the celestial equator, in degrees.
//! - **Radial distance (d)**       → distance from the reference center, typically in astronomical units (AU).
//!
//! Right Ascension is normalized to the [0°, 360°] range, and Declination to the [-90°, 90°] range.
//!
//! # Example
//! ```rust
//! use siderust::coordinates::spherical::position::GCRS;
//! use siderust::units::Degrees;
//!
//! let coord = GCRS::new(
//!     Degrees::new(120.0), Degrees::new(45.0), 1.0
//! );
//! println!("RA = {}, Dec = {}", coord.ra(), coord.dec());
//! ```

use super::{Position, Direction};
use crate::coordinates::spherical::SphericalCoord;
use crate::coordinates::{
    kinds::Kind,
    frames::ICRS,
    centers::*,
};
use crate::units::{Unit, Degrees};

impl<C: ReferenceCenter> Direction<C, ICRS> {
    /// Creates a new ICRS spherical coordinate with constant values.
    ///
    /// # Arguments
    /// - `ra`: Right Ascension (α), in degrees.
    /// - `dec`: Declination (δ), in degrees.
    ///
    /// # Returns
    /// A new `Direction` in the ICRS frame.
    pub const fn new_const(ra: Degrees, dec: Degrees) -> Self {
        Self::new_spherical_coord(dec, ra, None)
    }

    /// Constructs a new ICRS spherical coordinate with normalized input angular.
    ///
    /// Right Ascension is normalized to the [0°, 360°] range, and Declination to the [-90°, 90°] range.
    ///
    /// # Arguments
    /// - `ra`: Right Ascension (α), in degrees.
    /// - `dec`: Declination (δ), in degrees.
    ///
    /// # Returns
    /// A normalized `Direction` in the ICRS frame.
    pub fn new(ra: Degrees, dec: Degrees) -> Self {
        Self::new_const(
            ra.normalize(),
            dec.normalize_to_90_range()
        )
    }
}

impl<C: ReferenceCenter, U: Unit> Position<C, ICRS, U> {
    /// Creates a new ICRS spherical coordinate with constant values.
    ///
    /// # Arguments
    /// - `ra`: Right Ascension (α), in degrees.
    /// - `dec`: Declination (δ), in degrees.
    /// - `distance`: Distance to the object, typically in astronomical units (AU).
    ///
    /// # Returns
    /// A new `Position` in the ICRS frame.
    pub const fn new_const(ra: Degrees, dec: Degrees, distance: U) -> Self {
        Self::new_spherical_coord(dec, ra, Some(distance))
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
    /// A normalized `Position` in the ICRS frame.
    pub fn new(ra: Degrees, dec: Degrees, distance: U) -> Self {
        Self::new_const(
            ra.normalize(),
            dec.normalize_to_90_range(),
            distance)
    }
}

impl<C: ReferenceCenter, U: Unit, K: Kind> SphericalCoord<C, ICRS, U, K> {
    /// Returns the Declination (δ) in degrees.
    pub fn dec(&self) -> Degrees { self.polar }

    /// Returns the Right Ascension (α) in degrees.
    pub fn ra(&self) -> Degrees { self.azimuth }
}
