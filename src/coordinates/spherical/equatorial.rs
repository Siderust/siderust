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
//! # Example
//! ```rust
//! use siderust::coordinates::spherical::position::Equatorial;
//! use siderust::units::Degrees;
//!
//! let coord = Equatorial::new(
//!     Degrees::new(120.0), Degrees::new(45.0), 1.0
//! );
//! println!("RA = {}, Dec = {}", coord.ra(), coord.dec());
//! ```

use super::*;
use crate::coordinates::{
    frames::*,
    centers::*,
    kinds::Kind,
};
use crate::units::{LengthUnit, Degrees, Quantity};

impl<C: ReferenceCenter> Direction<C, Equatorial> {
    pub const fn new_const(ra: Degrees, dec: Degrees) -> Self {
        Self::new_spherical_coord(dec, ra, None)
    }

    /// Constructs a new equatorial spherical coordinate with normalized input angular.
    ///
    /// Right Ascension is normalized to the [0°, 360°] range, and Declination to the [-90°, 90°] range.
    ///
    /// # Arguments
    /// - `ra`: Right Ascension (α), in degrees.
    /// - `dec`: Declination (δ), in degrees.
    ///
    /// # Returns
    /// A normalized `Direction` in the equatorial frame.
    pub fn new(ra: Degrees, dec: Degrees) -> Self {
        Self::new_const(
            ra.normalize(),
            dec.normalize_to_90_range()
        )
    }
}

impl<C: ReferenceCenter, U: LengthUnit> Position<C, Equatorial, U> {
    pub const fn new_const(ra: Degrees, dec: Degrees, distance: Quantity<U>) -> Self {
        Self::new_spherical_coord(dec, ra, Some(distance))
    }

    /// Constructs a new equatorial spherical coordinate with normalized input angular.
    ///
    /// Right Ascension is normalized to the [0°, 360°] range, and Declination to the [-90°, 90°] range.
    ///
    /// # Arguments
    /// - `ra`: Right Ascension (α), in degrees.
    /// - `dec`: Declination (δ), in degrees.
    /// - `distance`: LengthUnit to the object, typically in astronomical units (AU).
    ///
    /// # Returns
    /// A normalized `Position` in the equatorial frame.
    pub fn new(ra: Degrees, dec: Degrees, distance: Quantity<U>) -> Self {
        Self::new_const(
            ra.normalize(),
            dec.normalize_to_90_range(),
            distance)
    }
}

impl<C: ReferenceCenter, U: LengthUnit, K: Kind> SphericalCoord<C, Equatorial, U, K> {
    /// Returns the Declination (δ) in degrees.
    pub fn dec(&self) -> Degrees { self.polar }

    /// Returns the Right Ascension (α) in degrees.
    pub fn ra(&self)  -> Degrees { self.azimuth }
}
