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
//! - **Radial distance (d)** → distance from the reference center, typically in astronomical units (AstronomicalUnits).
//!
//! Right Ascension is normalized to the [0°, 360°] range, and Declination to [-90°, 90°].
//!
//! # Example
//! ```rust
//! use siderust::coordinates::spherical::direction::Equatorial;
//! use qtty::*;
//!
//! let coord = Equatorial::new(
//!     120.0 * DEG, 45.0 * DEG
//! );
//! println!("RA = {}, Dec = {}", coord.ra(), coord.dec());
//! ```

use super::*;
use crate::coordinates::{centers::*, frames::*};
use qtty::*;

impl direction::Direction<Equatorial> {
    /// Constructs a new equatorial direction with normalized input angles.
    ///
    /// Right Ascension is normalized to [0°, 360°], Declination to [-90°, 90°].
    ///
    /// # Arguments
    /// - `ra`: Right Ascension (α), in degrees.
    /// - `dec`: Declination (δ), in degrees.
    pub fn new_equatorial(ra: Degrees, dec: Degrees) -> Self {
        Self::new(dec.wrap_quarter_fold(), ra.normalize())
    }
}

impl<C: ReferenceCenter<Params = ()>, U: LengthUnit> Position<C, Equatorial, U> {
    pub const fn new_const(ra: Degrees, dec: Degrees, distance: Quantity<U>) -> Self {
        Self::new_raw(dec, ra, distance)
    }

    /// Constructs a new equatorial spherical coordinate with normalized input angular.
    ///
    /// Right Ascension is normalized to the [0°, 360°] range, and Declination to the [-90°, 90°] range.
    ///
    /// # Arguments
    /// - `ra`: Right Ascension (α), in degrees.
    /// - `dec`: Declination (δ), in degrees.
    /// - `distance`: LengthUnit to the object, typically in astronomical units (AstronomicalUnits).
    ///
    /// # Returns
    /// A normalized `Position` in the equatorial frame.
    pub fn new<A, T>(ra: A, dec: A, distance: T) -> Self
    where
        A: Into<Degrees>,
        T: Into<Quantity<U>>,
    {
        Self::new_const(
            ra.into().normalize(),
            dec.into().wrap_quarter_fold(),
            distance.into(),
        )
    }
}

impl<C: ReferenceCenter, U: Unit> SphericalCoord<C, Equatorial, U> {
    /// Returns the Declination (δ) in degrees.
    pub fn dec(&self) -> Degrees {
        self.polar
    }

    /// Returns the Right Ascension (α) in degrees.
    pub fn ra(&self) -> Degrees {
        self.azimuth
    }
}
