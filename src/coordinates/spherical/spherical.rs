//! # Spherical Coordinates
//!
//! This module defines the generic [`SphericalCoord<C, F, K>`] type for representing positions or directions
//! in spherical coordinates, parameterized by astronomical reference centers and frames for strong type safety.
//!
//! ## Overview
//!
//! - **Generic over Center, Frame, and Kind:**
//!   - `C`: Reference center (e.g., `Heliocentric`, `Geocentric`, `Barycentric`).
//!   - `F`: Reference frame (e.g., `ICRS`, `Ecliptic`, `Equatorial`).
//!   - `K`: Kind marker (`PositionKind`, `DirectionKind`), enforcing semantic correctness.
//! - **Type Safety:** Operations are only allowed between coordinates with matching type parameters.
//! - **Units:** Angles are stored as [`Degrees`]; distance is optional and typically in AU or parsecs (see context).
//! - **Conversions:** Seamless conversion to and from [`CartesianCoord`] via `From`/`Into`.
//! - **Operations:** Compute Euclidean distance and angular separation between coordinates.
//!
//! ## Example
//! ```rust
//! use siderust::coordinates::spherical::Position;
//! use siderust::coordinates::centers::Heliocentric;
//! use siderust::coordinates::frames::Ecliptic;
//! use siderust::units::Degrees;
//!
//! // Create a heliocentric ecliptic spherical position
//! let sph = Position::<Heliocentric, Ecliptic>::new(Degrees::new(45.0), Degrees::new(7.0), 1.0);
//! println!("θ: {}, φ: {}, r: {:?}", sph.polar, sph.azimuth, sph.distance);
//! ```
//!
//! ## Type Aliases
//! - [`Position<C, F>`]: Spherical position (with distance).
//! - [`Direction<C, F>`]: Spherical direction (distance is typically `None`).
//!
//! ## Methods
//! - [`new(polar, azimuth, distance)`]: Construct a new coordinate.
//! - [`from_degrees(polar, azimuth, distance)`]: Construct from primitive values.
//! - [`distance_to(other)`]: Compute Euclidean distance to another coordinate.
//! - [`angular_separation(other)`]: Compute angular separation in degrees.
//!
//! ## Display
//! Implements `Display` for readable output including center, frame, angles, and distance.

use crate::units::{Degrees, Radians};

use crate::coordinates::{
    cartesian,
    frames::*,
    centers::*,
    kinds::*,
};
use std::marker::PhantomData;

/// Represents a point in spherical coordinates with a specific reference center, frame, and kind.
///
/// # Type Parameters
/// - `C`: The reference center (e.g., `Barycentric`, `Heliocentric`).
/// - `F`: The reference frame (e.g., `ICRS`, `Ecliptic`).
/// - `K`: The kind marker (`PositionKind`, `DirectionKind`).
#[derive(Debug, Clone, Copy)]
pub struct SphericalCoord<C: ReferenceCenter, F: ReferenceFrame, K: Kind> {
    pub polar: Degrees,      // θ (polar/latitude/declination)
    pub azimuth: Degrees,    // φ (azimuth/longitude/right ascension)
    pub distance: Option<f64>, // Optional distance (AU, parsec, etc.)

    _center: PhantomData<C>,
    _frame: PhantomData<F>,
    _kind: PhantomData<K>,
}

/// Spherical position (with distance).
pub type Position<C, F>  = SphericalCoord<C, F, PositionKind>;
/// Spherical direction (distance is typically `None`).
pub type Direction<C, F> = SphericalCoord<C, F, DirectionKind>;

impl<C, F, K> SphericalCoord<C, F, K>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
    K: Kind,
    cartesian::CartesianCoord<C, F, K>: for<'a> From<&'a Self>,
{
    /// The zero point (origin) in this coordinate system.
    pub const CENTER: Self = Self::from_degrees(0.0, 0.0, Some(0.0));

    /// Creates a new spherical coordinate from angle types and optional distance.
    pub const fn new_spherical_coord(polar: Degrees, azimuth: Degrees, distance: Option<f64>) -> Self {
        Self {
            polar,
            azimuth,
            distance,
            _center: PhantomData,
            _frame: PhantomData,
            _kind: PhantomData,
        }
    }

    /// Creates a new spherical coordinate from primitive values (degrees).
    pub const fn from_degrees(polar: f64, azimuth: f64, r: Option<f64>) -> Self{
        Self::new_spherical_coord(Degrees::new(polar), Degrees::new(azimuth), r)
    }

    /// Calculates the Euclidean distance to another spherical coordinate.
    ///
    /// # Arguments
    /// - `other`: The other spherical coordinate.
    ///
    /// # Returns
    /// The distance to the other coordinate.
    pub fn distance_to(&self, other: &Self) -> f64 {
        self.to_cartesian().distance_to(&other.to_cartesian())
    }

    /// Calculates the angular separation between this coordinate and another.
    ///
    /// # Arguments
    /// - `other`: The other spherical coordinate.
    ///
    /// # Returns
    /// The angular separation in degrees.
    pub fn angular_separation(&self, other: Self) -> Degrees {
        let az1 = self.azimuth.to_radians();
        let po1 = self.polar.to_radians();
        let az2 = other.azimuth.to_radians();
        let po2 = other.polar.to_radians();

        let x = (po1.cos() * po2.sin()) - (po1.sin() * po2.cos() * (az2 - az1).cos());
        let y = po2.cos() * (az2 - az1).sin();
        let z = (po1.sin() * po2.sin()) + (po1.cos() * po2.cos() * (az2 - az1).cos());

        let angle_rad = (x * x + y * y).sqrt().atan2(z);
        Radians::new(angle_rad).to_degrees()
    }
}

impl<C, F, K> std::fmt::Display for SphericalCoord<C, F, K>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
    K: Kind,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Center: {}, Frame: {}, θ: {:.6}, φ: {:.6}, r: {:.6}",
            C::center_name(),
            F::frame_name(),
            self.polar, self.azimuth, self.distance.unwrap_or(std::f64::NAN)
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::units::Degrees;
    
    #[test]
    fn test_spherical_coord_creation() {
        let coord = Position::<Barycentric, ICRS>::new(Degrees::new(45.0), Degrees::new(90.0), 1.0);
        assert_eq!(coord.ra().as_f64(), 45.0);
        assert_eq!(coord.dec().as_f64(), 90.0);
        assert_eq!(coord.distance.unwrap(), 1.0);
    }
    
    #[test]
    fn test_spherical_coord_to_string() {
        let coord = Position::<Geocentric, ICRS>::new(Degrees::new(30.0), Degrees::new(60.0), 1000.0);
        let coord_string = coord.to_string();
        assert!(coord_string.contains("θ: 60"));
        assert!(coord_string.contains("φ: 30"));
        assert!(coord_string.contains("r: 1000"));
    }
    
    #[test]
    fn test_spherical_coord_zero_values() {
        let coord = Position::<Heliocentric, ICRS>::new(Degrees::new(0.0), Degrees::new(0.0), 0.0);
        assert_eq!(coord.polar.as_f64(), 0.0);
        assert_eq!(coord.azimuth.as_f64(), 0.0);
        assert_eq!(coord.distance.unwrap(), 0.0);
    }

    #[test]
    fn test_spherical_coord_precision() {
        let coord = Position::<Barycentric, ICRS>::new(Degrees::new(90.654321), Degrees::new(45.123456), 1234.56789);
        assert!((coord.dec().as_f64() - 45.123456).abs() < 1e-6);
        assert!((coord.ra().as_f64() - 90.654321).abs() < 1e-6);
        assert!((coord.distance.unwrap() - 1234.56789).abs() < 1e-6);
    }
}
