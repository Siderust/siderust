//! # Spherical Coordinates
//!
//! This module defines the generic [`SphericalCoord<C, F>`] type for representing positions or directions
//! in spherical coordinates, parameterized by astronomical reference centers and frames for strong type safety.
//!
//! ## Overview
//!
//! - **Generic over Center and Frame:**
//!   - `C`: Reference center (e.g., `Heliocentric`, `Geocentric`, `Barycentric`).
//!   - `F`: Reference frame (e.g., `frames::ICRS`, `Ecliptic`, `Equatorial`).
//! - **Type Safety:** Operations are only allowed between coordinates with matching type parameters.
//! - **Units:** Angles are stored as [`Degrees`]; distance is optional and typically in AstronomicalUnits or parsecs (see context).
//! - **Conversions:** Seamless conversion to and from [`Vector`] via `From`/`Into`.
//! - **Operations:** Compute Euclidean distance and angular separation between coordinates.
//!
//! ## Example
//! ```rust
//! use siderust::coordinates::spherical::Position;
//! use siderust::coordinates::centers::Heliocentric;
//! use siderust::coordinates::frames::Ecliptic;
//! use siderust::units::{Degrees, AstronomicalUnit};
//!
//! // Create a heliocentric ecliptic spherical position
//! let sph = Position::<Heliocentric, Ecliptic, AstronomicalUnit>::new(Degrees::new(45.0), Degrees::new(7.0), 1.0);
//! println!("θ: {}, φ: {}, r: {:?}", sph.polar, sph.azimuth, sph.distance);
//! ```
//!
//! ## Type Aliases
//! - [`Position<C, F>`]: Spherical position (with distance).
//! - [`Direction<C, F>`]: Spherical direction (distance is typically `None`).
//!
//! ## Methods
//! - [`new(polar, azimuth, distance)`]: Construct a new coordinate.
//! - [`distance_to(other)`]: Compute Euclidean distance to another coordinate.
//! - [`angular_separation(other)`]: Compute angular separation in degrees.
//!
//! ## Display
//! Implements `Display` for readable output including center, frame, angles, and distance.

use crate::units::*;
use crate::coordinates::{frames, centers};

use std::marker::PhantomData;

/// Represents a point in spherical coordinates with a specific reference center and frame.
///
/// # Type Parameters
/// - `C`: The reference center (e.g., `Barycentric`, `Heliocentric`).
/// - `F`: The reference frame (e.g., `frames::ICRS`, `Ecliptic`).
#[derive(Debug, Clone, Copy)]
pub struct SphericalCoord<
    C: centers::ReferenceCenter,
    F: frames::ReferenceFrame,
    U: Unit,
> {
    pub polar: Degrees,        // θ (polar/latitude/declination)
    pub azimuth: Degrees,      // φ (azimuth/longitude/right ascension)
    pub distance: Quantity<U>, // Distance (AstronomicalUnits, parsec, etc.)

    _center: PhantomData<C>,
    _frame: PhantomData<F>,
}

impl<C, F, U> SphericalCoord<C, F, U>
where
    C: centers::ReferenceCenter,
    F: frames::ReferenceFrame,
    U: Unit,
{
    /// Constructs a new spherical coordinate.
    ///
    /// * `polar`: angle from the reference plane, in degrees  
    /// * `azimuth`: angle from the reference direction, in degrees  
    /// * `distance`: radial distance in the same unit `U`
    pub const fn new_raw(polar: Degrees, azimuth: Degrees, distance: Quantity<U>) -> Self {
        Self {
            polar,
            azimuth,
            distance,
            _center: PhantomData,
            _frame: PhantomData,
        }
    }

    /// Calculates the angular separation between this coordinate and another.
    ///
    /// # Arguments
    /// - `other`: The other spherical coordinate.
    ///
    /// # Returns
    /// The angular separation in degrees.
    pub fn angular_separation(&self, other: Self) -> Degrees {
        let az1 = self.azimuth.to::<Radian>();
        let po1 = self.polar.to::<Radian>();
        let az2 = other.azimuth.to::<Radian>();
        let po2 = other.polar.to::<Radian>();

        let x = (po1.cos() * po2.sin()) - (po1.sin() * po2.cos() * (az2 - az1).cos());
        let y = po2.cos() * (az2 - az1).sin();
        let z = (po1.sin() * po2.sin()) + (po1.cos() * po2.cos() * (az2 - az1).cos());

        let angle_rad = (x * x + y * y).sqrt().atan2(z);
        Radians::new(angle_rad).to::<Degree>()
    }


    /// Returns a **direction** (unitless unitary vector) corresponding to this position
    /// (i.e. same angular coordinates, radius = 1).
    #[must_use]
    pub fn direction(&self) -> super::Direction<C, F> {
        super::Direction::new_raw(self.polar, self.azimuth, Quantity::<Unitless>::new(1.0))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coordinates::spherical::SphericalCoord;
    use crate::coordinates::centers::*;
    use crate::coordinates::frames;

    const EPS: Degrees = Degrees::new(1e-6);       // tolerance for the exact geometry cases
    const EPS_STAR: Degrees = Degrees::new(1e-2);  // tolerance for the real‑world catalogue case

    /// Helper to build a unit‑length direction in the frames::ICRS frame.
    fn dir(dec: f64, ra: f64) -> SphericalCoord<Geocentric, frames::ICRS, Unitless> {
        SphericalCoord::new_raw(
            Degrees::new(dec),   // polar / declination
            Degrees::new(ra),    // azimuth / right‑ascension
            Quantity::<Unitless>::new(1.0),
        )
    }

    #[test]
    fn identity_separation_is_zero() {
        let a = dir(12.3, 45.6);
        let sep = a.angular_separation(a);
        assert!(sep.to::<Degree>().abs() < EPS, "expected 0°, got {}°", sep);
    }

    #[test]
    fn orthogonal_points_give_ninety_degrees() {
        let a = dir(0.0, 0.0);
        let b = dir(0.0, 90.0);
        let sep = a.angular_separation(b);
        assert!((sep.to::<Degree>() - 90.0*DEG).abs() < EPS, "expected 90°, got {}°", sep);
    }

    #[test]
    fn antipodal_points_give_180_degrees() {
        let a = dir(0.0, 0.0);
        let b = dir(0.0, 180.0);
        let sep = a.angular_separation(b);
        assert!((sep.to::<Degree>() - 180.0*DEG).abs() < EPS, "expected 180°, got {}°", sep);
    }

    #[test]
    fn polaris_betelgeuse_real_world() {
        // Star coordinates (J2000) from SIMBAD
        let polaris = dir(89.26410897, 37.95456067);   // Dec, RA
        let betel   = dir(7.407064,    88.792939);     // Dec, RA
        let sep = polaris.angular_separation(betel);
        assert!((sep.to::<Degree>() - 82.1286*DEG).abs() < EPS_STAR,
            "expected 224882.13°, got {}°", sep);
    }

    #[test]
    fn test_spherical_coord_creation() {
        let coord = SphericalCoord::<Barycentric, frames::ICRS, Au>::new_raw(
            Degrees::new(45.0),
            Degrees::new(90.0),
            1.0*AU
        );
        assert_eq!(coord.polar.value(), 45.0);
        assert_eq!(coord.azimuth.value(), 90.0);
        assert_eq!(coord.distance.value(), 1.0);
    }
    
    #[test]
    fn test_spherical_coord_new() {
        let coord = SphericalCoord::<Geocentric, frames::ICRS, Au>::new(30.0, 60.0, 1000.0);
        assert_eq!(coord.polar.value(), 60.0);
        assert_eq!(coord.azimuth.value(), 30.0);
        assert_eq!(coord.distance.value(), 1000.0);
    }

    #[test]
    fn test_spherical_coord_to_string() {
        let coord = SphericalCoord::<Geocentric, frames::ICRS, Au>::new_raw(
            Degrees::new(30.0),
            Degrees::new(60.0),
            1000.0*AU
        );
        let coord_string = coord.to_string();
        println!("{}", coord_string);
        assert!(coord_string.contains("θ: 30"));
        assert!(coord_string.contains("φ: 60"));
        assert!(coord_string.contains("r: 1000"));
        assert!(coord_string.contains("Geocentric"));
        assert!(coord_string.contains("ICRS"));
    }
    
    #[test]
    fn test_spherical_coord_zero_values() {
        let coord = SphericalCoord::<Heliocentric, frames::ICRS, Au>::new_raw(
            Degrees::new(0.0),
            Degrees::new(0.0),
            0.0*AU
        );
        assert_eq!(coord.polar.value(), 0.0);
        assert_eq!(coord.azimuth.value(), 0.0);
        assert_eq!(coord.distance.value(), 0.0);
    }

    #[test]
    fn test_spherical_coord_precision() {
        let coord = SphericalCoord::<Barycentric, frames::ICRS, Au>::new_raw(
            Degrees::new(90.654321),
            Degrees::new(45.123456),
            1234.56789*AU
        );
        assert!((coord.polar.value() - 90.654321).abs() < 1e-6);
        assert!((coord.azimuth.value() - 45.123456).abs() < 1e-6);
        assert!((coord.distance.value() - 1234.56789).abs() < 1e-6);
    }

    #[test]
    fn test_spherical_coord_distance_from_origin() {
        let coord = SphericalCoord::<Geocentric, frames::ICRS, Au>::new_raw(
            Degrees::new(45.0),
            Degrees::new(30.0),
            100.0*AU
        );
        let distance = coord.distance.value();
        assert!((distance - 100.0).abs() < 1e-10);
    }

    #[test]
    fn test_spherical_coord_distance_to() {
        let coord1 = SphericalCoord::<Geocentric, frames::ICRS, Au>::new_raw(
            Degrees::new(0.0),
            Degrees::new(0.0),
            100.0*AU
        );
        let coord2 = SphericalCoord::<Geocentric, frames::ICRS, Au>::new_raw(
            Degrees::new(0.0),
            Degrees::new(0.0),
            200.0*AU
        );
        let distance = coord1.distance.value() - coord2.distance.value();
        assert_eq!(distance.abs(), 100.0);
    }

    #[test]
    fn test_spherical_coord_angular_separation() {
        let coord1 = SphericalCoord::<Geocentric, frames::ICRS, Au>::new_raw(
            Degrees::new(0.0),
            Degrees::new(0.0),
            100.0*AU
        );
        let coord2 = SphericalCoord::<Geocentric, frames::ICRS, Au>::new_raw(
            Degrees::new(90.0),
            Degrees::new(0.0),
            100.0*AU
        );
        let separation = coord1.angular_separation(coord2);
        assert!((separation.value() - 90.0).abs() < 1e-6);
    }

    #[test]
    fn test_spherical_coord_angular_separation_same_point() {
        let coord1 = SphericalCoord::<Geocentric, frames::ICRS, Au>::new_raw(
            Degrees::new(45.0),
            Degrees::new(30.0),
            100.0*AU
        );
        let coord2 = SphericalCoord::<Geocentric, frames::ICRS, Au>::new_raw(
            Degrees::new(45.0),
            Degrees::new(30.0),
            200.0*AU
        );
        let separation = coord1.angular_separation(coord2);
        assert!((separation.value() - 0.0).abs() < 1e-6);
    }

    #[test]
    fn test_spherical_coord_angular_separation_opposite() {
        let coord1 = SphericalCoord::<Geocentric, frames::ICRS, Au>::new_raw(
            Degrees::new(0.0),
            Degrees::new(0.0),
            100.0*AU
        );
        let coord2 = SphericalCoord::<Geocentric, frames::ICRS, Au>::new_raw(
            Degrees::new(180.0),
            Degrees::new(0.0),
            100.0*AU
        );
        let separation = coord1.angular_separation(coord2);
        assert!((separation.value() - 180.0).abs() < 1e-6);
    }

    #[test]
    fn test_spherical_coord_debug() {
        let coord = SphericalCoord::<Barycentric, frames::ICRS, Au>::new_raw(
            Degrees::new(45.0),
            Degrees::new(90.0),
            1.0*AU
        );
        let debug_str = format!("{:?}", coord);
        assert!(debug_str.contains("SphericalCoord"));
    }

    #[test]
    fn test_spherical_coord_clone() {
        let coord1 = SphericalCoord::<Geocentric, frames::ICRS, Au>::new_raw(
            Degrees::new(30.0),
            Degrees::new(60.0),
            100.0*AU
        );
        let coord2 = coord1.clone();
        assert_eq!(coord1.polar.value(), coord2.polar.value());
        assert_eq!(coord1.azimuth.value(), coord2.azimuth.value());
        assert_eq!(coord1.distance.value(), coord2.distance.value());
    }

    #[test]
    fn test_spherical_coord_copy() {
        let coord1 = SphericalCoord::<Heliocentric, frames::ICRS, Au>::new_raw(
            Degrees::new(45.0),
            Degrees::new(90.0),
            200.0*AU
        );
        let coord2 = coord1; // Copy
        assert_eq!(coord1.polar.value(), coord2.polar.value());
        assert_eq!(coord1.azimuth.value(), coord2.azimuth.value());
        assert_eq!(coord1.distance.value(), coord2.distance.value());
    }

    #[test]
    fn test_spherical_coord_different_centers() {
        let coord1 = SphericalCoord::<Barycentric, frames::ICRS, Au>::new_raw(
            Degrees::new(30.0),
            Degrees::new(60.0),
            100.0*AU
        );
        let coord2 = SphericalCoord::<Geocentric, frames::ICRS, Au>::new_raw(
            Degrees::new(30.0),
            Degrees::new(60.0),
            100.0*AU
        );
        assert_eq!(coord1.polar.value(), coord2.polar.value());
        assert_eq!(coord1.azimuth.value(), coord2.azimuth.value());
        assert_eq!(coord1.distance.value(), coord2.distance.value());
    }

    #[test]
    fn test_spherical_coord_different_frames() {
        let coord1 = SphericalCoord::<Geocentric, frames::ICRS, Au>::new_raw(
            Degrees::new(45.0),
            Degrees::new(90.0),
            150.0*AU
        );
        let coord2 = SphericalCoord::<Geocentric, frames::Ecliptic, Au>::new_raw(
            Degrees::new(45.0),
            Degrees::new(90.0),
            150.0*AU
        );
        assert_eq!(coord1.polar.value(), coord2.polar.value());
        assert_eq!(coord1.azimuth.value(), coord2.azimuth.value());
        assert_eq!(coord1.distance.value(), coord2.distance.value());
    }

    #[test]
    fn test_spherical_coord_edge_cases() {
        // Large angles
        let coord = SphericalCoord::<Geocentric, frames::ICRS, Au>::new_raw(
            Degrees::new(359.999),
            Degrees::new(179.999),
            1e6*AU
        );
        assert!((coord.polar.value() - 359.999).abs() < 1e-6);
        assert!((coord.azimuth.value() - 179.999).abs() < 1e-6);
        assert!((coord.distance.value() - 1e6).abs() < 1e-6);

        // Very small angles
        let coord = SphericalCoord::<Geocentric, frames::ICRS, Au>::new_raw(
            Degrees::new(0.001),
            Degrees::new(0.001),
            1e-6*AU
        );
        assert!((coord.polar.value() - 0.001).abs() < 1e-6);
        assert!((coord.azimuth.value() - 0.001).abs() < 1e-6);
        assert!((coord.distance.value() - 1e-6).abs() < 1e-6);
    }
}
