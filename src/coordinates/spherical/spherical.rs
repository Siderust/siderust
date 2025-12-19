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
//! - **Center Parameters:** Each coordinate stores `C::Params` to support parameterized centers.
//!   For most centers, `Params = ()` (zero-cost). For `Topocentric`, it stores observer location.
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
//! use qtty::*;
//!
//! // Create a heliocentric ecliptic spherical position (Params=() is automatic)
//! let sph = Position::<Heliocentric, Ecliptic, AstronomicalUnit>::new(45.0 * DEG, 7.0 * DEG, 1.0);
//! println!("θ: {}, φ: {}, r: {:?}", sph.polar, sph.azimuth, sph.distance);
//! ```
//!
//! ## Type Aliases
//! - [`Position<C, F>`]: Spherical position (with distance and center).
//! - [`Direction<F>`]: Spherical direction (frame-only, no center or distance).
//!
//! ## Methods
//! - [`new(polar, azimuth, distance)`]: Construct a new coordinate.
//! - [`distance_to(other)`]: Compute Euclidean distance to another coordinate.
//! - [`angular_separation(other)`]: Compute angular separation in degrees.
//!
//! ## Display
//! Implements `Display` for readable output including center, frame, angles, and distance.

use crate::coordinates::{centers, frames, math};
use qtty::*;

use std::marker::PhantomData;

/// Represents a point in spherical coordinates with a specific reference center and frame.
///
/// # Type Parameters
/// - `C`: The reference center (e.g., `Barycentric`, `Heliocentric`, `Topocentric`).
/// - `F`: The reference frame (e.g., `frames::ICRS`, `Ecliptic`, `Horizontal`).
/// - `U`: The unit type for distance (e.g., `AstronomicalUnit`, `Kilometer`).
///
/// # Center Parameters
///
/// The coordinate stores `C::Params` to support parameterized centers:
/// - For `Barycentric`, `Heliocentric`, `Geocentric`: `Params = ()` (zero overhead)
/// - For `Topocentric`: `Params = ObserverSite` (stores observer location)
///
/// Use `new_raw()` with explicit center_params, or use convenience constructors
/// for centers with `Params = ()`.
#[derive(Debug, Clone, Copy)]
pub struct SphericalCoord<C: centers::ReferenceCenter, F: frames::ReferenceFrame, U: Unit> {
    pub polar: Degrees,        // θ (polar/latitude/declination)
    pub azimuth: Degrees,      // φ (azimuth/longitude/right ascension)
    pub distance: Quantity<U>, // Distance (AstronomicalUnits, parsec, etc.)

    center_params: C::Params,
    _frame: PhantomData<F>,
}

impl<C, F, U> SphericalCoord<C, F, U>
where
    C: centers::ReferenceCenter,
    F: frames::ReferenceFrame,
    U: Unit,
{
    /// Constructs a new spherical coordinate with explicit center parameters.
    ///
    /// * `center_params`: The center parameters (e.g., `()` or `ObserverSite`)
    /// * `polar`: angle from the reference plane, in degrees  
    /// * `azimuth`: angle from the reference direction, in degrees  
    /// * `distance`: radial distance in the same unit `U`
    pub const fn new_raw_with_params(
        center_params: C::Params,
        polar: Degrees,
        azimuth: Degrees,
        distance: Quantity<U>,
    ) -> Self {
        Self {
            polar,
            azimuth,
            distance,
            center_params,
            _frame: PhantomData,
        }
    }

    /// Returns a reference to the center parameters.
    ///
    /// For most centers this returns `&()`. For `Topocentric`, it returns `&ObserverSite`.
    pub fn center_params(&self) -> &C::Params {
        &self.center_params
    }

    /// Calculates the angular separation between this coordinate and another.
    ///
    /// # Arguments
    /// - `other`: The other spherical coordinate.
    ///
    /// # Returns
    /// The angular separation in degrees.
    pub fn angular_separation(&self, other: Self) -> Degrees {
        let az1 = self.azimuth.to::<Radian>().value();
        let po1 = self.polar.to::<Radian>().value();
        let az2 = other.azimuth.to::<Radian>().value();
        let po2 = other.polar.to::<Radian>().value();

        let angle_rad = math::geometry::angular_separation(az1, po1, az2, po2);
        Radians::new(angle_rad).to::<Degree>()
    }

    /// Returns a **direction** (unitless unit vector) corresponding to this position
    /// (i.e. same angular coordinates, radius = 1).
    ///
    /// Note: Directions are frame-only types (no center). This extracts the
    /// normalized direction regardless of the position's center.
    #[must_use]
    pub fn direction(&self) -> super::direction::Direction<F> {
        super::direction::Direction::new(self.polar, self.azimuth)
    }
}

// =============================================================================
// Convenience constructors for centers with Params = ()
// =============================================================================

impl<C, F, U> SphericalCoord<C, F, U>
where
    C: centers::ReferenceCenter<Params = ()>,
    F: frames::ReferenceFrame,
    U: Unit,
{
    /// Constructs a new spherical coordinate for centers with `Params = ()`.
    ///
    /// This is a convenience constructor that doesn't require passing `()` explicitly.
    /// Use this for `Barycentric`, `Heliocentric`, `Geocentric`, etc.
    ///
    /// * `polar`: angle from the reference plane, in degrees  
    /// * `azimuth`: angle from the reference direction, in degrees  
    /// * `distance`: radial distance in the same unit `U`
    pub const fn new_raw(polar: Degrees, azimuth: Degrees, distance: Quantity<U>) -> Self {
        Self::new_raw_with_params((), polar, azimuth, distance)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coordinates::centers::*;
    use crate::coordinates::frames;
    use crate::coordinates::spherical::SphericalCoord;

    const EPS: Degrees = Degrees::new(1e-6); // tolerance for the exact geometry cases
    const EPS_STAR: Degrees = Degrees::new(1e-2); // tolerance for the real‑world catalogue case

    /// Helper to build a unit‑length direction in the frames::ICRS frame.
    fn dir(dec: f64, ra: f64) -> SphericalCoord<Geocentric, frames::ICRS, Unitless> {
        SphericalCoord::new_raw(
            Degrees::new(dec), // polar / declination
            Degrees::new(ra),  // azimuth / right‑ascension
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
        assert!(
            (sep.to::<Degree>() - 90.0 * DEG).abs() < EPS,
            "expected 90°, got {}°",
            sep
        );
    }

    #[test]
    fn antipodal_points_give_180_degrees() {
        let a = dir(0.0, 0.0);
        let b = dir(0.0, 180.0);
        let sep = a.angular_separation(b);
        assert!(
            (sep.to::<Degree>() - 180.0 * DEG).abs() < EPS,
            "expected 180°, got {}°",
            sep
        );
    }

    #[test]
    fn polaris_betelgeuse_real_world() {
        // Star coordinates (J2000) from SIMBAD
        let polaris = dir(89.26410897, 37.95456067); // Dec, RA
        let betel = dir(7.407064, 88.792939); // Dec, RA
        let sep = polaris.angular_separation(betel);
        assert!(
            (sep.to::<Degree>() - 82.1286 * DEG).abs() < EPS_STAR,
            "expected 224882.13°, got {}°",
            sep
        );
    }

    #[test]
    fn test_spherical_coord_creation() {
        let coord = SphericalCoord::<Barycentric, frames::ICRS, Au>::new_raw(
            Degrees::new(45.0),
            Degrees::new(90.0),
            1.0 * AU,
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
            1000.0 * AU,
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
            0.0 * AU,
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
            1234.56789 * AU,
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
            100.0 * AU,
        );
        let distance = coord.distance.value();
        assert!((distance - 100.0).abs() < 1e-10);
    }

    #[test]
    fn test_spherical_coord_distance_to() {
        let coord1 = SphericalCoord::<Geocentric, frames::ICRS, Au>::new_raw(
            Degrees::new(0.0),
            Degrees::new(0.0),
            100.0 * AU,
        );
        let coord2 = SphericalCoord::<Geocentric, frames::ICRS, Au>::new_raw(
            Degrees::new(0.0),
            Degrees::new(0.0),
            200.0 * AU,
        );
        let distance = coord1.distance.value() - coord2.distance.value();
        assert_eq!(distance.abs(), 100.0);
    }

    #[test]
    fn test_spherical_coord_angular_separation() {
        let coord1 = SphericalCoord::<Geocentric, frames::ICRS, Au>::new_raw(
            Degrees::new(0.0),
            Degrees::new(0.0),
            100.0 * AU,
        );
        let coord2 = SphericalCoord::<Geocentric, frames::ICRS, Au>::new_raw(
            Degrees::new(90.0),
            Degrees::new(0.0),
            100.0 * AU,
        );
        let separation = coord1.angular_separation(coord2);
        assert!((separation.value() - 90.0).abs() < 1e-6);
    }

    #[test]
    fn test_spherical_coord_angular_separation_same_point() {
        let coord1 = SphericalCoord::<Geocentric, frames::ICRS, Au>::new_raw(
            Degrees::new(45.0),
            Degrees::new(30.0),
            100.0 * AU,
        );
        let coord2 = SphericalCoord::<Geocentric, frames::ICRS, Au>::new_raw(
            Degrees::new(45.0),
            Degrees::new(30.0),
            200.0 * AU,
        );
        let separation = coord1.angular_separation(coord2);
        assert!((separation.value() - 0.0).abs() < 1e-6);
    }

    #[test]
    fn test_spherical_coord_angular_separation_opposite() {
        let coord1 = SphericalCoord::<Geocentric, frames::ICRS, Au>::new_raw(
            Degrees::new(0.0),
            Degrees::new(0.0),
            100.0 * AU,
        );
        let coord2 = SphericalCoord::<Geocentric, frames::ICRS, Au>::new_raw(
            Degrees::new(180.0),
            Degrees::new(0.0),
            100.0 * AU,
        );
        let separation = coord1.angular_separation(coord2);
        assert!((separation.value() - 180.0).abs() < 1e-6);
    }

    #[test]
    fn test_spherical_coord_debug() {
        let coord = SphericalCoord::<Barycentric, frames::ICRS, Au>::new_raw(
            Degrees::new(45.0),
            Degrees::new(90.0),
            1.0 * AU,
        );
        let debug_str = format!("{:?}", coord);
        assert!(debug_str.contains("SphericalCoord"));
    }

    #[test]
    fn test_spherical_coord_clone() {
        let coord1 = SphericalCoord::<Geocentric, frames::ICRS, Au>::new_raw(
            Degrees::new(30.0),
            Degrees::new(60.0),
            100.0 * AU,
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
            200.0 * AU,
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
            100.0 * AU,
        );
        let coord2 = SphericalCoord::<Geocentric, frames::ICRS, Au>::new_raw(
            Degrees::new(30.0),
            Degrees::new(60.0),
            100.0 * AU,
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
            150.0 * AU,
        );
        let coord2 = SphericalCoord::<Geocentric, frames::Ecliptic, Au>::new_raw(
            Degrees::new(45.0),
            Degrees::new(90.0),
            150.0 * AU,
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
            1e6 * AU,
        );
        assert!((coord.polar.value() - 359.999).abs() < 1e-6);
        assert!((coord.azimuth.value() - 179.999).abs() < 1e-6);
        assert!((coord.distance.value() - 1e6).abs() < 1e-6);

        // Very small angles
        let coord = SphericalCoord::<Geocentric, frames::ICRS, Au>::new_raw(
            Degrees::new(0.001),
            Degrees::new(0.001),
            1e-6 * AU,
        );
        assert!((coord.polar.value() - 0.001).abs() < 1e-6);
        assert!((coord.azimuth.value() - 0.001).abs() < 1e-6);
        assert!((coord.distance.value() - 1e-6).abs() < 1e-6);
    }
}
