use crate::units::{Degrees, Radians};

use crate::coordinates::{
    frames::*,
    centers::*,
};
use std::marker::PhantomData;

/// Represents a point in spherical coordinates with a specific reference center and frame.
///
/// # Type Parameters
/// - `Center`: The reference center (e.g., Barycentric, Heliocentric).
/// - `Frame`: The reference frame (e.g., ICRS, Ecliptic).
#[derive(Debug, Clone, Copy)]
pub struct SphericalCoord<Center: ReferenceCenter, Frame: ReferenceFrame> {
    /// The polar angle (θ), in degrees.
    pub polar: Degrees,
    /// The azimuthal angle (φ), in degrees.
    pub azimuth: Degrees,
    /// The radial distance (r).
    pub radial_distance: f64,

    _center: PhantomData<Center>,
    _frame: PhantomData<Frame>,
}

impl<Center, Frame> SphericalCoord<Center, Frame>
where
    Center: ReferenceCenter,
    Frame: ReferenceFrame,
    SphericalCoord<Center, Frame>: SphericalBuilder<Center, Frame>,
{

    pub const fn new_spherical_coord(polar: Degrees, azimuth: Degrees, radial_distance: f64) -> Self {
        Self {
            polar,
            azimuth,
            radial_distance,
            _center: PhantomData,
            _frame: PhantomData,
        }
    }
    
    pub const fn from_degrees(polar: f64, azimuth: f64, r: f64) -> Self{
        Self::new_spherical_coord(Degrees::new(polar), Degrees::new(azimuth), r)
    }

    /// Calculates the Euclidean distance from the origin.
    ///
    /// # Returns
    /// The distance from the origin.
    pub fn distance_from_origin(&self) -> f64 {
        self.to_cartesian().distance_from_origin()
    }

    /// Calculates the Euclidean distance to another spherical coordinate.
    ///
    /// # Arguments
    /// - `other`: The other spherical coordinate.
    ///
    /// # Returns
    /// The distance to the other coordinate.
    pub fn distance_to(&self, other: &SphericalCoord<Center, Frame>) -> f64 {
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

impl<Center, Frame> std::fmt::Display for SphericalCoord<Center, Frame>
where
    Center: crate::coordinates::centers::ReferenceCenter,
    Frame: crate::coordinates::frames::ReferenceFrame,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Center: {}, Frame: {}, θ: {:.6}, φ: {:.6}, r: {:.6}",
            Center::center_name(),
            Frame::frame_name(),
            self.polar, self.azimuth, self.radial_distance
        )
    }
}


/// A trait for building spherical coordinates in a structured way.
///
/// # Type Parameters
/// - `Center`: The reference center (e.g., Barycentric, Heliocentric).
/// - `Frame`: The reference frame (e.g., ICRS, Ecliptic).
pub trait SphericalBuilder<Center: ReferenceCenter, Frame: ReferenceFrame> {
    /// Builds a new spherical coordinate.
    ///
    /// # Arguments
    /// - `polar`: The polar angle (θ), in degrees.
    /// - `azimuth`: The azimuthal angle (φ), in degrees.
    /// - `r`: The radial distance.
    ///
    /// # Returns
    /// A new `SphericalCoord` instance.
    fn build(
        polar: Degrees,
        azimuth: Degrees,
        r: f64
    ) -> SphericalCoord<Center, Frame>;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::units::Degrees;
    
    #[test]
    fn test_spherical_coord_creation() {
        let coord = SphericalCoord::<Barycentric, ICRS>::new_spherical_coord(
            Degrees::new(45.0),
            Degrees::new(90.0),
            1.0
        );
        assert_eq!(coord.polar.as_f64(), 45.0);
        assert_eq!(coord.azimuth.as_f64(), 90.0);
        assert_eq!(coord.radial_distance, 1.0);
    }
    
    #[test]
    fn test_spherical_coord_from_degrees() {
        let coord = SphericalCoord::<Geocentric, ICRS>::from_degrees(30.0, 60.0, 1000.0);
        assert_eq!(coord.polar.as_f64(), 30.0);
        assert_eq!(coord.azimuth.as_f64(), 60.0);
        assert_eq!(coord.radial_distance, 1000.0);
    }

    #[test]
    fn test_spherical_coord_to_string() {
        let coord = SphericalCoord::<Geocentric, ICRS>::new_spherical_coord(
            Degrees::new(30.0),
            Degrees::new(60.0),
            1000.0
        );
        let coord_string = coord.to_string();
        assert!(coord_string.contains("θ: 30.000000"));
        assert!(coord_string.contains("φ: 60.000000"));
        assert!(coord_string.contains("r: 1000.000000"));
        assert!(coord_string.contains("Geocentric"));
        assert!(coord_string.contains("ICRS"));
    }
    
    #[test]
    fn test_spherical_coord_zero_values() {
        let coord = SphericalCoord::<Heliocentric, ICRS>::new_spherical_coord(
            Degrees::new(0.0),
            Degrees::new(0.0),
            0.0
        );
        assert_eq!(coord.polar.as_f64(), 0.0);
        assert_eq!(coord.azimuth.as_f64(), 0.0);
        assert_eq!(coord.radial_distance, 0.0);
    }

    #[test]
    fn test_spherical_coord_precision() {
        let coord = SphericalCoord::<Barycentric, ICRS>::new_spherical_coord(
            Degrees::new(90.654321),
            Degrees::new(45.123456),
            1234.56789
        );
        assert!((coord.polar.as_f64() - 90.654321).abs() < 1e-6);
        assert!((coord.azimuth.as_f64() - 45.123456).abs() < 1e-6);
        assert!((coord.radial_distance - 1234.56789).abs() < 1e-6);
    }

    #[test]
    fn test_spherical_coord_distance_from_origin() {
        let coord = SphericalCoord::<Geocentric, ICRS>::new_spherical_coord(
            Degrees::new(45.0),
            Degrees::new(30.0),
            100.0
        );
        let distance = coord.distance_from_origin();
        assert!((distance - 100.0).abs() < 1e-10);
    }

    #[test]
    fn test_spherical_coord_distance_to() {
        let coord1 = SphericalCoord::<Geocentric, ICRS>::new_spherical_coord(
            Degrees::new(0.0),
            Degrees::new(0.0),
            100.0
        );
        let coord2 = SphericalCoord::<Geocentric, ICRS>::new_spherical_coord(
            Degrees::new(0.0),
            Degrees::new(0.0),
            200.0
        );
        let distance = coord1.distance_to(&coord2);
        assert!((distance - 100.0).abs() < 1e-10);
    }

    #[test]
    fn test_spherical_coord_angular_separation() {
        let coord1 = SphericalCoord::<Geocentric, ICRS>::new_spherical_coord(
            Degrees::new(0.0),
            Degrees::new(0.0),
            100.0
        );
        let coord2 = SphericalCoord::<Geocentric, ICRS>::new_spherical_coord(
            Degrees::new(90.0),
            Degrees::new(0.0),
            100.0
        );
        let separation = coord1.angular_separation(coord2);
        assert!((separation.as_f64() - 90.0).abs() < 1e-6);
    }

    #[test]
    fn test_spherical_coord_angular_separation_same_point() {
        let coord1 = SphericalCoord::<Geocentric, ICRS>::new_spherical_coord(
            Degrees::new(45.0),
            Degrees::new(30.0),
            100.0
        );
        let coord2 = SphericalCoord::<Geocentric, ICRS>::new_spherical_coord(
            Degrees::new(45.0),
            Degrees::new(30.0),
            200.0
        );
        let separation = coord1.angular_separation(coord2);
        assert!((separation.as_f64() - 0.0).abs() < 1e-6);
    }

    #[test]
    fn test_spherical_coord_angular_separation_opposite() {
        let coord1 = SphericalCoord::<Geocentric, ICRS>::new_spherical_coord(
            Degrees::new(0.0),
            Degrees::new(0.0),
            100.0
        );
        let coord2 = SphericalCoord::<Geocentric, ICRS>::new_spherical_coord(
            Degrees::new(180.0),
            Degrees::new(0.0),
            100.0
        );
        let separation = coord1.angular_separation(coord2);
        assert!((separation.as_f64() - 180.0).abs() < 1e-6);
    }

    #[test]
    fn test_spherical_coord_debug() {
        let coord = SphericalCoord::<Barycentric, ICRS>::new_spherical_coord(
            Degrees::new(45.0),
            Degrees::new(90.0),
            1.0
        );
        let debug_str = format!("{:?}", coord);
        assert!(debug_str.contains("SphericalCoord"));
    }

    #[test]
    fn test_spherical_coord_clone() {
        let coord1 = SphericalCoord::<Geocentric, ICRS>::new_spherical_coord(
            Degrees::new(30.0),
            Degrees::new(60.0),
            100.0
        );
        let coord2 = coord1.clone();
        assert_eq!(coord1.polar.as_f64(), coord2.polar.as_f64());
        assert_eq!(coord1.azimuth.as_f64(), coord2.azimuth.as_f64());
        assert_eq!(coord1.radial_distance, coord2.radial_distance);
    }

    #[test]
    fn test_spherical_coord_copy() {
        let coord1 = SphericalCoord::<Heliocentric, ICRS>::new_spherical_coord(
            Degrees::new(45.0),
            Degrees::new(90.0),
            200.0
        );
        let coord2 = coord1; // Copy
        assert_eq!(coord1.polar.as_f64(), coord2.polar.as_f64());
        assert_eq!(coord1.azimuth.as_f64(), coord2.azimuth.as_f64());
        assert_eq!(coord1.radial_distance, coord2.radial_distance);
    }

    #[test]
    fn test_spherical_coord_different_centers() {
        let coord1 = SphericalCoord::<Barycentric, ICRS>::new_spherical_coord(
            Degrees::new(30.0),
            Degrees::new(60.0),
            100.0
        );
        let coord2 = SphericalCoord::<Geocentric, ICRS>::new_spherical_coord(
            Degrees::new(30.0),
            Degrees::new(60.0),
            100.0
        );
        assert_eq!(coord1.polar.as_f64(), coord2.polar.as_f64());
        assert_eq!(coord1.azimuth.as_f64(), coord2.azimuth.as_f64());
        assert_eq!(coord1.radial_distance, coord2.radial_distance);
    }

    #[test]
    fn test_spherical_coord_different_frames() {
        let coord1 = SphericalCoord::<Geocentric, ICRS>::new_spherical_coord(
            Degrees::new(45.0),
            Degrees::new(90.0),
            150.0
        );
        let coord2 = SphericalCoord::<Geocentric, Ecliptic>::new_spherical_coord(
            Degrees::new(45.0),
            Degrees::new(90.0),
            150.0
        );
        assert_eq!(coord1.polar.as_f64(), coord2.polar.as_f64());
        assert_eq!(coord1.azimuth.as_f64(), coord2.azimuth.as_f64());
        assert_eq!(coord1.radial_distance, coord2.radial_distance);
    }

    #[test]
    fn test_spherical_coord_edge_cases() {
        // Large angles
        let coord = SphericalCoord::<Geocentric, ICRS>::new_spherical_coord(
            Degrees::new(359.999),
            Degrees::new(179.999),
            1e6
        );
        assert!((coord.polar.as_f64() - 359.999).abs() < 1e-6);
        assert!((coord.azimuth.as_f64() - 179.999).abs() < 1e-6);
        assert!((coord.radial_distance - 1e6).abs() < 1e-6);

        // Very small angles
        let coord = SphericalCoord::<Geocentric, ICRS>::new_spherical_coord(
            Degrees::new(0.001),
            Degrees::new(0.001),
            1e-6
        );
        assert!((coord.polar.as_f64() - 0.001).abs() < 1e-6);
        assert!((coord.azimuth.as_f64() - 0.001).abs() < 1e-6);
        assert!((coord.radial_distance - 1e-6).abs() < 1e-6);
    }
}
