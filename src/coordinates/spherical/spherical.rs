//! # Spherical Coordinates (Positions)
//!
//! This module defines [`SphericalCoord<C, F, U>`], the core spherical **position** type: it carries
//! a reference **center** (`C`), a reference **frame** (`F`), and a radial **distance** (`U`).
//!
//! Spherical **directions** (unit vectors) are represented separately by [`crate::coordinates::spherical::Direction<F>`]
//! and intentionally have **no** reference center.
//!
//! ## Overview
//! Legacy module removed.
//!
//! The crate used to define `SphericalCoord` here. It has been removed in favor of:
//! - `crate::coordinates::spherical::Position<C, F, U>` for spherical **positions** (center + frame + distance)
//! - `crate::coordinates::spherical::Direction<F>` for spherical **directions** (frame-only)
//!
//! This file remains only to provide a clear error if someone tries to
//! re-introduce `mod spherical; pub use spherical::*;`.

compile_error!(
    "`SphericalCoord` has been removed. Use `crate::coordinates::spherical::Position` (positions) or `crate::coordinates::spherical::Direction` (directions)."
);

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
