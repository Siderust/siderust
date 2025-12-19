//! # Coordinate System Conversions
//!
//! Pure functions for converting between spherical and Cartesian coordinate
//! representations. All angles are in radians for internal calculations.
//!
//! ## Conventions
//!
//! - **Spherical coordinates**: (azimuth, polar, radius)
//!   - `azimuth` (φ): angle in the XY plane from +X toward +Y (longitude, RA)
//!   - `polar` (θ): angle from the XY plane toward +Z (latitude, declination)
//!   - `radius` (r): distance from origin
//!
//! - **Cartesian coordinates**: (x, y, z)
//!   - Standard right-handed coordinate system

/// Result of a spherical-to-Cartesian conversion.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CartesianResult {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl CartesianResult {
    /// Creates a new Cartesian result.
    #[inline]
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }
}

/// Result of a Cartesian-to-spherical conversion.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SphericalResult {
    /// Azimuthal angle in radians (φ): angle in XY plane from +X toward +Y
    pub azimuth_rad: f64,
    /// Polar angle in radians (θ): angle from XY plane toward +Z
    pub polar_rad: f64,
    /// Radial distance from origin
    pub radius: f64,
}

impl SphericalResult {
    /// Creates a new spherical result.
    #[inline]
    pub const fn new(azimuth_rad: f64, polar_rad: f64, radius: f64) -> Self {
        Self {
            azimuth_rad,
            polar_rad,
            radius,
        }
    }

    /// Returns azimuth in degrees.
    #[inline]
    pub fn azimuth_deg(&self) -> f64 {
        self.azimuth_rad.to_degrees()
    }

    /// Returns polar angle in degrees.
    #[inline]
    pub fn polar_deg(&self) -> f64 {
        self.polar_rad.to_degrees()
    }
}

/// Converts spherical coordinates to Cartesian coordinates.
///
/// # Arguments
/// * `azimuth_rad` - Azimuthal angle in radians (φ): angle in XY plane from +X
/// * `polar_rad` - Polar angle in radians (θ): angle from XY plane
/// * `radius` - Radial distance from origin
///
/// # Formulas
/// ```text
/// x = r × cos(θ) × cos(φ)
/// y = r × cos(θ) × sin(φ)
/// z = r × sin(θ)
/// ```
///
/// # Returns
/// Cartesian (x, y, z) coordinates
#[inline]
pub fn spherical_to_cartesian(azimuth_rad: f64, polar_rad: f64, radius: f64) -> CartesianResult {
    let (sin_polar, cos_polar) = polar_rad.sin_cos();
    let (sin_azimuth, cos_azimuth) = azimuth_rad.sin_cos();

    CartesianResult::new(
        radius * cos_polar * cos_azimuth,
        radius * cos_polar * sin_azimuth,
        radius * sin_polar,
    )
}

/// Converts Cartesian coordinates to spherical coordinates.
///
/// # Arguments
/// * `x`, `y`, `z` - Cartesian coordinates
///
/// # Formulas
/// ```text
/// r = √(x² + y² + z²)
/// θ = arcsin(z / r)      (polar angle from XY plane)
/// φ = atan2(y, x)        (azimuth angle in XY plane)
/// ```
///
/// # Returns
/// Spherical coordinates with angles in radians.
/// Returns (0, 0, 0) for a zero-length input vector.
#[inline]
pub fn cartesian_to_spherical(x: f64, y: f64, z: f64) -> SphericalResult {
    let radius = (x * x + y * y + z * z).sqrt();

    if radius == 0.0 {
        return SphericalResult::new(0.0, 0.0, 0.0);
    }

    let polar_rad = (z / radius).asin();
    let azimuth_rad = y.atan2(x);

    SphericalResult::new(azimuth_rad, polar_rad, radius)
}

/// Converts Cartesian direction components to spherical angles.
///
/// This is a convenience function for unit vectors where radius is implicitly 1.
///
/// # Arguments
/// * `x`, `y`, `z` - Cartesian direction components (should be normalized)
///
/// # Returns
/// (azimuth_rad, polar_rad) angles in radians
#[inline]
#[allow(dead_code)]
pub fn direction_to_angles(x: f64, y: f64, z: f64) -> (f64, f64) {
    let result = cartesian_to_spherical(x, y, z);
    (result.azimuth_rad, result.polar_rad)
}

/// Converts spherical angles to Cartesian direction components.
///
/// This is a convenience function for unit vectors where radius is implicitly 1.
///
/// # Arguments
/// * `azimuth_rad` - Azimuthal angle in radians
/// * `polar_rad` - Polar angle in radians
///
/// # Returns
/// (x, y, z) direction components (unit vector)
#[inline]
#[allow(dead_code)]
pub fn angles_to_direction(azimuth_rad: f64, polar_rad: f64) -> CartesianResult {
    spherical_to_cartesian(azimuth_rad, polar_rad, 1.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI};

    const EPSILON: f64 = 1e-12;

    fn assert_approx_eq(a: f64, b: f64, msg: &str) {
        assert!(
            (a - b).abs() < EPSILON,
            "{msg}: {a} vs {b}, diff={}",
            (a - b).abs()
        );
    }

    #[test]
    fn test_spherical_to_cartesian_x_axis() {
        // (azimuth=0, polar=0, r=1) -> (1, 0, 0)
        let result = spherical_to_cartesian(0.0, 0.0, 1.0);
        assert_approx_eq(result.x, 1.0, "x");
        assert_approx_eq(result.y, 0.0, "y");
        assert_approx_eq(result.z, 0.0, "z");
    }

    #[test]
    fn test_spherical_to_cartesian_y_axis() {
        // (azimuth=π/2, polar=0, r=1) -> (0, 1, 0)
        let result = spherical_to_cartesian(FRAC_PI_2, 0.0, 1.0);
        assert_approx_eq(result.x, 0.0, "x");
        assert_approx_eq(result.y, 1.0, "y");
        assert_approx_eq(result.z, 0.0, "z");
    }

    #[test]
    fn test_spherical_to_cartesian_z_axis() {
        // (azimuth=0, polar=π/2, r=1) -> (0, 0, 1)
        let result = spherical_to_cartesian(0.0, FRAC_PI_2, 1.0);
        assert_approx_eq(result.x, 0.0, "x");
        assert_approx_eq(result.y, 0.0, "y");
        assert_approx_eq(result.z, 1.0, "z");
    }

    #[test]
    fn test_cartesian_to_spherical_round_trip() {
        let original = (1.0, 2.0, 3.0);
        let sph = cartesian_to_spherical(original.0, original.1, original.2);
        let cart = spherical_to_cartesian(sph.azimuth_rad, sph.polar_rad, sph.radius);

        assert_approx_eq(cart.x, original.0, "x round trip");
        assert_approx_eq(cart.y, original.1, "y round trip");
        assert_approx_eq(cart.z, original.2, "z round trip");
    }

    #[test]
    fn test_zero_vector() {
        let result = cartesian_to_spherical(0.0, 0.0, 0.0);
        assert_eq!(result.radius, 0.0, "zero vector should have zero radius");
    }

    #[test]
    fn test_45_degree_angles() {
        // 45° in XY plane, 45° elevation
        let azimuth = FRAC_PI_4;
        let polar = FRAC_PI_4;
        let radius = 2.0_f64.sqrt();

        let cart = spherical_to_cartesian(azimuth, polar, radius);
        let back = cartesian_to_spherical(cart.x, cart.y, cart.z);

        assert_approx_eq(back.azimuth_rad, azimuth, "azimuth round trip");
        assert_approx_eq(back.polar_rad, polar, "polar round trip");
        assert_approx_eq(back.radius, radius, "radius round trip");
    }

    #[test]
    fn test_negative_x_axis() {
        // (azimuth=π, polar=0, r=1) -> (-1, 0, 0)
        let result = spherical_to_cartesian(PI, 0.0, 1.0);
        assert_approx_eq(result.x, -1.0, "x");
        assert_approx_eq(result.y, 0.0, "y");
        assert_approx_eq(result.z, 0.0, "z");
    }
}
