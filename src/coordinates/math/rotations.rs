//! # Frame Rotation Matrices
//!
//! Pure rotation functions for transforming 3D vectors between reference frames.
//! These functions implement the mathematical rotations without any type system
//! overhead.
//!
//! ## Obliquity of the Ecliptic
//!
//! The mean obliquity of the ecliptic at J2000.0 is approximately 23.439281°.
//! This is the angle between Earth's equatorial plane and the ecliptic plane.
//!
//! ## Rotation Convention
//!
//! All rotations are right-handed about the specified axis:
//! - Ecliptic ↔ Equatorial: Rotation about the X-axis by obliquity ε

/// Mean obliquity of the ecliptic at J2000.0 in radians.
///
/// ε = 23.439281° converted to radians
/// Computed as: 23.439281_f64.to_radians()
pub const OBLIQUITY_J2000_RAD: f64 = 0.40909262775014904;

/// Mean obliquity of the ecliptic at J2000.0 in degrees.
#[allow(dead_code)]
pub const OBLIQUITY_J2000_DEG: f64 = 23.439281;

/// Result of a 3D rotation operation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RotationResult {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl RotationResult {
    /// Creates a new rotation result.
    #[inline]
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }
}

/// Rotates a vector from the Ecliptic frame to the Equatorial frame.
///
/// This is a right-hand rotation about the +X axis by the obliquity ε:
///
/// ```text
/// ┌    ┐   ┌ 1    0        0     ┐ ┌   ┐
/// │ x' │   │ 0  cos(ε)  -sin(ε)  │ │ x │
/// │ y' │ = │ 0  sin(ε)   cos(ε)  │ │ y │
/// │ z' │   └                     ┘ │ z │
/// └    ┘                           └   ┘
/// ```
///
/// # Arguments
/// * `x`, `y`, `z` - Input vector components in the Ecliptic frame
///
/// # Returns
/// The rotated vector components in the Equatorial frame
#[inline]
pub fn rotate_ecliptic_to_equatorial(x: f64, y: f64, z: f64) -> RotationResult {
    let (sin_e, cos_e) = OBLIQUITY_J2000_RAD.sin_cos();
    RotationResult::new(x, cos_e * y - sin_e * z, sin_e * y + cos_e * z)
}

/// Rotates a vector from the Equatorial frame to the Ecliptic frame.
///
/// This is the inverse rotation (transpose of the forward matrix):
///
/// ```text
/// ┌    ┐   ┌ 1    0        0     ┐ ┌   ┐
/// │ x' │   │ 0  cos(ε)   sin(ε)  │ │ x │
/// │ y' │ = │ 0 -sin(ε)   cos(ε)  │ │ y │
/// │ z' │   └                     ┘ │ z │
/// └    ┘                           └   ┘
/// ```
///
/// # Arguments
/// * `x`, `y`, `z` - Input vector components in the Equatorial frame
///
/// # Returns
/// The rotated vector components in the Ecliptic frame
#[inline]
pub fn rotate_equatorial_to_ecliptic(x: f64, y: f64, z: f64) -> RotationResult {
    let (sin_e, cos_e) = OBLIQUITY_J2000_RAD.sin_cos();
    RotationResult::new(x, cos_e * y + sin_e * z, -sin_e * y + cos_e * z)
}

/// Rotates a vector about the X-axis by the given angle (in radians).
///
/// This is a general rotation about +X:
///
/// ```text
/// ┌    ┐   ┌ 1    0         0      ┐ ┌   ┐
/// │ x' │   │ 0  cos(θ)  -sin(θ)    │ │ x │
/// │ y' │ = │ 0  sin(θ)   cos(θ)    │ │ y │
/// │ z' │   └                       ┘ │ z │
/// └    ┘                             └   ┘
/// ```
#[inline]
#[allow(dead_code)]
pub fn rotate_about_x(x: f64, y: f64, z: f64, angle_rad: f64) -> RotationResult {
    let (sin_a, cos_a) = angle_rad.sin_cos();
    RotationResult::new(x, cos_a * y - sin_a * z, sin_a * y + cos_a * z)
}

/// Rotates a vector about the Y-axis by the given angle (in radians).
///
/// ```text
/// ┌    ┐   ┌  cos(θ)  0  sin(θ)  ┐ ┌   ┐
/// │ x' │   │    0     1    0     │ │ x │
/// │ y' │ = │ -sin(θ)  0  cos(θ)  │ │ y │
/// │ z' │   └                     ┘ │ z │
/// └    ┘                           └   ┘
/// ```
#[inline]
#[allow(dead_code)]
pub fn rotate_about_y(x: f64, y: f64, z: f64, angle_rad: f64) -> RotationResult {
    let (sin_a, cos_a) = angle_rad.sin_cos();
    RotationResult::new(cos_a * x + sin_a * z, y, -sin_a * x + cos_a * z)
}

/// Rotates a vector about the Z-axis by the given angle (in radians).
///
/// ```text
/// ┌    ┐   ┌ cos(θ)  -sin(θ)  0  ┐ ┌   ┐
/// │ x' │   │ sin(θ)   cos(θ)  0  │ │ x │
/// │ y' │ = │   0        0     1  │ │ y │
/// │ z' │   └                     ┘ │ z │
/// └    ┘                           └   ┘
/// ```
#[inline]
#[allow(dead_code)]
pub fn rotate_about_z(x: f64, y: f64, z: f64, angle_rad: f64) -> RotationResult {
    let (sin_a, cos_a) = angle_rad.sin_cos();
    RotationResult::new(cos_a * x - sin_a * y, sin_a * x + cos_a * y, z)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::FRAC_PI_2;

    const EPSILON: f64 = 1e-12;

    fn assert_approx_eq(result: RotationResult, expected: (f64, f64, f64), msg: &str) {
        assert!(
            (result.x - expected.0).abs() < EPSILON,
            "{msg}: x mismatch: {} vs {}",
            result.x,
            expected.0
        );
        assert!(
            (result.y - expected.1).abs() < EPSILON,
            "{msg}: y mismatch: {} vs {}",
            result.y,
            expected.1
        );
        assert!(
            (result.z - expected.2).abs() < EPSILON,
            "{msg}: z mismatch: {} vs {}",
            result.z,
            expected.2
        );
    }

    #[test]
    fn test_ecliptic_equatorial_round_trip() {
        let original = (1.0, 2.0, 3.0);

        // Forward transform
        let eq = rotate_ecliptic_to_equatorial(original.0, original.1, original.2);

        // Inverse transform
        let back = rotate_equatorial_to_ecliptic(eq.x, eq.y, eq.z);

        assert_approx_eq(back, original, "round trip ecliptic->equatorial->ecliptic");
    }

    #[test]
    fn test_x_axis_unchanged_by_obliquity_rotation() {
        // X-axis is the rotation axis, should be unchanged
        let result = rotate_ecliptic_to_equatorial(1.0, 0.0, 0.0);
        assert_approx_eq(result, (1.0, 0.0, 0.0), "X-axis should be unchanged");
    }

    #[test]
    fn test_rotate_about_z_90_degrees() {
        // Rotate (1, 0, 0) by 90° about Z should give (0, 1, 0)
        let result = rotate_about_z(1.0, 0.0, 0.0, FRAC_PI_2);
        assert_approx_eq(result, (0.0, 1.0, 0.0), "90° rotation about Z");
    }

    #[test]
    fn test_rotate_about_y_90_degrees() {
        // Rotate (1, 0, 0) by 90° about Y should give (0, 0, -1)
        let result = rotate_about_y(1.0, 0.0, 0.0, FRAC_PI_2);
        assert_approx_eq(result, (0.0, 0.0, -1.0), "90° rotation about Y");
    }

    #[test]
    fn test_rotate_about_x_90_degrees() {
        // Rotate (0, 1, 0) by 90° about X should give (0, 0, 1)
        let result = rotate_about_x(0.0, 1.0, 0.0, FRAC_PI_2);
        assert_approx_eq(result, (0.0, 0.0, 1.0), "90° rotation about X");
    }
}
