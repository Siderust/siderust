//! # Geometric Operations
//!
//! Pure functions for common geometric calculations on 3D vectors.
//! These operate on raw `f64` values without any type system overhead.

/// Calculates the Euclidean magnitude (length) of a 3D vector.
///
/// # Formula
/// ```text
/// |v| = √(x² + y² + z²)
/// ```
#[inline]
pub fn magnitude(x: f64, y: f64, z: f64) -> f64 {
    (x * x + y * y + z * z).sqrt()
}

/// Calculates the squared magnitude of a 3D vector.
///
/// This is more efficient than `magnitude()` when only comparisons are needed.
///
/// # Formula
/// ```text
/// |v|² = x² + y² + z²
/// ```
#[inline]
#[allow(dead_code)]
pub fn magnitude_squared(x: f64, y: f64, z: f64) -> f64 {
    x * x + y * y + z * z
}

/// Normalizes a 3D vector to unit length.
///
/// # Arguments
/// * `x`, `y`, `z` - Vector components
///
/// # Returns
/// Normalized (x, y, z) components. Returns (0, 0, 0) for zero-length input.
#[inline]
pub fn normalize(x: f64, y: f64, z: f64) -> (f64, f64, f64) {
    let mag = magnitude(x, y, z);
    if mag == 0.0 {
        return (0.0, 0.0, 0.0);
    }
    (x / mag, y / mag, z / mag)
}

/// Calculates the Euclidean distance between two points.
///
/// # Arguments
/// * `x1`, `y1`, `z1` - First point
/// * `x2`, `y2`, `z2` - Second point
///
/// # Returns
/// The distance between the two points.
#[inline]
#[allow(dead_code)]
pub fn distance(x1: f64, y1: f64, z1: f64, x2: f64, y2: f64, z2: f64) -> f64 {
    magnitude(x2 - x1, y2 - y1, z2 - z1)
}

/// Calculates the dot product of two 3D vectors.
///
/// # Formula
/// ```text
/// a · b = ax×bx + ay×by + az×bz
/// ```
#[inline]
#[allow(dead_code)]
pub fn dot_product(ax: f64, ay: f64, az: f64, bx: f64, by: f64, bz: f64) -> f64 {
    ax * bx + ay * by + az * bz
}

/// Calculates the cross product of two 3D vectors.
///
/// # Formula
/// ```text
/// a × b = (ay×bz - az×by, az×bx - ax×bz, ax×by - ay×bx)
/// ```
#[inline]
#[allow(dead_code)]
pub fn cross_product(ax: f64, ay: f64, az: f64, bx: f64, by: f64, bz: f64) -> (f64, f64, f64) {
    (
        ay * bz - az * by,
        az * bx - ax * bz,
        ax * by - ay * bx,
    )
}

/// Calculates the angular separation between two directions using the
/// robust vincenty formula (generalized haversine).
///
/// This is numerically stable for all angular separations, including
/// very small angles and angles near 180°.
///
/// # Arguments
/// * `az1_rad`, `po1_rad` - First direction (azimuth, polar) in radians
/// * `az2_rad`, `po2_rad` - Second direction (azimuth, polar) in radians
///
/// # Returns
/// Angular separation in radians.
///
/// # Formula
/// Uses the Vincenty formula for angular distance on a sphere:
/// ```text
/// Δσ = atan2(
///     √[(cos(θ₂)×sin(Δφ))² + (cos(θ₁)×sin(θ₂) - sin(θ₁)×cos(θ₂)×cos(Δφ))²],
///     sin(θ₁)×sin(θ₂) + cos(θ₁)×cos(θ₂)×cos(Δφ)
/// )
/// ```
/// where θ is polar angle and φ is azimuth.
#[inline]
pub fn angular_separation(az1_rad: f64, po1_rad: f64, az2_rad: f64, po2_rad: f64) -> f64 {
    let (sin_po1, cos_po1) = po1_rad.sin_cos();
    let (sin_po2, cos_po2) = po2_rad.sin_cos();
    let (sin_daz, cos_daz) = (az2_rad - az1_rad).sin_cos();

    // Vincenty formula numerator components
    let x = cos_po1 * sin_po2 - sin_po1 * cos_po2 * cos_daz;
    let y = cos_po2 * sin_daz;

    // Vincenty formula denominator
    let z = sin_po1 * sin_po2 + cos_po1 * cos_po2 * cos_daz;

    (x * x + y * y).sqrt().atan2(z)
}

/// Calculates the angular separation between two unit vectors in Cartesian form.
///
/// # Arguments
/// * `x1`, `y1`, `z1` - First unit vector
/// * `x2`, `y2`, `z2` - Second unit vector
///
/// # Returns
/// Angular separation in radians.
#[inline]
#[allow(dead_code)]
pub fn angular_separation_cartesian(
    x1: f64,
    y1: f64,
    z1: f64,
    x2: f64,
    y2: f64,
    z2: f64,
) -> f64 {
    // For unit vectors, dot product gives cos(angle)
    // Use atan2 for numerical stability
    let dot = dot_product(x1, y1, z1, x2, y2, z2);
    let (cx, cy, cz) = cross_product(x1, y1, z1, x2, y2, z2);
    let cross_mag = magnitude(cx, cy, cz);
    cross_mag.atan2(dot)
}

/// Subtracts two 3D vectors.
///
/// # Returns
/// (a - b) as (dx, dy, dz)
#[inline]
#[allow(dead_code)]
pub fn subtract(ax: f64, ay: f64, az: f64, bx: f64, by: f64, bz: f64) -> (f64, f64, f64) {
    (ax - bx, ay - by, az - bz)
}

/// Adds two 3D vectors.
///
/// # Returns
/// (a + b) as (x, y, z)
#[inline]
#[allow(dead_code)]
pub fn add(ax: f64, ay: f64, az: f64, bx: f64, by: f64, bz: f64) -> (f64, f64, f64) {
    (ax + bx, ay + by, az + bz)
}

/// Scales a 3D vector by a scalar.
#[inline]
#[allow(dead_code)]
pub fn scale(x: f64, y: f64, z: f64, s: f64) -> (f64, f64, f64) {
    (x * s, y * s, z * s)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::{FRAC_PI_2, PI};

    const EPSILON: f64 = 1e-12;

    fn assert_approx_eq(a: f64, b: f64, msg: &str) {
        assert!(
            (a - b).abs() < EPSILON,
            "{msg}: {a} vs {b}, diff={}",
            (a - b).abs()
        );
    }

    #[test]
    fn test_magnitude() {
        assert_approx_eq(magnitude(3.0, 4.0, 0.0), 5.0, "3-4-5 triangle");
        assert_approx_eq(magnitude(1.0, 1.0, 1.0), 3.0_f64.sqrt(), "unit cube diagonal");
        assert_approx_eq(magnitude(0.0, 0.0, 0.0), 0.0, "zero vector");
    }

    #[test]
    fn test_normalize() {
        let (x, y, z) = normalize(3.0, 0.0, 4.0);
        assert_approx_eq(magnitude(x, y, z), 1.0, "normalized magnitude");
        assert_approx_eq(x, 0.6, "x component");
        assert_approx_eq(z, 0.8, "z component");
    }

    #[test]
    fn test_normalize_zero() {
        let (x, y, z) = normalize(0.0, 0.0, 0.0);
        assert_eq!((x, y, z), (0.0, 0.0, 0.0), "zero vector normalization");
    }

    #[test]
    fn test_distance() {
        assert_approx_eq(
            distance(0.0, 0.0, 0.0, 1.0, 0.0, 0.0),
            1.0,
            "unit distance along x"
        );
        assert_approx_eq(
            distance(1.0, 2.0, 3.0, 1.0, 2.0, 3.0),
            0.0,
            "same point distance"
        );
    }

    #[test]
    fn test_dot_product() {
        assert_approx_eq(
            dot_product(1.0, 0.0, 0.0, 0.0, 1.0, 0.0),
            0.0,
            "perpendicular vectors"
        );
        assert_approx_eq(
            dot_product(1.0, 0.0, 0.0, 1.0, 0.0, 0.0),
            1.0,
            "parallel vectors"
        );
        assert_approx_eq(
            dot_product(1.0, 2.0, 3.0, 4.0, 5.0, 6.0),
            32.0,
            "general case"
        );
    }

    #[test]
    fn test_cross_product() {
        // x × y = z
        let (cx, cy, cz) = cross_product(1.0, 0.0, 0.0, 0.0, 1.0, 0.0);
        assert_approx_eq(cx, 0.0, "x component");
        assert_approx_eq(cy, 0.0, "y component");
        assert_approx_eq(cz, 1.0, "z component");
    }

    #[test]
    fn test_angular_separation_same_direction() {
        let sep = angular_separation(0.0, 0.0, 0.0, 0.0);
        assert_approx_eq(sep, 0.0, "same direction");
    }

    #[test]
    fn test_angular_separation_90_degrees() {
        // Along x-axis (az=0, po=0) vs along z-axis (az=0, po=π/2)
        let sep = angular_separation(0.0, 0.0, 0.0, FRAC_PI_2);
        assert_approx_eq(sep, FRAC_PI_2, "90 degree separation");
    }

    #[test]
    fn test_angular_separation_180_degrees() {
        // Opposite directions on equator
        let sep = angular_separation(0.0, 0.0, PI, 0.0);
        assert_approx_eq(sep, PI, "180 degree separation");
    }

    #[test]
    fn test_angular_separation_cartesian() {
        // Unit vectors along x and y axes
        let sep = angular_separation_cartesian(1.0, 0.0, 0.0, 0.0, 1.0, 0.0);
        assert_approx_eq(sep, FRAC_PI_2, "90 degrees x to y");
    }
}
