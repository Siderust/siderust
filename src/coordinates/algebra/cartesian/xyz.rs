//! # XYZ: Shared Storage for Cartesian Types
//!
//! This module provides the internal `XYZ<T>` type, a thin `#[repr(transparent)]`
//! wrapper around `nalgebra::Vector3<T>`.
//!
//! ## Purpose
//!
//! `XYZ<T>` centralizes all component-wise mathematical operations in one place:
//! - Addition, subtraction, scaling
//! - Dot product, magnitude, normalization
//!
//! The semantic types (`Position`, `Vector`, `Direction`) are thin wrappers over
//! `XYZ<T>` with `PhantomData` markers for center, frame, and unit information.
//! This eliminates code duplication while maintaining zero-cost abstractions.
//!
//! ## Design Notes
//!
//! - Uses `#[repr(transparent)]` for zero-cost ABI compatibility with `Vector3<T>`
//! - All operations are `#[inline]` for optimal performance
//! - Generic over the component type `T` (typically `f64` or `Quantity<U>`)

use nalgebra::Vector3;
use std::ops::{Add, Mul, Neg, Sub};

/// Internal shared storage for 3D Cartesian coordinates.
///
/// This is a thin wrapper around `nalgebra::Vector3<T>` that provides
/// centralized implementations of component-wise operations.
///
/// # Type Parameter
/// - `T`: The component type (e.g., `f64`, `Quantity<U>`)
#[derive(Debug, Clone, Copy, PartialEq)]
#[repr(transparent)]
pub struct XYZ<T>(pub(crate) Vector3<T>);

// =============================================================================
// Constructors
// =============================================================================

impl<T> XYZ<T> {
    /// Creates a new XYZ from individual components.
    #[inline]
    pub const fn new(x: T, y: T, z: T) -> Self {
        Self(Vector3::new(x, y, z))
    }

    /// Creates an XYZ from a nalgebra Vector3.
    #[inline]
    pub const fn from_vec3(vec: Vector3<T>) -> Self {
        Self(vec)
    }

    /// Returns the underlying nalgebra Vector3.
    #[inline]
    pub const fn as_vec3(&self) -> &Vector3<T> {
        &self.0
    }

    /// Consumes self and returns the underlying Vector3.
    #[inline]
    pub fn into_vec3(self) -> Vector3<T> {
        self.0
    }
}

// =============================================================================
// Component Access
// =============================================================================

impl<T: Copy> XYZ<T> {
    /// Returns the x-component.
    #[inline]
    pub fn x(&self) -> T {
        self.0[0]
    }

    /// Returns the y-component.
    #[inline]
    pub fn y(&self) -> T {
        self.0[1]
    }

    /// Returns the z-component.
    #[inline]
    pub fn z(&self) -> T {
        self.0[2]
    }

    /// Returns components as a tuple.
    #[inline]
    pub fn components(&self) -> (T, T, T) {
        (self.x(), self.y(), self.z())
    }
}

// =============================================================================
// Arithmetic Operations for f64 components
// =============================================================================

impl XYZ<f64> {
    /// The zero vector.
    pub const ZERO: Self = Self::new(0.0, 0.0, 0.0);

    /// Computes the Euclidean magnitude (length) of the vector.
    #[inline]
    pub fn magnitude(&self) -> f64 {
        self.0.magnitude()
    }

    /// Computes the squared magnitude (avoids sqrt).
    #[inline]
    pub fn magnitude_squared(&self) -> f64 {
        self.0.magnitude_squared()
    }

    /// Computes the dot product with another vector.
    #[inline]
    pub fn dot(&self, other: &Self) -> f64 {
        self.0.dot(&other.0)
    }

    /// Computes the cross product with another vector.
    #[inline]
    pub fn cross(&self, other: &Self) -> Self {
        Self(self.0.cross(&other.0))
    }

    /// Returns a normalized (unit length) version, or None if magnitude is zero.
    #[inline]
    pub fn try_normalize(&self) -> Option<Self> {
        let mag = self.magnitude();
        if mag.abs() < f64::EPSILON {
            None
        } else {
            Some(Self(self.0 / mag))
        }
    }

    /// Returns a normalized version, assuming non-zero magnitude.
    ///
    /// # Safety
    /// Caller must ensure the vector is non-zero.
    #[inline]
    pub fn normalize_unchecked(&self) -> Self {
        Self(self.0.normalize())
    }

    /// Scales the vector by a scalar factor.
    #[inline]
    pub fn scale(&self, factor: f64) -> Self {
        Self(self.0 * factor)
    }

    /// Component-wise addition.
    #[inline]
    pub fn add(&self, other: &Self) -> Self {
        Self(self.0 + other.0)
    }

    /// Component-wise subtraction.
    #[inline]
    pub fn sub(&self, other: &Self) -> Self {
        Self(self.0 - other.0)
    }

    /// Negates all components.
    #[inline]
    pub fn neg(&self) -> Self {
        Self(-self.0)
    }
}

// =============================================================================
// Operator Implementations for XYZ<f64>
// =============================================================================

impl Add for XYZ<f64> {
    type Output = Self;

    #[inline]
    fn add(self, other: Self) -> Self::Output {
        Self(self.0 + other.0)
    }
}

impl Sub for XYZ<f64> {
    type Output = Self;

    #[inline]
    fn sub(self, other: Self) -> Self::Output {
        Self(self.0 - other.0)
    }
}

impl Neg for XYZ<f64> {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self::Output {
        Self(-self.0)
    }
}

impl Mul<f64> for XYZ<f64> {
    type Output = Self;

    #[inline]
    fn mul(self, scalar: f64) -> Self::Output {
        Self(self.0 * scalar)
    }
}

impl Mul<XYZ<f64>> for f64 {
    type Output = XYZ<f64>;

    #[inline]
    fn mul(self, xyz: XYZ<f64>) -> Self::Output {
        XYZ(xyz.0 * self)
    }
}

// =============================================================================
// Quantity-aware operations
// =============================================================================

use qtty::{Quantity, Unit};

impl<U: Unit> XYZ<Quantity<U>> {
    /// The zero vector for this unit type.
    pub const ZERO: Self = Self::new(
        Quantity::<U>::new(0.0),
        Quantity::<U>::new(0.0),
        Quantity::<U>::new(0.0),
    );

    /// Computes the Euclidean magnitude in the same unit.
    #[inline]
    pub fn magnitude(&self) -> Quantity<U> {
        let mag = Vector3::new(
            self.0[0].value(),
            self.0[1].value(),
            self.0[2].value(),
        )
        .magnitude();
        Quantity::new(mag)
    }

    /// Computes the squared magnitude.
    #[inline]
    pub fn magnitude_squared(&self) -> f64 {
        let x = self.0[0].value();
        let y = self.0[1].value();
        let z = self.0[2].value();
        x * x + y * y + z * z
    }

    /// Extracts raw f64 values as an XYZ<f64>.
    #[inline]
    pub fn to_raw(&self) -> XYZ<f64> {
        XYZ::new(
            self.0[0].value(),
            self.0[1].value(),
            self.0[2].value(),
        )
    }

    /// Creates from raw f64 values.
    #[inline]
    pub fn from_raw(raw: XYZ<f64>) -> Self {
        Self::new(
            Quantity::new(raw.x()),
            Quantity::new(raw.y()),
            Quantity::new(raw.z()),
        )
    }

    /// Component-wise addition.
    #[inline]
    pub fn add(&self, other: &Self) -> Self {
        Self(self.0 + other.0)
    }

    /// Component-wise subtraction.
    #[inline]
    pub fn sub(&self, other: &Self) -> Self {
        Self(self.0 - other.0)
    }
}

impl<U: Unit> Add for XYZ<Quantity<U>> {
    type Output = Self;

    #[inline]
    fn add(self, other: Self) -> Self::Output {
        Self(self.0 + other.0)
    }
}

impl<U: Unit> Sub for XYZ<Quantity<U>> {
    type Output = Self;

    #[inline]
    fn sub(self, other: Self) -> Self::Output {
        Self(self.0 - other.0)
    }
}

impl<U: Unit> Neg for XYZ<Quantity<U>> {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self::Output {
        Self::new(-self.0[0], -self.0[1], -self.0[2])
    }
}

// =============================================================================
// Default Implementation
// =============================================================================

impl<T: Default> Default for XYZ<T> {
    fn default() -> Self {
        Self(Vector3::new(T::default(), T::default(), T::default()))
    }
}

// =============================================================================
// Display
// =============================================================================

impl<T: std::fmt::Display + Copy> std::fmt::Display for XYZ<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "({:.6}, {:.6}, {:.6})", self.x(), self.y(), self.z())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_xyz_basic_ops() {
        let a = XYZ::new(1.0, 2.0, 3.0);
        let b = XYZ::new(4.0, 5.0, 6.0);

        // Addition
        let sum = a + b;
        assert!((sum.x() - 5.0).abs() < f64::EPSILON);
        assert!((sum.y() - 7.0).abs() < f64::EPSILON);
        assert!((sum.z() - 9.0).abs() < f64::EPSILON);

        // Subtraction
        let diff = b - a;
        assert!((diff.x() - 3.0).abs() < f64::EPSILON);
        assert!((diff.y() - 3.0).abs() < f64::EPSILON);
        assert!((diff.z() - 3.0).abs() < f64::EPSILON);

        // Scaling
        let scaled = a.scale(2.0);
        assert!((scaled.x() - 2.0).abs() < f64::EPSILON);
        assert!((scaled.y() - 4.0).abs() < f64::EPSILON);
        assert!((scaled.z() - 6.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_xyz_magnitude_and_normalize() {
        let v = XYZ::new(3.0, 4.0, 0.0);
        assert!((v.magnitude() - 5.0).abs() < f64::EPSILON);

        let n = v.try_normalize().unwrap();
        assert!((n.magnitude() - 1.0).abs() < f64::EPSILON);
        assert!((n.x() - 0.6).abs() < f64::EPSILON);
        assert!((n.y() - 0.8).abs() < f64::EPSILON);
    }

    #[test]
    fn test_xyz_dot_and_cross() {
        let a = XYZ::new(1.0, 0.0, 0.0);
        let b = XYZ::new(0.0, 1.0, 0.0);

        // Dot product of perpendicular vectors
        assert!(a.dot(&b).abs() < f64::EPSILON);

        // Cross product
        let c = a.cross(&b);
        assert!((c.x()).abs() < f64::EPSILON);
        assert!((c.y()).abs() < f64::EPSILON);
        assert!((c.z() - 1.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_xyz_zero_normalize() {
        let zero = XYZ::<f64>::ZERO;
        assert!(zero.try_normalize().is_none());
    }
}
