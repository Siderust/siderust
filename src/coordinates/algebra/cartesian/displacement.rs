//! # Cartesian Displacement Vector
//!
//! This module defines [`Displacement<F, U>`], a **free displacement vector** in 3D space.
//!
//! ## Mathematical Model
//!
//! A displacement vector represents a directed magnitude in space. It is:
//! - **Frame-dependent**: The orientation is relative to a reference frame `F`
//! - **Center-independent**: Displacement is translation-invariant
//! - **Dimensioned**: Has a length unit `U`
//!
//! Displacements form a vector space and support:
//! - Addition: `Displacement + Displacement -> Displacement`
//! - Subtraction: `Displacement - Displacement -> Displacement`
//! - Scalar multiplication: `Displacement * f64 -> Displacement`
//! - Normalization: `normalize(Displacement) -> Option<Direction>`
//!
//! ## Relationship to Other Types
//!
//! - `Position - Position -> Displacement` (displacement between points)
//! - `Position + Displacement -> Position` (translating a point)
//! - `Direction * Length -> Displacement` (scaling a direction)
//! - `normalize(Displacement) -> Direction` (extracting orientation)
//!
//! ## Example
//!
//! ```rust
//! use siderust::coordinates::algebra::cartesian::Displacement;
//! use siderust::coordinates::algebra::frames::Ecliptic;
//! use qtty::*;
//!
//! // Create a displacement vector
//! let displacement = Displacement::<Ecliptic, AstronomicalUnit>::new(1.0, 2.0, 3.0);
//!
//! // Displacement arithmetic
//! let doubled = displacement + displacement;
//! assert!((doubled.x().value() - 2.0).abs() < 1e-12);
//!
//! // Get the magnitude
//! let magnitude = displacement.magnitude();
//! ```

use super::xyz::XYZ;
use crate::coordinates::algebra::frames::ReferenceFrame;
use qtty::{LengthUnit, Quantity};

use std::marker::PhantomData;
use std::ops::{Add, Neg, Sub};

/// A free displacement vector in 3D Cartesian coordinates.
///
/// Displacements are frame-dependent but center-independent. They represent
/// directed magnitudes (displacements) in space.
///
/// # Type Parameters
/// - `F`: The reference frame (e.g., `ICRS`, `Ecliptic`, `Equatorial`)
/// - `U`: The length unit (e.g., `AstronomicalUnit`, `Kilometer`)
///
/// # Zero-Cost Abstraction
///
/// This type uses `#[repr(transparent)]` over the internal storage,
/// ensuring no runtime overhead compared to raw `Vector3<Quantity<U>>`.
#[derive(Debug, Clone, Copy)]
#[repr(transparent)]
pub struct Displacement<F: ReferenceFrame, U: LengthUnit> {
    storage: DisplacementStorage<F, U>,
}

/// Internal storage for Displacement with PhantomData marker.
#[derive(Debug, Clone, Copy)]
struct DisplacementStorage<F: ReferenceFrame, U: LengthUnit> {
    xyz: XYZ<Quantity<U>>,
    _frame: PhantomData<F>,
}

// =============================================================================
// Constructors
// =============================================================================

impl<F: ReferenceFrame, U: LengthUnit> Displacement<F, U> {
    /// Creates a new displacement vector from component values.
    ///
    /// # Arguments
    /// - `x`, `y`, `z`: Component values (converted to `Quantity<U>`)
    #[inline]
    pub fn new<T: Into<Quantity<U>>>(x: T, y: T, z: T) -> Self {
        Self {
            storage: DisplacementStorage {
                xyz: XYZ::new(x.into(), y.into(), z.into()),
                _frame: PhantomData,
            },
        }
    }

    /// Creates a vector from a nalgebra Vector3.
    #[inline]
    pub fn from_vec3(vec3: nalgebra::Vector3<Quantity<U>>) -> Self {
        Self {
            storage: DisplacementStorage {
                xyz: XYZ::from_vec3(vec3),
                _frame: PhantomData,
            },
        }
    }

    /// Creates a vector from the internal XYZ storage.
    #[inline]
    pub(crate) fn from_xyz(xyz: XYZ<Quantity<U>>) -> Self {
        Self {
            storage: DisplacementStorage {
                xyz,
                _frame: PhantomData,
            },
        }
    }

    /// The zero vector.
    pub const ZERO: Self = Self {
        storage: DisplacementStorage {
            xyz: XYZ::new(
                Quantity::<U>::new(0.0),
                Quantity::<U>::new(0.0),
                Quantity::<U>::new(0.0),
            ),
            _frame: PhantomData,
        },
    };
}

// =============================================================================
// Component Access
// =============================================================================

impl<F: ReferenceFrame, U: LengthUnit> Displacement<F, U> {
    /// Returns the x-component.
    #[inline]
    pub fn x(&self) -> Quantity<U> {
        self.storage.xyz.x()
    }

    /// Returns the y-component.
    #[inline]
    pub fn y(&self) -> Quantity<U> {
        self.storage.xyz.y()
    }

    /// Returns the z-component.
    #[inline]
    pub fn z(&self) -> Quantity<U> {
        self.storage.xyz.z()
    }

    /// Returns the underlying nalgebra Vector3.
    #[inline]
    pub fn as_vec3(&self) -> &nalgebra::Vector3<Quantity<U>> {
        self.storage.xyz.as_vec3()
    }

    /// Returns the internal XYZ storage.
    #[inline]
    pub(crate) fn xyz(&self) -> &XYZ<Quantity<U>> {
        &self.storage.xyz
    }
}

// =============================================================================
// Vector Space Operations
// =============================================================================

impl<F: ReferenceFrame, U: LengthUnit> Displacement<F, U> {
    /// Computes the Euclidean magnitude (length) of the vector.
    #[inline]
    pub fn magnitude(&self) -> Quantity<U> {
        self.storage.xyz.magnitude()
    }

    /// Computes the squared magnitude (avoids sqrt, useful for comparisons).
    #[inline]
    pub fn magnitude_squared(&self) -> f64 {
        self.storage.xyz.magnitude_squared()
    }

    /// Normalizes this vector to a unit direction.
    ///
    /// Returns `None` if the vector has zero (or near-zero) magnitude.
    ///
    /// # Example
    /// ```rust
    /// use siderust::coordinates::algebra::cartesian::Displacement;
    /// use siderust::coordinates::algebra::frames::Ecliptic;
    /// use qtty::*;
    ///
    /// let v = Displacement::<Ecliptic, AstronomicalUnit>::new(3.0, 4.0, 0.0);
    /// let dir = v.normalize().expect("non-zero vector");
    /// // dir is now a unit Direction<Ecliptic>
    /// ```
    #[inline]
    pub fn normalize(&self) -> Option<super::Direction<F>> {
        self.storage.xyz.to_raw().try_normalize().map(super::Direction::from_xyz_unchecked)
    }

    /// Returns a unit direction, assuming non-zero magnitude.
    ///
    /// # Panics
    /// May produce NaN if the vector has zero magnitude.
    #[inline]
    pub fn normalize_unchecked(&self) -> super::Direction<F> {
        super::Direction::from_xyz_unchecked(self.storage.xyz.to_raw().normalize_unchecked())
    }

    /// Scales the vector by a scalar factor.
    #[inline]
    pub fn scale(&self, factor: f64) -> Self {
        Self::from_xyz(XYZ::from_raw(self.storage.xyz.to_raw().scale(factor)))
    }

    /// Computes the dot product with another vector (returns dimensionless f64).
    ///
    /// Note: For proper dimensional analysis, the result would be U²,
    /// but we return the raw f64 value for convenience.
    #[inline]
    pub fn dot(&self, other: &Self) -> f64 {
        self.storage.xyz.to_raw().dot(&other.storage.xyz.to_raw())
    }

    /// Computes the cross product with another vector.
    #[inline]
    pub fn cross(&self, other: &Self) -> Self {
        Self::from_xyz(XYZ::from_raw(
            self.storage.xyz.to_raw().cross(&other.storage.xyz.to_raw()),
        ))
    }

    /// Returns the negation of this vector.
    #[inline]
    pub fn negate(&self) -> Self {
        Self::from_xyz(-self.storage.xyz.clone())
    }
}

// =============================================================================
// Operator Implementations
// =============================================================================

impl<F: ReferenceFrame, U: LengthUnit> Add for Displacement<F, U> {
    type Output = Self;

    /// Displacement + Displacement -> Vector
    #[inline]
    fn add(self, other: Self) -> Self::Output {
        Self::from_xyz(self.storage.xyz + other.storage.xyz)
    }
}

impl<F: ReferenceFrame, U: LengthUnit> Sub for Displacement<F, U> {
    type Output = Self;

    /// Displacement - Displacement -> Vector
    #[inline]
    fn sub(self, other: Self) -> Self::Output {
        Self::from_xyz(self.storage.xyz - other.storage.xyz)
    }
}

impl<F: ReferenceFrame, U: LengthUnit> Neg for Displacement<F, U> {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self::Output {
        self.negate()
    }
}

// =============================================================================
// Display
// =============================================================================

impl<F: ReferenceFrame, U: LengthUnit> std::fmt::Display for Displacement<F, U>
where
    Quantity<U>: std::fmt::Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Position<{}> X: {:.6}, Y: {:.6}, Z: {:.6}",
            F::frame_name(),
            self.x(),
            self.y(),
            self.z()
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coordinates::algebra::frames::Ecliptic;
    use qtty::AstronomicalUnit;

    type DispAu = Displacement<Ecliptic, AstronomicalUnit>;

    #[test]
    fn test_displacement_add_sub() {
        let a = DispAu::new(1.0, 2.0, 3.0);
        let b = DispAu::new(4.0, 5.0, 6.0);

        let sum = a + b;
        assert!((sum.x().value() - 5.0).abs() < 1e-12);
        assert!((sum.y().value() - 7.0).abs() < 1e-12);
        assert!((sum.z().value() - 9.0).abs() < 1e-12);

        let diff = b - a;
        assert!((diff.x().value() - 3.0).abs() < 1e-12);
        assert!((diff.y().value() - 3.0).abs() < 1e-12);
        assert!((diff.z().value() - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_displacement_magnitude() {
        let v = DispAu::new(3.0, 4.0, 0.0);
        assert!((v.magnitude().value() - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_displacement_normalize() {
        let v = DispAu::new(3.0, 4.0, 0.0);
        let dir = v.normalize().expect("non-zero displacement");
        let norm = (dir.x() * dir.x() + dir.y() * dir.y() + dir.z() * dir.z()).sqrt();
        assert!((norm - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_zero_displacement_normalize() {
        let zero = DispAu::ZERO;
        assert!(zero.normalize().is_none());
    }
}
