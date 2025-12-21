//! # Cartesian Free Vector
//!
//! This module defines [`Vector<F, U>`], a **free vector** in 3D Cartesian space.
//!
//! ## Mathematical Model
//!
//! A free vector represents a directed magnitude in space. It is:
//! - **Frame-dependent**: The orientation is relative to a reference frame `F`
//! - **Center-independent**: Free vectors are translation-invariant
//! - **Dimensioned**: Has a unit `U` (length, velocity, acceleration, etc.)
//!
//! Free vectors form a vector space and support:
//! - Addition: `Vector + Vector -> Vector`
//! - Subtraction: `Vector - Vector -> Vector`
//! - Scalar multiplication: `Vector * f64 -> Vector`
//! - Normalization (for length units): `normalize(Vector) -> Option<Direction>`
//!
//! ## Semantic Type Aliases
//!
//! This module provides semantic aliases for clarity:
//! - [`Displacement<F, U>`] — displacement vector with length unit
//! - [`Velocity<F, U>`] — velocity vector with velocity unit
//!
//! ## Example
//!
//! ```rust
//! use siderust::coordinates::algebra::cartesian::{Vector, Displacement, Velocity};
//! use siderust::coordinates::algebra::frames::Ecliptic;
//! use qtty::*;
//!
//! // Displacement with length unit
//! let displacement = Displacement::<Ecliptic, AstronomicalUnit>::new(1.0, 2.0, 3.0);
//!
//! // Velocity with velocity unit
//! type AuPerDay = Per<AstronomicalUnit, Day>;
//! let velocity = Velocity::<Ecliptic, AuPerDay>::new(0.01, 0.02, 0.0);
//!
//! // Both are just Vector<F, U> with different units
//! let v: Vector<Ecliptic, AuPerDay> = velocity;
//! ```

use super::xyz::XYZ;
use crate::coordinates::algebra::frames::ReferenceFrame;
use qtty::{LengthUnit, Quantity, Unit};

use std::marker::PhantomData;
use std::ops::{Add, Neg, Sub};

/// A free vector in 3D Cartesian coordinates.
///
/// Free vectors are frame-dependent but center-independent. They represent
/// directed magnitudes (displacement, velocity, acceleration, etc.) in space.
///
/// # Type Parameters
/// - `F`: The reference frame (e.g., `ICRS`, `Ecliptic`, `Equatorial`)
/// - `U`: The unit (e.g., `AstronomicalUnit`, `Per<Kilometer, Second>`)
///
/// # Zero-Cost Abstraction
///
/// This type uses `#[repr(transparent)]` over the internal storage,
/// ensuring no runtime overhead compared to raw `Vector3<Quantity<U>>`.
#[derive(Debug, Clone, Copy)]
#[repr(transparent)]
pub struct Vector<F: ReferenceFrame, U: Unit> {
    storage: VectorStorage<F, U>,
}

/// Internal storage for Vector with PhantomData marker.
#[derive(Debug, Clone, Copy)]
struct VectorStorage<F: ReferenceFrame, U: Unit> {
    xyz: XYZ<Quantity<U>>,
    _frame: PhantomData<F>,
}

// =============================================================================
// Constructors
// =============================================================================

impl<F: ReferenceFrame, U: Unit> Vector<F, U> {
    /// Creates a new vector from component values.
    ///
    /// # Arguments
    /// - `x`, `y`, `z`: Component values (converted to `Quantity<U>`)
    #[inline]
    pub fn new<T: Into<Quantity<U>>>(x: T, y: T, z: T) -> Self {
        Self {
            storage: VectorStorage {
                xyz: XYZ::new(x.into(), y.into(), z.into()),
                _frame: PhantomData,
            },
        }
    }

    /// Creates a vector from a nalgebra Vector3.
    #[inline]
    pub fn from_vec3(vec3: nalgebra::Vector3<Quantity<U>>) -> Self {
        Self {
            storage: VectorStorage {
                xyz: XYZ::from_vec3(vec3),
                _frame: PhantomData,
            },
        }
    }

    /// Creates a vector from the internal XYZ storage.
    #[inline]
    pub(crate) fn from_xyz(xyz: XYZ<Quantity<U>>) -> Self {
        Self {
            storage: VectorStorage {
                xyz,
                _frame: PhantomData,
            },
        }
    }

    /// The zero vector.
    pub const ZERO: Self = Self {
        storage: VectorStorage {
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

impl<F: ReferenceFrame, U: Unit> Vector<F, U> {
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
}

// =============================================================================
// Vector Space Operations (for all units)
// =============================================================================

impl<F: ReferenceFrame, U: Unit> Vector<F, U> {
    /// Computes the Euclidean magnitude of the vector.
    #[inline]
    pub fn magnitude(&self) -> Quantity<U> {
        self.storage.xyz.magnitude()
    }

    /// Computes the squared magnitude (avoids sqrt, useful for comparisons).
    #[inline]
    pub fn magnitude_squared(&self) -> f64 {
        self.storage.xyz.magnitude_squared()
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
// Length-Specific Operations (only for LengthUnit)
// =============================================================================

impl<F: ReferenceFrame, U: LengthUnit> Vector<F, U> {
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
        self.storage
            .xyz
            .to_raw()
            .try_normalize()
            .map(super::Direction::from_xyz_unchecked)
    }

    /// Returns a unit direction, assuming non-zero magnitude.
    ///
    /// # Panics
    /// May produce NaN if the vector has zero magnitude.
    #[inline]
    pub fn normalize_unchecked(&self) -> super::Direction<F> {
        super::Direction::from_xyz_unchecked(self.storage.xyz.to_raw().normalize_unchecked())
    }
}

// =============================================================================
// Operator Implementations
// =============================================================================

impl<F: ReferenceFrame, U: Unit> Add for Vector<F, U> {
    type Output = Self;

    /// Vector + Vector -> Vector
    #[inline]
    fn add(self, other: Self) -> Self::Output {
        Self::from_xyz(self.storage.xyz + other.storage.xyz)
    }
}

impl<F: ReferenceFrame, U: Unit> Sub for Vector<F, U> {
    type Output = Self;

    /// Vector - Vector -> Vector
    #[inline]
    fn sub(self, other: Self) -> Self::Output {
        Self::from_xyz(self.storage.xyz - other.storage.xyz)
    }
}

impl<F: ReferenceFrame, U: Unit> Neg for Vector<F, U> {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self::Output {
        self.negate()
    }
}

// =============================================================================
// PartialEq
// =============================================================================

impl<F: ReferenceFrame, U: Unit> PartialEq for Vector<F, U> {
    fn eq(&self, other: &Self) -> bool {
        self.storage.xyz == other.storage.xyz
    }
}

// =============================================================================
// Display
// =============================================================================

impl<F: ReferenceFrame, U: Unit> std::fmt::Display for Vector<F, U>
where
    Quantity<U>: std::fmt::Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Vector<{}> X: {:.6}, Y: {:.6}, Z: {:.6}",
            F::frame_name(),
            self.x(),
            self.y(),
            self.z()
        )
    }
}

// =============================================================================
// Type Aliases
// =============================================================================

/// A displacement vector (free vector with length unit).
///
/// This is a semantic alias for [`Vector<F, U>`] where `U` is a length unit.
/// Displacements represent directed distances in space.
pub type Displacement<F, U> = Vector<F, U>;

/// A velocity vector (free vector with velocity unit).
///
/// This is a semantic alias for [`Vector<F, U>`] where `U` is a velocity unit.
/// Velocities represent rates of change of position.
pub type Velocity<F, U> = Vector<F, U>;

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coordinates::algebra::frames::Ecliptic;
    use qtty::{AstronomicalUnit, Day, Per};

    type DispAu = Displacement<Ecliptic, AstronomicalUnit>;
    type AuPerDay = Per<AstronomicalUnit, Day>;
    type VelAuDay = Velocity<Ecliptic, AuPerDay>;

    #[test]
    fn test_vector_add_sub() {
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
    fn test_vector_magnitude() {
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

    #[test]
    fn test_velocity_add_sub() {
        let v1 = VelAuDay::new(
            Quantity::<AuPerDay>::new(1.0),
            Quantity::<AuPerDay>::new(2.0),
            Quantity::<AuPerDay>::new(3.0),
        );
        let v2 = VelAuDay::new(
            Quantity::<AuPerDay>::new(0.5),
            Quantity::<AuPerDay>::new(1.0),
            Quantity::<AuPerDay>::new(1.5),
        );

        let sum = v1 + v2;
        assert!((sum.x().value() - 1.5).abs() < 1e-12);
        assert!((sum.y().value() - 3.0).abs() < 1e-12);
        assert!((sum.z().value() - 4.5).abs() < 1e-12);

        let diff = v1 - v2;
        assert!((diff.x().value() - 0.5).abs() < 1e-12);
        assert!((diff.y().value() - 1.0).abs() < 1e-12);
        assert!((diff.z().value() - 1.5).abs() < 1e-12);
    }

    #[test]
    fn test_velocity_magnitude() {
        let v = VelAuDay::new(
            Quantity::<AuPerDay>::new(3.0),
            Quantity::<AuPerDay>::new(4.0),
            Quantity::<AuPerDay>::new(0.0),
        );
        assert!((v.magnitude().value() - 5.0).abs() < 1e-12);
    }
}
