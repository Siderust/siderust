//! # Cartesian Velocity (Free Vectors)
//!
//! This module defines [`Velocity<F, U>`], a **free velocity vector** representing
//! rate of change of position in 3D space.
//!
//! ## Mathematical Model
//!
//! A velocity is a free vector representing the derivative of position with respect
//! to time. Like displacement vectors, velocities are:
//!
//! - **Frame-dependent**: Direction is relative to frame `F`
//! - **Center-independent**: Velocities are translation-invariant
//! - **Dimensioned**: Has velocity unit `U` (e.g., `Per<Kilometer, Second>`)
//!
//! ## Supported Operations
//!
//! | Operation | Result | Meaning |
//! |-----------|--------|---------|
//! | `Velocity + Velocity` | `Velocity` | Add velocities |
//! | `Velocity - Velocity` | `Velocity` | Velocity difference |
//! | `Velocity * scalar` | `Velocity` | Scale velocity |
//!
//! ## Example
//!
//! ```rust
//! use siderust::coordinates::algebra::cartesian::Velocity;
//! use siderust::coordinates::algebra::frames::Ecliptic;
//! use qtty::*;
//!
//! type AuPerDay = Per<AstronomicalUnit, Day>;
//!
//! let v1 = Velocity::<Ecliptic, AuPerDay>::new(0.01, 0.02, 0.0);
//! let v2 = Velocity::<Ecliptic, AuPerDay>::new(0.005, 0.01, 0.0);
//!
//! let sum = v1 + v2;
//! ```

use super::xyz::XYZ;
use crate::coordinates::algebra::frames::ReferenceFrame;
use qtty::{Quantity, Unit};

use std::marker::PhantomData;
use std::ops::{Add, Neg, Sub};

/// A velocity vector in 3D Cartesian coordinates.
///
/// Velocities are frame-dependent but center-independent (free vectors).
/// They represent rates of change and cannot undergo center transformations.
///
/// # Type Parameters
/// - `F`: The reference frame (e.g., `ICRS`, `Ecliptic`, `Equatorial`)
/// - `U`: The velocity unit (e.g., `Per<AstronomicalUnit, Day>`)
#[derive(Debug, Clone, Copy)]
#[repr(transparent)]
pub struct Velocity<F: ReferenceFrame, U: Unit> {
    storage: VelocityStorage<F, U>,
}

/// Internal storage for Velocity with PhantomData marker.
#[derive(Debug, Clone, Copy)]
struct VelocityStorage<F: ReferenceFrame, U: Unit> {
    xyz: XYZ<Quantity<U>>,
    _frame: PhantomData<F>,
}

// =============================================================================
// Constructors
// =============================================================================

impl<F: ReferenceFrame, U: Unit> Velocity<F, U> {
    /// Creates a new velocity vector from component values.
    #[inline]
    pub fn new<T: Into<Quantity<U>>>(x: T, y: T, z: T) -> Self {
        Self {
            storage: VelocityStorage {
                xyz: XYZ::new(x.into(), y.into(), z.into()),
                _frame: PhantomData,
            },
        }
    }

    /// Creates a velocity from a nalgebra Vector3.
    #[inline]
    pub fn from_vec3(vec3: nalgebra::Vector3<Quantity<U>>) -> Self {
        Self {
            storage: VelocityStorage {
                xyz: XYZ::from_vec3(vec3),
                _frame: PhantomData,
            },
        }
    }

    /// Creates a velocity from internal XYZ storage.
    #[inline]
    fn from_xyz(xyz: XYZ<Quantity<U>>) -> Self {
        Self {
            storage: VelocityStorage {
                xyz,
                _frame: PhantomData,
            },
        }
    }

    /// The zero velocity.
    pub const ZERO: Self = Self {
        storage: VelocityStorage {
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

impl<F: ReferenceFrame, U: Unit> Velocity<F, U> {
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
// Operations
// =============================================================================

impl<F: ReferenceFrame, U: Unit> Velocity<F, U> {
    /// Computes the magnitude (speed) of this velocity vector.
    #[inline]
    pub fn magnitude(&self) -> Quantity<U> {
        self.storage.xyz.magnitude()
    }

    /// Computes the squared magnitude.
    #[inline]
    pub fn magnitude_squared(&self) -> f64 {
        self.storage.xyz.magnitude_squared()
    }

    /// Scales the velocity by a scalar factor.
    #[inline]
    pub fn scale(&self, factor: f64) -> Self {
        Self::from_xyz(XYZ::from_raw(self.storage.xyz.to_raw().scale(factor)))
    }

    /// Returns the negation of this velocity.
    #[inline]
    pub fn negate(&self) -> Self {
        Self::from_xyz(-self.storage.xyz.clone())
    }
}

// =============================================================================
// Operator Implementations
// =============================================================================

impl<F: ReferenceFrame, U: Unit> Add for Velocity<F, U> {
    type Output = Self;

    #[inline]
    fn add(self, other: Self) -> Self::Output {
        Self::from_xyz(self.storage.xyz + other.storage.xyz)
    }
}

impl<F: ReferenceFrame, U: Unit> Sub for Velocity<F, U> {
    type Output = Self;

    #[inline]
    fn sub(self, other: Self) -> Self::Output {
        Self::from_xyz(self.storage.xyz - other.storage.xyz)
    }
}

impl<F: ReferenceFrame, U: Unit> Neg for Velocity<F, U> {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self::Output {
        self.negate()
    }
}

// =============================================================================
// Display
// =============================================================================

impl<F: ReferenceFrame, U: Unit> std::fmt::Display for Velocity<F, U>
where
    Quantity<U>: std::fmt::Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Velocity<{}> Vx: {:.6}, Vy: {:.6}, Vz: {:.6}",
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
    use qtty::{AstronomicalUnit, Day, Per};

    type AuPerDay = Per<AstronomicalUnit, Day>;
    type VelAuDay = Velocity<Ecliptic, AuPerDay>;

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
