//! # Cartesian Velocity (Free Vectors)
//!
//! This module defines velocity types - free vectors representing rates of change
//! in space without any spatial origin (center).
//!
//! ## Mathematical Foundations
//!
//! Velocities are **free vectors**: they represent rates of change (derivatives
//! of position with respect to time) and are translation-invariant:
//!
//! - Rotating a velocity is valid (frame transformation)
//! - "Translating" a velocity to a different origin is **meaningless**
//!
//! Therefore, `Velocity<F, U>` has no center parameter, only a frame `F` and
//! a velocity unit `U`.
//!
//! ## Example
//!
//! ```rust
//! use siderust::coordinates::cartesian::Velocity;
//! use siderust::coordinates::frames::Ecliptic;
//! use qtty::*;
//!
//! type AuPerDay = qtty::Per<AstronomicalUnit, Day>;
//!
//! let vel = Velocity::<Ecliptic, AuPerDay>::new(
//!     qtty::velocity::Velocity::<AstronomicalUnit, Day>::new(0.01),
//!     qtty::velocity::Velocity::<AstronomicalUnit, Day>::new(0.02),
//!     qtty::velocity::Velocity::<AstronomicalUnit, Day>::new(0.0),
//! );
//! ```

use crate::coordinates::frames::{self, MutableFrame};
use qtty::{Quantity, Unit};

use nalgebra::Vector3;
use std::marker::PhantomData;
use std::ops::{Add, Sub};

/// A velocity vector in 3D space.
///
/// Velocities are frame-dependent but center-independent (free vectors).
/// They cannot undergo center transformations, only frame transformations.
///
/// # Type Parameters
/// - `F`: The reference frame (e.g., `ICRS`, `Ecliptic`, `Equatorial`)
/// - `U`: The velocity unit (e.g., `Per<AstronomicalUnit, Day>`)
#[derive(Debug, Clone, Copy)]
pub struct Velocity<F: frames::ReferenceFrame, U: Unit> {
    xyz: Vector3<Quantity<U>>,
    _frame: PhantomData<F>,
}

impl<F, U> Velocity<F, U>
where
    F: frames::ReferenceFrame,
    U: Unit,
{
    /// Creates a new velocity vector.
    pub fn new<T: Into<Quantity<U>>>(x: T, y: T, z: T) -> Self {
        Self {
            xyz: Vector3::new(x.into(), y.into(), z.into()),
            _frame: PhantomData,
        }
    }

    /// Creates a velocity from a nalgebra Vector3.
    pub fn from_vec3(vec3: Vector3<Quantity<U>>) -> Self {
        Self {
            xyz: vec3,
            _frame: PhantomData,
        }
    }

    /// Returns the underlying nalgebra Vector3.
    pub const fn as_vec3(&self) -> Vector3<Quantity<U>> {
        self.xyz
    }

    /// Gets the x-component.
    pub fn x(&self) -> Quantity<U> {
        self.xyz[0]
    }

    /// Gets the y-component.
    pub fn y(&self) -> Quantity<U> {
        self.xyz[1]
    }

    /// Gets the z-component.
    pub fn z(&self) -> Quantity<U> {
        self.xyz[2]
    }

    /// Calculates the magnitude (speed) of this velocity vector.
    pub fn magnitude(&self) -> Quantity<U> {
        let mag =
            Vector3::<f64>::new(self.x().value(), self.y().value(), self.z().value()).magnitude();
        Quantity::new(mag)
    }
}

impl<F, U> Add for Velocity<F, U>
where
    F: frames::ReferenceFrame,
    U: Unit,
{
    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        Self::from_vec3(self.xyz + other.xyz)
    }
}

impl<F, U> Sub for Velocity<F, U>
where
    F: frames::ReferenceFrame,
    U: Unit,
{
    type Output = Self;

    fn sub(self, other: Self) -> Self::Output {
        Self::from_vec3(self.xyz - other.xyz)
    }
}

impl<F, U> std::fmt::Display for Velocity<F, U>
where
    F: frames::ReferenceFrame,
    U: Unit,
    Quantity<U>: std::fmt::Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Frame: {}, Vx: {:.6}, Vy: {:.6}, Vz: {:.6}",
            F::frame_name(),
            self.x(),
            self.y(),
            self.z()
        )
    }
}

/// Type aliases for common velocity systems.
pub type Ecliptic<U> = Velocity<frames::Ecliptic, U>;
pub type Equatorial<U> = Velocity<frames::Equatorial, U>;
pub type Horizontal<U> = Velocity<frames::Horizontal, U>;
pub type ICRS<U> = Velocity<frames::ICRS, U>;
