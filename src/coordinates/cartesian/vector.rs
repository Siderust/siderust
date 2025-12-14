//! # Cartesian Coordinates
//!
//! This module defines the generic [`Vector<C, F>`] type for representing 3D positions or directions
//! in astronomical reference frames and centers, with strong compile-time type safety.
//!
//! ## Overview
//!
//! - **Generic over Center, Frame, and Kind:**
//!   - `C`: Reference center (e.g., `Heliocentric`, `Geocentric`).
//!   - `F`: Reference frame (e.g., `ICRS`, `Ecliptic`, `Equatorial`).
//! - **Type Safety:** Operations are only allowed between coordinates with matching type parameters.
//! - **Units:** Coordinates are expressed in units (Length and Velocity).
//! - **Vector Operations:** Supports addition, subtraction, scaling, and distance calculation.
//!
//! ## Example
//! ```rust
//! use siderust::coordinates::cartesian::position::Ecliptic;
//! use qtty::AstronomicalUnit;
//!
//! // Create a heliocentric ecliptic position
//! let pos = Ecliptic::<AstronomicalUnit>::new(1.0, 0.0, 0.0);
//! println!("X: {}, Y: {}, Z: {}", pos.x(), pos.y(), pos.z());
//! ```
//!
//! ## Type Aliases
//! We define type aliases for common systems, e.g.,
//! `type Distance<C: ReferenceCenter, F: ReferenceFrame> = Vector<C, F, Unitless>;`
//! `type Velocity<C: ReferenceCenter, F: ReferenceFrame, U: VelocityUnit> = Vector<C, F, U>;`
//! `type Position<C: ReferenceCenter, F: ReferenceFrame, U: LengthUnit>   = Vector<C, F, U>;`

use crate::coordinates::{centers, frames};
use qtty::*;

use nalgebra::Vector3;
use std::marker::PhantomData;
use std::ops::{Add, Sub};

#[derive(Debug, Clone, Copy)]
pub struct Vector<C: centers::ReferenceCenter, F: frames::ReferenceFrame, U: Unit> {
    xyz: nalgebra::Vector3<Quantity<U>>,
    _center: PhantomData<C>,
    _frame: PhantomData<F>,
}

impl<C, F, U> Vector<C, F, U>
where
    C: centers::ReferenceCenter,
    F: frames::ReferenceFrame,
    U: Unit,
{
    /// Creates a new Cartesian coordinate.
    ///
    /// # Arguments
    /// - `x`: The x-coordinate in AstronomicalUnits.
    /// - `y`: The y-coordinate in AstronomicalUnits.
    /// - `z`: The z-coordinate in AstronomicalUnits.
    ///
    /// # Returns
    /// A new `Vector<Center, Frame>`.
    pub const fn new_const(x: Quantity<U>, y: Quantity<U>, z: Quantity<U>) -> Self {
        Self::from_vec3(Vector3::new(x, y, z))
    }

    pub fn new<T>(x: T, y: T, z: T) -> Self
    where
        T: Into<Quantity<U>>,
    {
        Self::from_vec3(Vector3::new(x.into(), y.into(), z.into()))
    }

    pub const fn from_vec3(vec3: Vector3<Quantity<U>>) -> Self {
        Vector {
            xyz: vec3,
            _center: PhantomData,
            _frame: PhantomData,
        }
    }

    pub const fn as_vec3(&self) -> Vector3<Quantity<U>> {
        self.xyz
    }

    /// Gets the x-coordinate in AstronomicalUnits.
    pub fn x(&self) -> Quantity<U> {
        self.xyz[0]
    }

    /// Gets the y-coordinate in AstronomicalUnits.
    pub fn y(&self) -> Quantity<U> {
        self.xyz[1]
    }

    /// Gets the z-coordinate in AstronomicalUnits.
    pub fn z(&self) -> Quantity<U> {
        self.xyz[2]
    }

    pub fn sub(&self, other: &Self) -> Self
    where
        U: std::cmp::PartialEq + std::fmt::Debug,
    {
        Self::from_vec3(self.as_vec3() - other.as_vec3())
    }

    /// Calculates the Euclidean distance with respect to the ReferenceCenter.
    ///
    /// # Returns
    /// The distance from the ReferenceCenter in units of U.
    pub fn distance(&self) -> Quantity<U> {
        let distance =
            Vector3::<f64>::new(self.x().value(), self.y().value(), self.z().value()).magnitude();
        Quantity::new(distance)
    }

    /// Computes the Euclidean distance to another vector of the same type.
    pub fn distance_to(&self, other: &Self) -> Quantity<U>
    where
        U: std::cmp::PartialEq + std::fmt::Debug,
    {
        self.sub(other).distance()
    }
}

impl<C, F, U> std::fmt::Display for Vector<C, F, U>
where
    C: centers::ReferenceCenter,
    F: frames::ReferenceFrame,
    U: Unit,
    Quantity<U>: std::fmt::Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Center: {}, Frame: {}, X: {:.6}, Y: {:.6}, Z: {:.6}",
            C::center_name(),
            F::frame_name(),
            self.x(),
            self.y(),
            self.z()
        )
    }
}

impl<C, F, U> Add for Vector<C, F, U>
where
    C: centers::ReferenceCenter,
    F: frames::ReferenceFrame,
    U: Unit,
{
    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        Self::from_vec3(self.xyz + other.xyz)
    }
}

impl<C, F, U> Sub for Vector<C, F, U>
where
    C: centers::ReferenceCenter,
    F: frames::ReferenceFrame,
    U: Unit,
{
    type Output = Self;
    fn sub(self, other: Self) -> Self::Output {
        Self::sub(&self, &other)
    }
}
