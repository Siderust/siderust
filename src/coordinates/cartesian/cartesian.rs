//! # Cartesian Coordinates
//!
//! This module defines the generic [`Vector<C, F, K>`] type for representing 3D positions or directions
//! in astronomical reference frames and centers, with strong compile-time type safety.
//!
//! ## Overview
//!
//! - **Generic over Center, Frame, and Kind:**
//!   - `C`: Reference center (e.g., `Heliocentric`, `Geocentric`).
//!   - `F`: Reference frame (e.g., `ICRS`, `Ecliptic`, `Equatorial`).
//!   - `K`: Kind marker (`Position`, `Direction`), enforcing semantic correctness.
//! - **Type Safety:** Operations are only allowed between coordinates with matching type parameters.
//! - **Units:** Coordinates are expressed in astronomical units (AU) by convention, but may represent other units if documented.
//! - **Vector Operations:** Supports addition, subtraction, scaling, and distance calculation.
//! - **Interoperability:** Seamless conversion to and from `nalgebra::Vector3<U>`.
//!
//! ## Example
//! ```rust
//! use siderust::coordinates::cartesian::Vector;
//! use siderust::coordinates::centers::Heliocentric;
//! use siderust::coordinates::frames::Ecliptic;
//! use siderust::coordinates::kinds::PositionKind;
//!
//! // Create a heliocentric ecliptic position
//! let pos = Vector::<Heliocentric, Ecliptic, f64, PositionKind>::new(1.0, 0.0, 0.0);
//! println!("X: {}, Y: {}, Z: {}", pos.x(), pos.y(), pos.z());
//! ```
//!
//! ## Type Aliases
//! You may define type aliases for common systems, e.g.,
//! `type HeliocentricEcliptic = Vector<Heliocentric, Ecliptic, PositionKind>;`

use crate::coordinates::{
    frames, centers,
    kinds::*,
};

use std::marker::PhantomData;
use nalgebra::Vector3;
use std::ops::{Add, Sub};
use crate::units::{Quantity, LengthUnit};

// Refactoriza la definición de Vector:
#[derive(Debug, Clone, Copy)]
pub struct Vector<
    C: centers::ReferenceCenter,
    F: frames::ReferenceFrame,
    U: LengthUnit,
    K: Kind = PositionKind,
> {
    xyz: Vector3<Quantity<U>>,
    _center: PhantomData<C>,
    _frame: PhantomData<F>,
    _kind: PhantomData<K>,
}

// Refactoriza los métodos y traits:
impl<C, F, U, K> Vector<C, F, U, K>
where
    C: centers::ReferenceCenter,
    F: frames::ReferenceFrame,
    U: LengthUnit,
    K: Kind,
{
    /// Creates a new Cartesian coordinate.
    ///
    /// # Arguments
    /// - `x`: The x-coordinate in AU.
    /// - `y`: The y-coordinate in AU.
    /// - `z`: The z-coordinate in AU.
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
        Vector { xyz: vec3, _center: PhantomData, _frame: PhantomData, _kind: PhantomData }
    }

    pub const fn as_vec3(&self) -> Vector3<Quantity<U>> { self.xyz }

    /// Gets the x-coordinate in AU.
    pub fn x(&self) -> Quantity<U> { self.xyz[0] }

    /// Gets the y-coordinate in AU.
    pub fn y(&self) -> Quantity<U> { self.xyz[1] }

    /// Gets the z-coordinate in AU.
    pub fn z(&self) -> Quantity<U> { self.xyz[2] }

    pub fn sub(&self, other: &Self) -> Self
    where
        U: std::cmp::PartialEq + std::fmt::Debug
    {
        Self::from_vec3(self.as_vec3() - other.as_vec3())
    }

    /// Calculates the Euclidean distance with respect to the ReferenceCenter.
    ///
    /// # Returns
    /// The distance from the ReferenceCenter in units of U.
    pub fn distance(&self) -> Quantity<U> {
        let distance = Vector3::<f64>::new(
            self.x().value(),
            self.y().value(),
            self.z().value()
        ).magnitude();
        Quantity::new(distance)
    }

    /// Computes the Euclidean distance to another Cartesian coordinate of the same type.
    pub fn distance_to(&self, other: &Self) -> Quantity<U>
    where
        U: std::cmp::PartialEq + std::fmt::Debug
    {
        self.sub(other).distance()
    }

}


impl<C, F, U, K> std::fmt::Display for Vector<C, F, U, K>
where
    C: centers::ReferenceCenter,
    F: frames::ReferenceFrame,
    U: LengthUnit,
    K: Kind,
    Quantity<U>: std::fmt::Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Center: {}, Frame: {}, X: {:.6}, Y: {:.6}, Z: {:.6}",
            C::center_name(),
            F::frame_name(),
            self.x(), self.y(), self.z()
        )
    }
}

impl<C, F, U, K> Add for Vector<C, F, U, K>
where
    C: centers::ReferenceCenter,
    F: frames::ReferenceFrame,
    U: LengthUnit + std::cmp::PartialEq + std::fmt::Debug,
    K: Kind,
{
    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        Self::from_vec3(self.xyz + other.xyz)
    }
}

impl<C, F, U, K> Sub for Vector<C, F, U, K>
where
    C: centers::ReferenceCenter,
    F: frames::ReferenceFrame,
    U: LengthUnit + std::cmp::PartialEq + std::fmt::Debug,
    K: Kind,
{
    type Output = Self;
    fn sub(self, other: Self) -> Self::Output { (&self).sub(&other) }
}

/*impl<C, F, U, K> Div<U> for Vector<C, F, U, K>
where
    C: centers::ReferenceCenter,
    F: frames::ReferenceFrame,
    U: LengthUnit,
    K: Kind,
{
    type Output = Self;

    fn div(self, rhs: U) -> Self::Output {
        Self::from_vec3(self.xyz / rhs.as_vec3())
    }
}

impl<C, F, U, K> Mul<U> for Vector<C, F, U, K>
where
    C: centers::ReferenceCenter,
    F: frames::ReferenceFrame,
    U: LengthUnit,
    K: Kind,
{
    type Output = Self;

    fn mul(self, rhs: U) -> Self::Output {
        Self::from_vec3(self.xyz * rhs.into())
    }
}
*/
