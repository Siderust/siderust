//! Defines the generic [`CartesianCoord<Center, Frame>`] type for representing 3D+time positions
//! in various astronomical reference frames and centers.
//!
//! This module provides:
//! - A generic Cartesian coordinate struct with time (Julian Date).
//! - Type aliases for commonly used systems like ICRS, HCRS, Ecliptic, etc.
//! - Distance calculation methods and string formatting.
//!
//! Coordinates are expressed in abstract units (KM, LY, AU, ...) and Julian Days (JD).
//! The type system enforces correct pairing of reference centers and frames.

use crate::coordinates::{
    frames, centers,
    kinds::*,
};

use std::marker::PhantomData;
use nalgebra::Vector3;
use std::ops::{Add, Sub, Div, Mul};

/// A Cartesian coordinate representation with a specific reference center and frame.
///
/// # Type Parameters
/// - `Center`: The reference center (e.g., Barycentric, Heliocentric).
/// - `Frame`: The reference frame (e.g., ICRS, Ecliptic).
#[derive(Debug, Clone, Copy)]
pub struct CartesianCoord<
    C: centers::ReferenceCenter,
    F: frames::ReferenceFrame,
    K : Kind
> {

    xyz: Vector3<f64>,
    _center: PhantomData<C>,
    _frame : PhantomData<F>,
    _kind  : PhantomData<K>,
}

impl<C, F, K> CartesianCoord<C, F, K>
where
    C: centers::ReferenceCenter,
    F: frames::ReferenceFrame,
    K: Kind,
{
    /// Creates a new Cartesian coordinate.
    ///
    /// # Arguments
    /// - `x`: The x-coordinate in AU.
    /// - `y`: The y-coordinate in AU.
    /// - `z`: The z-coordinate in AU.
    /// - `t`: The time coordinate (e.g., Julian Date).
    ///
    /// # Returns
    /// A new `CartesianCoord<Center, Frame>`.
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        // El inliner eliminará por completo la llamada para PositionKind
        K::validate(x, y, z);
        Self::from_vec3(Vector3::new(x, y, z))
    }

    pub const fn from_vec3(vec3: Vector3<f64>) -> Self {
        CartesianCoord { xyz: vec3, _center: PhantomData, _frame: PhantomData, _kind : PhantomData }
    }

    pub const fn as_vec3(&self) -> Vector3<f64> { self.xyz }

    /// Gets the x-coordinate in AU.
    pub fn x(&self) -> f64 { self.xyz[0] }

    /// Gets the y-coordinate in AU.
    pub fn y(&self) -> f64 { self.xyz[1] }

    /// Gets the z-coordinate in AU.
    pub fn z(&self) -> f64 { self.xyz[2] }

    /// Computes the Euclidean distance to another Cartesian coordinate of the same type.
    pub fn distance_to(
        &self,
        other: &Self,
    ) -> f64 {
        ((self.x() - other.x()).powi(2) + (self.y() - other.y()).powi(2) + (self.z() - other.z()).powi(2)).sqrt()
    }
}

impl<C, F, K> std::fmt::Display for CartesianCoord<C, F, K>
where
    C: centers::ReferenceCenter,
    F: frames::ReferenceFrame,
    K: Kind,
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


impl<C, F, K> Add for CartesianCoord<C, F, K>
where
    C: centers::ReferenceCenter,
    F: frames::ReferenceFrame,
    K: Kind,
{
    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        Self::from_vec3(self.xyz + other.xyz)
    }
}

impl<C, F, K> Sub for CartesianCoord<C, F, K>
where
    C: centers::ReferenceCenter,
    F: frames::ReferenceFrame,
    K: Kind,
{
    type Output = Self;

    fn sub(self, other: Self) -> Self::Output {
        Self::from_vec3(self.xyz - other.xyz)
    }
}

impl<C, F, K> Div<f64> for CartesianCoord<C, F, K>
where
    C: centers::ReferenceCenter,
    F: frames::ReferenceFrame,
    K: Kind,
{
    type Output = Self;

    fn div(self, rhs: f64) -> Self::Output {
        Self::from_vec3(self.xyz / rhs)
    }
}

impl<C, F, K> Mul<f64> for CartesianCoord<C, F, K>
where
    C: centers::ReferenceCenter,
    F: frames::ReferenceFrame,
    K: Kind,
{
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        Self::from_vec3(self.xyz * rhs)
    }
}
