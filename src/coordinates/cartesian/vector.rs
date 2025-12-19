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
//! - **Center Parameters:** Each vector stores `C::Params` to support parameterized centers.
//!   For most centers, `Params = ()` (zero-cost). For `Topocentric`, it stores observer location.
//! - **Type Safety:** Operations are only allowed between coordinates with matching type parameters.
//! - **Units:** Coordinates are expressed in units (Length and Velocity).
//! - **Vector Operations:** Supports addition, subtraction, scaling, and distance calculation.
//!
//! ## Example
//! ```rust
//! use siderust::coordinates::cartesian::position::Ecliptic;
//! use qtty::AstronomicalUnit;
//!
//! // Create a heliocentric ecliptic position (Params=() is automatic)
//! let pos = Ecliptic::<AstronomicalUnit>::new(1.0, 0.0, 0.0);
//! println!("X: {}, Y: {}, Z: {}", pos.x(), pos.y(), pos.z());
//! ```
//!
//! ## Type Aliases
//! We define type aliases for common systems, e.g.,
//! `type Direction<F: ReferenceFrame> = Vector<NoCenter, F, DirectionUnit>;`  (free vector)
//! `type Velocity<F: ReferenceFrame, U: VelocityUnit> = Vector<NoCenter, F, U>;`  (free vector)
//! `type Position<C: ReferenceCenter, F: ReferenceFrame, U: LengthUnit> = Vector<C, F, U>;`  (affine)

use crate::coordinates::{centers, frames, math};
use qtty::*;

use nalgebra::Vector3;
use std::marker::PhantomData;
use std::ops::{Add, Sub};

/// A 3D vector in Cartesian coordinates with a specific reference center and frame.
///
/// # Type Parameters
/// - `C`: The reference center (e.g., `Barycentric`, `Heliocentric`, `Geocentric`, `Topocentric`).
/// - `F`: The reference frame (e.g., `ICRS`, `Ecliptic`, `Equatorial`).
/// - `U`: The unit type (e.g., `AstronomicalUnit`, `Kilometer`).
///
/// # Center Parameters
///
/// The vector stores `C::Params` to support parameterized centers:
/// - For `Barycentric`, `Heliocentric`, `Geocentric`: `Params = ()` (zero overhead)
/// - For `Topocentric`: `Params = ObserverSite` (stores observer location)
///
/// Use `new()` / `from_vec3()` with explicit center_params, or use the convenience
/// constructors `new_origin()` / `from_vec3_origin()` for centers with `Params = ()`.
#[derive(Debug, Clone, Copy)]
pub struct Vector<C: centers::ReferenceCenter, F: frames::ReferenceFrame, U: Unit> {
    xyz: nalgebra::Vector3<Quantity<U>>,
    center_params: C::Params,
    _frame: PhantomData<F>,
}

impl<C, F, U> Vector<C, F, U>
where
    C: centers::ReferenceCenter,
    F: frames::ReferenceFrame,
    U: Unit,
{
    /// Creates a new Cartesian coordinate with explicit center parameters.
    ///
    /// # Arguments
    /// - `center_params`: The center parameters (e.g., `()` for most centers, `ObserverSite` for Topocentric).
    /// - `x`: The x-coordinate.
    /// - `y`: The y-coordinate.
    /// - `z`: The z-coordinate.
    ///
    /// # Returns
    /// A new `Vector<Center, Frame, U>`.
    pub const fn new_const(
        center_params: C::Params,
        x: Quantity<U>,
        y: Quantity<U>,
        z: Quantity<U>,
    ) -> Self {
        Self::from_vec3(center_params, Vector3::new(x, y, z))
    }

    /// Creates a new Cartesian coordinate with explicit center parameters.
    ///
    /// # Arguments
    /// - `center_params`: The center parameters.
    /// - `x`, `y`, `z`: Coordinate values (converted to `Quantity<U>`).
    pub fn new_with_params<T>(center_params: C::Params, x: T, y: T, z: T) -> Self
    where
        T: Into<Quantity<U>>,
    {
        Self::from_vec3(center_params, Vector3::new(x.into(), y.into(), z.into()))
    }

    /// Creates a new Cartesian coordinate from a nalgebra Vector3 with explicit center parameters.
    pub const fn from_vec3(center_params: C::Params, vec3: Vector3<Quantity<U>>) -> Self {
        Vector {
            xyz: vec3,
            center_params,
            _frame: PhantomData,
        }
    }

    /// Returns a reference to the center parameters.
    ///
    /// For most centers this returns `&()`. For `Topocentric`, it returns `&ObserverSite`.
    pub fn center_params(&self) -> &C::Params {
        &self.center_params
    }

    pub const fn as_vec3(&self) -> Vector3<Quantity<U>> {
        self.xyz
    }

    /// Gets the x-coordinate.
    pub fn x(&self) -> Quantity<U> {
        self.xyz[0]
    }

    /// Gets the y-coordinate.
    pub fn y(&self) -> Quantity<U> {
        self.xyz[1]
    }

    /// Gets the z-coordinate.
    pub fn z(&self) -> Quantity<U> {
        self.xyz[2]
    }

    pub fn sub(&self, other: &Self) -> Self
    where
        U: std::cmp::PartialEq + std::fmt::Debug,
    {
        debug_assert!(
            self.center_params == other.center_params,
            "Cannot subtract vectors with different center parameters"
        );
        Self::from_vec3(self.center_params.clone(), self.as_vec3() - other.as_vec3())
    }

    /// Calculates the Euclidean distance with respect to the ReferenceCenter.
    ///
    /// # Returns
    /// The distance from the ReferenceCenter in units of U.
    pub fn distance(&self) -> Quantity<U> {
        let d = math::geometry::magnitude(self.x().value(), self.y().value(), self.z().value());
        Quantity::new(d)
    }

    /// Computes the Euclidean distance to another vector of the same type.
    pub fn distance_to(&self, other: &Self) -> Quantity<U>
    where
        U: std::cmp::PartialEq + std::fmt::Debug,
    {
        self.sub(other).distance()
    }
}

// =============================================================================
// Convenience constructors for centers with Params = ()
// =============================================================================

impl<C, F, U> Vector<C, F, U>
where
    C: centers::ReferenceCenter<Params = ()>,
    F: frames::ReferenceFrame,
    U: Unit,
{
    /// Creates a new Cartesian coordinate for centers with `Params = ()`.
    ///
    /// This is a convenience constructor that doesn't require passing `()` explicitly.
    /// Use this for `Barycentric`, `Heliocentric`, `Geocentric`, etc.
    ///
    /// # Arguments
    /// - `x`, `y`, `z`: Coordinate values (converted to `Quantity<U>`).
    ///
    /// # Example
    /// ```rust
    /// use siderust::coordinates::cartesian::position::Ecliptic;
    /// use qtty::AstronomicalUnit;
    ///
    /// let pos = Ecliptic::<AstronomicalUnit>::new(1.0, 0.0, 0.0);
    /// ```
    pub fn new<T>(x: T, y: T, z: T) -> Self
    where
        T: Into<Quantity<U>>,
    {
        Self::new_with_params((), x, y, z)
    }

    /// Creates a new Cartesian coordinate from a nalgebra Vector3 for centers with `Params = ()`.
    ///
    /// This is a convenience constructor that doesn't require passing `()` explicitly.
    pub fn from_vec3_origin(vec3: Vector3<Quantity<U>>) -> Self {
        Self::from_vec3((), vec3)
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
        debug_assert!(
            self.center_params == other.center_params,
            "Cannot add vectors with different center parameters"
        );
        Self::from_vec3(self.center_params.clone(), self.xyz + other.xyz)
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
