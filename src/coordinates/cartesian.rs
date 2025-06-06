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

use super::{frames, centers};

use std::marker::PhantomData;
use nalgebra::Vector3;
use std::ops::{Add, Sub, Div, Mul};

pub trait Kind {}
pub struct PositionKind;
pub struct DirectionKind; // ||v|| == 1
impl Kind for PositionKind {}
impl Kind for DirectionKind {}

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

pub type Position<C, F>  = CartesianCoord<C, F, PositionKind>;
pub type Direction<C, F> = CartesianCoord<C, F, DirectionKind>;

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
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self::from_vec3(Vector3::<f64>::new(x, y, z))
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

    /// Calculates the Euclidean distance with respect to the ReferenceCenter.
    ///
    /// # Returns
    /// The distance from the ReferenceCenter in AU.
    pub fn distance_from_origin(&self) -> f64 {
        (self.x().powi(2) + self.y().powi(2) + self.z().powi(2)).sqrt()
    }

    pub fn normalize(&self) -> Self {
        let r = self.distance_from_origin();
        Self::new(self.x() / r, self.y() / r, self.z() / r)
    }

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

// === ICRS-based Cartesian coordinate types ===
pub type ICRS = Position<centers::Barycentric,  frames::ICRS>;
pub type HCRS = Position<centers::Heliocentric, frames::ICRS>;
pub type GCRS = Position<centers::Geocentric,   frames::ICRS>;
pub type TCRS = Position<centers::Topocentric,  frames::ICRS>;

// === Ecliptic frame ===
pub type EclipticBarycentricCartesianCoord  = Position<centers::Barycentric,  frames::Ecliptic>;
pub type EclipticHeliocentricCartesianCoord = Position<centers::Heliocentric, frames::Ecliptic>;
pub type EclipticGeocentricCartesianCoord   = Position<centers::Geocentric,   frames::Ecliptic>;
pub type EclipticTopocentricCartesianCoord  = Position<centers::Topocentric,  frames::Ecliptic>;

// === Equatorial frame ===
pub type EquatorialBarycentricCartesianCoord  = Position<centers::Barycentric,  frames::Equatorial>;
pub type EquatorialHeliocentricCartesianCoord = Position<centers::Heliocentric, frames::Equatorial>;
pub type EquatorialGeocentricCartesianCoord   = Position<centers::Geocentric,   frames::Equatorial>;
pub type EquatorialTopocentricCartesianCoord  = Position<centers::Topocentric,  frames::Equatorial>;

// === Horizontal and Earth-fixed frame ===
pub type HorizontalTopocentricCartesianCoord  = Position<centers::Topocentric,  frames::Horizontal>;

// === Geographic and Earth-fixed frames ===
pub type GeographicCartesianCoord = Position<centers::Geocentric, frames::ECEF>;
