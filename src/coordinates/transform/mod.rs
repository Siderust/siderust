//! # Coordinate Transformations Module
//!
//! This module provides a unified and extensible framework for transforming astronomical coordinates
//! between different reference centers (e.g., Barycentric, Heliocentric, Geocentric, Topocentric)
//! and reference frames (e.g., Ecliptic, Equatorial, ICRS, Horizontal).
//!
//! ## Core Concepts
//!
//! - **Transform Trait**: The central abstraction is the [`Transform`] trait, which defines a method
//!   for converting a coordinate of one type into another, possibly using additional context such as
//!   the Julian Date (for time-dependent transformations).
//!
//! - **Center vs Frame Transforms**:
//!   - **Center transforms** (translations) apply only to **positions**. Changing a center
//!     moves the origin from which positions are measured.
//!   - **Frame transforms** (rotations) apply to positions, directions, and velocities.
//!
//! ## Mathematical Foundations
//!
//! - **Positions** are affine points - they can undergo both center and frame transforms.
//! - **Directions** are free vectors (unit vectors) - they can only undergo frame transforms.
//! - **Velocities** are free vectors - they can only undergo frame transforms.
//!
//! Attempting to center-transform a direction or velocity is mathematically undefined and
//! prevented at the type level.
//!
//! ## Observer-Dependent Directions (Line of Sight)
//!
//! To compute the direction to a target as seen from an observer, use the
//! [`line_of_sight`](crate::coordinates::cartesian::line_of_sight) function:
//!
//! ```rust
//! use siderust::coordinates::cartesian::{line_of_sight, Position};
//! use siderust::coordinates::centers::Geocentric;
//! use siderust::coordinates::frames::Equatorial;
//! use qtty::*;
//!
//! let observer = Position::<Geocentric, Equatorial, AstronomicalUnit>::new(0.0, 0.0, 0.0);
//! let target = Position::<Geocentric, Equatorial, AstronomicalUnit>::new(1.0, 1.0, 1.0);
//!
//! let direction = line_of_sight(&observer, &target);
//! ```
//!
//! ## Usage Example
//!
//! ```rust
//! use siderust::coordinates::{cartesian::Position, frames::*, centers::*};
//! use siderust::coordinates::transform::{Transform, TransformFrame};
//! use qtty::AstronomicalUnit;
//! use siderust::astro::JulianDate;
//!
//! let cart_eq = Position::<Geocentric, Equatorial, AstronomicalUnit>::new(1.0, 2.0, 3.0);
//! let jd = JulianDate::J2000;
//! // Transform to Geocentric Ecliptic coordinates (frame transform)
//! let cart_geo_ecl: Position<Geocentric, Ecliptic, AstronomicalUnit> = cart_eq.to_frame();
//! // Transform to Heliocentric Ecliptic coordinates (center transform)
//! let cart_helio_ecl: Position<Heliocentric, Ecliptic, AstronomicalUnit> = cart_geo_ecl.transform(jd);
//! ```
//!
//! ## Related Modules
//!
//! - [`centers`]: Transformations between reference centers (positions only).
//! - [`frames`]: Transformations between reference frames (all coordinate types).

mod centers;
mod frames;
mod to_cartesian;
mod to_spherical;

pub use centers::TransformCenter;
pub use frames::TransformFrame;

use crate::astro::JulianDate;
use crate::coordinates::{
    cartesian, cartesian::Vector, centers::ReferenceCenter, frames::MutableFrame, spherical,
};
use qtty::LengthUnit;

/// Trait for transforming coordinates between different centers and/or frames.
///
/// This trait is primarily used for **position** transformations that may involve
/// both center changes (translations) and frame changes (rotations).
pub trait Transform<Coord> {
    /// Transform this coordinate to a different center and/or frame.
    ///
    /// # Arguments
    ///
    /// - `jd`: The Julian Date at which to perform the transformation.
    fn transform(&self, jd: crate::astro::JulianDate) -> Coord;
}

/// Blanket implementation for Position transformations (center + frame changes).
///
/// This implementation allows converting a [`Vector`] in from one
/// reference center and frame (`C1`, `F1`) to another (`C2`, `F2`) by applying two
/// transformations:
/// 1. Frame transformation (within the same center)
/// 2. Center transformation (within the new frame)
impl<C1, C2, F1, F2, U> Transform<Vector<C2, F2, U>> for Vector<C1, F1, U>
where
    Vector<C1, F1, U>: TransformFrame<Vector<C1, F2, U>>,
    Vector<C1, F2, U>: TransformCenter<Vector<C2, F2, U>>,
    C1: ReferenceCenter,
    C2: ReferenceCenter,
    F1: MutableFrame,
    F2: MutableFrame,
    U: LengthUnit,
{
    fn transform(&self, jd: JulianDate) -> Vector<C2, F2, U> {
        self.to_frame().to_center(jd)
    }
}

/// Blanket implementation to allow chaining two consecutive `From` operations.
///
/// This implementation allows converting a [`spherical::Position`] from one
/// reference center and frame (`C1`, `F1`) to another (`C2`, `F2`) by applying two
/// transformations:
/// 1. Frame transformation (within the same center)
/// 2. Center transformation (within the new frame)
impl<C1, C2, F1, F2, U> Transform<spherical::Position<C2, F2, U>>
    for spherical::Position<C1, F1, U>
where
    cartesian::Vector<C1, F1, U>: Transform<cartesian::Vector<C2, F2, U>>,
    C1: ReferenceCenter,
    C2: ReferenceCenter,
    F1: MutableFrame,
    F2: MutableFrame,
    U: LengthUnit,
{
    fn transform(&self, jd: JulianDate) -> spherical::Position<C2, F2, U> {
        self.to_cartesian().transform(jd).to_spherical()
    }
}

/// Blanket implementation to allow chaining two consecutive `From` operations.
///
/// This implementation allows converting a [`Vector`] in from one
/// reference center and frame (`C1`, `F1`) to another (`C2`, `F2`) by applying two
/// transformations:
/// 1. Frame transformation (within the same center)
/// 2. Center transformation (within the new frame)
impl<C1, F1, C2, F2, U> From<&cartesian::Vector<C1, F1, U>> for cartesian::Vector<C2, F2, U>
where
    cartesian::Vector<C1, F1, U>: TransformFrame<cartesian::Vector<C1, F2, U>>,
    cartesian::Vector<C1, F2, U>: Transform<cartesian::Vector<C2, F2, U>>,
    C1: ReferenceCenter,
    C2: ReferenceCenter,
    F1: MutableFrame,
    F2: MutableFrame,
    U: LengthUnit,
{
    fn from(orig: &cartesian::Vector<C1, F1, U>) -> Self {
        orig.to_frame().transform(JulianDate::J2000)
    }
}

/// Blanket implementation for transforming [`spherical::Position`],
/// involving frame and center changes. Internally uses Cartesian conversions.
///
/// The transformation follows these steps:
/// 1. Convert spherical coordinates to Cartesian.
/// 2. Apply frame transformation.
/// 3. Apply center transformation.
/// 4. Convert back to spherical coordinates.
impl<C1, F1, C2, F2, U> From<&spherical::Position<C1, F1, U>> for spherical::Position<C2, F2, U>
where
    Vector<C1, F1, U>: TransformFrame<Vector<C1, F2, U>>, // transform frame
    Vector<C1, F2, U>: Transform<Vector<C2, F2, U>>,      // transform center
    C1: ReferenceCenter,
    C2: ReferenceCenter,
    F1: MutableFrame,
    F2: MutableFrame,
    U: LengthUnit,
{
    fn from(orig: &spherical::Position<C1, F1, U>) -> Self {
        // Step 1: Convert spherical to Cartesian
        // Step 2: Transform to new frame
        // Step 3: Transform to new center
        // Step 4: Convert back to spherical
        orig.to_cartesian()
            .to_frame()
            .transform(JulianDate::J2000)
            .to_spherical()
    }
}
