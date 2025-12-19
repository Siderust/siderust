//! # Cartesian Coordinates Module
//!
//! This module provides Cartesian coordinate types for astronomical calculations:
//!
//! - [`Position`]: Affine points with a center and frame (can be center-transformed)
//! - [`Direction`]: Free unit vectors with only a frame (frame-transformable only)
//! - [`Velocity`]: Free velocity vectors with only a frame (frame-transformable only)
//!
//! ## Mathematical Model
//!
//! Positions are points in affine space - they have a meaningful origin (center).
//! Directions and velocities are free vectors - they are translation-invariant
//! and do not have a center.
//!
//! ## Line of Sight
//!
//! To compute observer-dependent directions (line of sight to a target), use
//! [`line_of_sight`] instead of attempting to transform a direction between centers.

pub mod direction;
pub mod position;
pub mod vector;
pub mod velocity;

pub use direction::Direction;
pub use position::Position;
pub use vector::*;
pub use velocity::Velocity;

use crate::coordinates::centers::ReferenceCenter;
use crate::coordinates::frames::ReferenceFrame;
use crate::coordinates::math;
use crate::coordinates::spherical::direction::DirectionUnit;
use qtty::{LengthUnit, Quantity};

/// Computes the line-of-sight direction from an observer to a target.
///
/// This is the mathematically correct way to obtain an observer-dependent
/// direction. The result is a free direction (unit vector) pointing from
/// the observer position toward the target position.
///
/// # Requirements
///
/// Both positions must be in the **same center and frame**. If they are not,
/// convert them first using center and/or frame transformations.
///
/// # Mathematical Definition
///
/// ```text
/// direction = normalize(target - observer)
/// ```
///
/// # Example
///
/// ```rust
/// use siderust::coordinates::cartesian::{line_of_sight, Position, Direction};
/// use siderust::coordinates::centers::Geocentric;
/// use siderust::coordinates::frames::Equatorial;
/// use qtty::*;
///
/// let observer = Position::<Geocentric, Equatorial, AstronomicalUnit>::new(0.0, 0.0, 0.0);
/// let target = Position::<Geocentric, Equatorial, AstronomicalUnit>::new(1.0, 1.0, 0.0);
///
/// let los: Direction<Equatorial> = line_of_sight(&observer, &target);
/// // los points from observer toward target
/// ```
///
/// # Panics
///
/// Panics if the observer and target positions are identical (zero-length vector).
pub fn line_of_sight<C, F, U>(
    observer: &Position<C, F, U>,
    target: &Position<C, F, U>,
) -> Direction<F>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
    U: LengthUnit,
{
    let diff = target.sub(observer);

    // Compute normalized direction components using math module
    let (x, y, z) = math::geometry::normalize(
        diff.x().value(),
        diff.y().value(),
        diff.z().value(),
    );

    Direction::<F>::from_vec3(nalgebra::Vector3::new(
        Quantity::<DirectionUnit>::new(x),
        Quantity::<DirectionUnit>::new(y),
        Quantity::<DirectionUnit>::new(z),
    ))
}

/// Computes the line-of-sight direction and distance from an observer to a target.
///
/// Returns both the direction and the distance (magnitude of the separation vector).
/// This is useful when you need both the pointing direction and the range.
///
/// # Example
///
/// ```rust
/// use siderust::coordinates::cartesian::{line_of_sight_with_distance, Position};
/// use siderust::coordinates::centers::Geocentric;
/// use siderust::coordinates::frames::Equatorial;
/// use qtty::*;
///
/// let observer = Position::<Geocentric, Equatorial, AstronomicalUnit>::new(0.0, 0.0, 0.0);
/// let target = Position::<Geocentric, Equatorial, AstronomicalUnit>::new(1.0, 0.0, 0.0);
///
/// let (dir, dist) = line_of_sight_with_distance(&observer, &target);
/// assert!((dist.value() - 1.0).abs() < 1e-10);
/// ```
pub fn line_of_sight_with_distance<C, F, U>(
    observer: &Position<C, F, U>,
    target: &Position<C, F, U>,
) -> (Direction<F>, Quantity<U>)
where
    C: ReferenceCenter,
    F: ReferenceFrame,
    U: LengthUnit,
{
    let diff = target.sub(observer);
    let d = diff.distance();

    // Compute normalized direction components using math module
    let (x, y, z) = math::geometry::normalize(
        diff.x().value(),
        diff.y().value(),
        diff.z().value(),
    );

    let dir = Direction::<F>::from_vec3(nalgebra::Vector3::new(
        Quantity::<DirectionUnit>::new(x),
        Quantity::<DirectionUnit>::new(y),
        Quantity::<DirectionUnit>::new(z),
    ));

    (dir, d)
}
