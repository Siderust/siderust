// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Coordinate Transformations Module
//!
//! This module provides a unified and extensible framework for transforming astronomical coordinates
//! between different reference centers (e.g., Barycentric, Heliocentric, Geocentric, Topocentric)
//! and reference frames (e.g., EclipticMeanJ2000, EquatorialMeanJ2000, ICRS, Horizontal).
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
//! use siderust::coordinates::frames::EquatorialMeanJ2000;
//! use qtty::*;
//!
//! let observer =
//!     Position::<Geocentric, EquatorialMeanJ2000, AstronomicalUnit>::new(0.0, 0.0, 0.0);
//! let target =
//!     Position::<Geocentric, EquatorialMeanJ2000, AstronomicalUnit>::new(1.0, 1.0, 1.0);
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
//! use siderust::time::JulianDate;
//!
//! let cart_eq = Position::<Geocentric, EquatorialMeanJ2000, AstronomicalUnit>::new(1.0, 2.0, 3.0);
//! let jd = JulianDate::J2000;
//! // Transform to Geocentric EclipticMeanJ2000 coordinates (frame transform)
//! let cart_geo_ecl: Position<Geocentric, EclipticMeanJ2000, AstronomicalUnit> = cart_eq.to_frame();
//! // Transform to Heliocentric EclipticMeanJ2000 coordinates (center transform)
//! let cart_helio_ecl: Position<Heliocentric, EclipticMeanJ2000, AstronomicalUnit> = cart_geo_ecl.transform(jd);
//! ```
//!
//! ## Related Modules
//!
//! - [`centers`]: Transformations between reference centers (positions only).
//! - [`crate::coordinates::frames`]: Transformations between reference frames (all coordinate types).
//! - [`context`]: Astronomical context for transformation configuration.
//! - [`providers`]: Provider traits for computing time-dependent operators.
//! - [`ext`]: Extension traits for ergonomic method-style transforms.

pub mod centers;
pub mod context;
pub mod ecliptic_of_date;
pub mod ext;
mod frames;
pub mod horizontal;
pub mod providers;
mod to_cartesian;
mod to_spherical;

pub use centers::{to_topocentric_with_ctx, IntoTransformArgs, TransformCenter};
pub use context::AstroContext;
pub use ecliptic_of_date::{FromEclipticTrueOfDate, ToEclipticTrueOfDate};
pub use ext::{DirectionAstroExt, PositionAstroExt, SphericalDirectionAstroExt, VectorAstroExt};
pub use ext::{UsingEngine, WithEngine};
pub use frames::TransformFrame;
pub use horizontal::{
    FromHorizontal, ToHorizontal, TopocentricEquatorialExt, TopocentricHorizontalExt,
};
pub use providers::{center_shift, frame_rotation, CenterShiftProvider, FrameRotationProvider};

use crate::coordinates::{
    cartesian, cartesian::Position, centers::ReferenceCenter, frames::MutableFrame, spherical,
};
use crate::time::JulianDate;
use affn::Rotation3;
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
    fn transform(&self, jd: crate::time::JulianDate) -> Coord;
}

/// Blanket implementation for Position transformations (center + frame changes).
///
/// Applies two transformations:
/// 1. Frame transformation (within the same center)
/// 2. Center transformation (within the rotated frame)
///
/// Restricted to standard centers (`Params = ()`) because combined transforms
/// for Bodycentric or Topocentric require explicit parameter values.
/// Use `to_center` directly for those cases.
impl<C1, C2, F1, F2, U> Transform<Position<C2, F2, U>> for Position<C1, F1, U>
where
    C1: ReferenceCenter<Params = ()>,
    C2: ReferenceCenter<Params = ()>,
    Position<C1, F2, U>: TransformCenter<C2, F2, U>,
    (): FrameRotationProvider<F1, F2>,
    F1: MutableFrame,
    F2: MutableFrame,
    U: LengthUnit,
{
    fn transform(&self, jd: JulianDate) -> Position<C2, F2, U> {
        // Apply the frame rotation at the requested epoch, then shift centers.
        let rot: Rotation3 = frame_rotation::<F1, F2>(jd, &AstroContext::default());
        let [x, y, z] = rot * [self.x(), self.y(), self.z()];
        let rotated = {
            self.center_params();
            Position::<C1, F2, U>::from_vec3((), nalgebra::Vector3::new(x, y, z))
        };
        rotated.to_center(jd)
    }
}

/// Blanket implementation to allow chaining two consecutive `From` operations.
///
/// This implementation allows converting a [`spherical::Position`] from one
/// reference center and frame (`C1`, `F1`) to another (`C2`, `F2`) by applying two
/// transformations:
/// 1. Frame transformation (within the same center)
/// 2. Center transformation (within the new frame)
impl<C1, C2, F1, F2, U> Transform<spherical::Position<C2, F2, U>> for spherical::Position<C1, F1, U>
where
    cartesian::Position<C1, F1, U>: Transform<cartesian::Position<C2, F2, U>>,
    C1: ReferenceCenter,
    C2: ReferenceCenter,
    F1: MutableFrame,
    F2: MutableFrame,
    U: LengthUnit,
{
    fn transform(&self, jd: JulianDate) -> spherical::Position<C2, F2, U> {
        let rotated: cartesian::Position<C2, F2, U> = self.to_cartesian().transform(jd);
        spherical::Position::from_cartesian(&rotated)
    }
}

// Note: Frame/center transformations using `From` trait were removed because they
// violate Rust's orphan rules when using affn types directly.
//
// Use the extension traits instead:
// - `position.to_frame::<NewFrame>(&jd, &ctx)` - for frame transforms
// - `position.to_center::<NewCenter>(&jd, &ctx)` - for center transforms
// - `position.to::<NewCenter, NewFrame>(&jd, &ctx)` - for combined transforms
//
// Or use the `Transform` trait:
// - `position.transform(jd)` - uses type inference for target
