// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Coordinates Module
//!
//! This module defines strongly typed spherical and cartesian coordinate systems used in astronomy.
//! The coordinate systems are implemented using Rust's type system with phantom types
//! to enforce compile-time safety. These phantom types represent the **frame** and **center** of
//! the coordinate system, ensuring that operations between incompatible coordinate systems are
//! disallowed unless explicitly converted. Moreover, thanks to the Unit module we can distinguish
//! the different vector types such as Directions (unitless), Position (Lenght Units) and Velocity
//! (Velocity Units), that enforce the compiler to validate any transformation of coordinates.

//! ## Key Concepts
//! - **Position, Direction and Velocity Types**: Both spherical and cartesian coordinates are parameterized
//!   by a reference center (e.g., `Heliocentric`, `Geocentric`), a reference frame (e.g., `Ecliptic`, `EquatorialMeanJ2000`, `ICRS`),
//!   and a measure unit (`Unitless`, `LengthUnit`, `VelocityUnit`). This ensures that only compatible coordinates can be used together.
//! - **Phantom Types**: The `Center`, `Frame` and `Unit`types are zero-cost markers that encode coordinate semantics at compile time.
//! - **Type Safety**: Operations between coordinates are only allowed when their type parameters match, preventing accidental mixing of frames, centers or magnitude.
//! - **Conversions**: Seamless conversion between spherical and cartesian forms, and between different frames and centers, is provided via `From`/`Into` and the `Transform` trait.
//!
//! ## Module Organization
//!
//! The coordinate system is organized into several modules:
//!
//! ### Core Modules
//!
//! - **Reference frames and centers**: Trait definitions for orientation and origin in [`frames`] and [`centers`]
//! - **Cartesian types**: Vector, Direction, Position, Velocity in [`cartesian`]
//! - **Spherical types**: Direction, Position with astronomical extensions in [`spherical`]
//!
//! ### Additional Modules
//!
//! - **transform**: Generic transformations between coordinate systems and frames
//! - **observation**: Observer-dependent effects like aberration
//! - **horizontal**: Convention conversion helpers for horizontal (alt-az) coordinates
//!
//! The coordinate types are built on top of the `affn` crate (the pure geometry kernel).
//! Astronomy-specific frames, centers, and convenience methods are defined in this module.
//!
//! ## Coordinate Transform Architecture
//!
//! The coordinate system maintains a clean separation of concerns:
//!
//! - **Center transforms** (translations): Apply only to positions. Moving from geocentric to
//!   heliocentric is a pure vector subtraction. No observation effects.
//!
//! - **Frame transforms** (rotations): Apply to positions, directions, and velocities.
//!   Changing from ecliptic to equatorial is a pure rotation matrix.
//!
//! - **Observation transforms** (in [`observation`] module): Observer-dependent effects like
//!   aberration. These require explicit `ObserverState` and produce directions with explicit
//!   observational state (`Astrometric` or `Apparent`).
//!
//! ## Supported Reference Frames and Centers
//! - **Frames**: `EquatorialMeanJ2000`, `EquatorialMeanOfDate`, `EquatorialTrueOfDate`, `Ecliptic`, `Horizontal`, `ICRS`, `ECEF`
//! - **Centers**: `Heliocentric`, `Geocentric`, `Barycentric`, `Topocentric`, `Bodycentric`
//!
//! ## Example
//! ```rust
//! use siderust::coordinates::spherical;
//! use siderust::coordinates::cartesian;
//! use siderust::coordinates::frames::Ecliptic;
//! use qtty::*;
//!
//! // Create an ecliptic spherical direction (frame-only, no center)
//! let spherical = spherical::Direction::<Ecliptic>::new(
//!     45.0 * DEG, 7.0 * DEG
//! );
//!
//! // Convert to cartesian coordinates
//! let cartesian: cartesian::Direction<Ecliptic> = spherical.to_cartesian();
//!
//! // Convert back to spherical coordinates
//! let spherical_converted: spherical::Direction<Ecliptic> =
//!     spherical::Direction::from_cartesian(&cartesian);
//!
//! println!("Spherical -> Cartesian -> Spherical: {:?}", spherical_converted);
//! ```
//!
//! ## Submodules
//! - **frames**: Reference frame definitions (Ecliptic, EquatorialMeanJ2000, ICRS, etc.)
//! - **centers**: Reference center definitions (Heliocentric, Geocentric, etc.)
//! - **cartesian**: Cartesian coordinate types and astronomical type aliases
//! - **spherical**: Spherical coordinate types and astronomical extensions
//! - **transform**: Generic transformations between coordinate systems and frames
//! - **observation**: Observational state types (`Astrometric`, `Apparent`) and aberration
//!
//! ## Prelude
//!
//! For ergonomic imports of extension traits:
//! ```rust,ignore
//! use siderust::coordinates::prelude::*;
//! ```

pub mod cartesian;
pub mod centers;
pub mod frames;
pub mod horizontal;
pub mod observation;
pub mod spherical;
pub mod transform;

/// Prelude module for convenient imports.
///
/// Import this to get access to all coordinate extension traits:
///
/// ```rust
/// use siderust::coordinates::prelude::*;
/// ```
///
/// This includes:
/// - [`DirectionAstroExt`](transform::DirectionAstroExt) - Frame transforms for directions
/// - [`VectorAstroExt`](transform::VectorAstroExt) - Frame transforms for vectors
/// - [`PositionAstroExt`](transform::PositionAstroExt) - Frame and center transforms for positions
/// - [`AstroContext`](transform::AstroContext) - Context for transformations
pub mod prelude {
    pub use super::transform::{AstroContext, DirectionAstroExt, PositionAstroExt, VectorAstroExt};
}
