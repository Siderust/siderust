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
//!   by a reference center (e.g., `Heliocentric`, `Geocentric`), a reference frame (e.g., `Ecliptic`, `Equatorial`, `ICRS`),
//!   and a measure unit (`Unitless`, `LengthUnit`, `VelocityUnit`). This ensures that only compatible coordinates can be used together.
//! - **Phantom Types**: The `Center`, `Frame` and `Unit`types are zero-cost markers that encode coordinate semantics at compile time.
//! - **Type Safety**: Operations between coordinates are only allowed when their type parameters match, preventing accidental mixing of frames, centers or magnitude.
//! - **Conversions**: Seamless conversion between spherical and cartesian forms, and between different frames and centers, is provided via `From`/`Into` and the `Transform` trait.
//!
//! ## Supported Reference Frames and Centers
//! - **Frames**: `Equatorial`, `Ecliptic`, `Horizontal`, `ICRS`, `ECEF`
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
//! let spherical_converted: spherical::Direction<Ecliptic> = cartesian.to_spherical();
//!
//! println!("Spherical -> Cartesian -> Spherical: {:?}", spherical_converted);
//! ```
//!
//! ## Submodules
//! - **transform**: Generic transformations between coordinate systems and frames.
//! - **cartesian**: Cartesian coordinate types and operations.
//! - **spherical**: Spherical coordinate types and operations.
//! - **frames**: Reference frame marker types (e.g., `Ecliptic`, `Equatorial`, `ICRS`).
//! - **centers**: Reference center marker types (e.g., `Heliocentric`, `Geocentric`).

pub mod cartesian;
pub mod centers;
pub mod frames;
pub mod spherical;
pub mod transform;
