//! # Coordinates Module
//!
//! This module defines strongly typed spherical and cartesian coordinate systems used in astronomy.
//! The coordinate systems are implemented using Rust's type system with phantom types
//! to enforce compile-time safety. These phantom types represent the **frame** and **center** of
//! the coordinate system, ensuring that operations between incompatible coordinate systems are
//! disallowed unless explicitly converted.
//!
//! ## Key Concepts
//! - **Position and Direction Types**: Both spherical and cartesian coordinates are parameterized
//!   by a reference center (e.g., `Heliocentric`, `Geocentric`), a reference frame (e.g., `Ecliptic`, `Equatorial`, `ICRS`),
//!   and a kind (`Position`, `Direction`). This ensures that only compatible coordinates can be used together.
//! - **Phantom Types**: The `Center`, `Frame`, and `Kind` types are zero-cost markers that encode coordinate semantics at compile time.
//! - **Type Safety**: Operations between coordinates are only allowed when their type parameters match, preventing accidental mixing of frames or centers.
//! - **Conversions**: Seamless conversion between spherical and cartesian forms, and between different frames and centers, is provided via `From`/`Into` and the `Transform` trait.
//!
//! ## Supported Reference Frames and Centers
//! - **Frames**: `Equatorial`, `Ecliptic`, `Horizontal`, `ICRS`, `ECEF`
//! - **Centers**: `Heliocentric`, `Geocentric`, `Barycentric`, `Topocentric`
//!
//! ## Example
//! ```rust
//! use siderust::coordinates::spherical;
//! use siderust::coordinates::cartesian;
//! use siderust::coordinates::centers::Heliocentric;
//! use siderust::coordinates::frames::Ecliptic;
//! use siderust::units::Degrees;
//!
//! // Create a heliocentric ecliptic spherical position
//! let spherical = spherical::Position::<Heliocentric, Ecliptic>::new(
//!     Degrees::new(45.0), Degrees::new(7.0), 1.0
//! );
//!
//! // Convert to cartesian coordinates
//! let cartesian: cartesian::Position<Heliocentric, Ecliptic> = (&spherical).into();
//!
//! // Convert back to spherical coordinates
//! let spherical_converted: spherical::Position<Heliocentric, Ecliptic> = (&cartesian).into();
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
//! - **kinds**: Marker types for distinguishing positions and directions.

pub mod transform;
pub mod cartesian;
pub mod spherical;
pub mod frames;
pub mod centers;
pub mod kinds;
