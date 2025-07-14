//! Defines the generic [`SphericalCoord<Center, Frame>`] type for representing
//! strongly-typed spherical coordinates in astronomy.
//!
//! ## Features
//! - **Generic Type Safety:** The type system enforces correct pairing of reference centers and frames.
//! - **Strongly Typed Angles:** Uses [`crate::units::Degrees`] and [`crate::units::Radians`] for angular components to avoid unit confusion.
//! - **LengthUnit and Angular Separation:** Methods for computing Euclidean distances and angular separations.
//! - **String Representation:** Human-readable formatting for debugging and logging.
//!
//! ## Coordinate Convention
//! The `SphericalCoord<Center, Frame>` type represents a point in 3D space using:
//!
//! - **Polar angle (θ):** `polar` — angle from a reference plane (e.g., latitude, declination, altitude), in degrees.
//! - **Azimuthal angle (φ):** `azimuth` — angle from a reference direction (e.g., longitude, right ascension, azimuth), in degrees.
//! - **Radial distance (r):** `distance` — distance from the reference center, in astronomical units (AstronomicalUnits) or meters (for ECEF).
//!
//! The interpretation of these components depends on the chosen `Frame` and `Center`.
//!
//! ## Type Safety
//! By encoding the reference center and frame as type parameters, this module prevents mixing incompatible coordinate systems at compile time.
//!
//! ## Example Usage
//! ```rust,ignore
//! use siderust::units::Degrees;
//! use siderust::coordinates::spherical::direction::ICRS;
//!
//! let coord = ICRS::new(
//!     Degrees::new(120.0), Degrees::new(45.0)
//! );
//! println!("{}", coord.to_string());
//! ```
//!
//! ## See Also
//! - [`crate::coordinates::cartesian::Vector`] for 3D Cartesian representation.
//! - Frame-specific modules for specialized constructors and accessors.
//!
//! ---
//!
//! Defines the generic `SphericalCoord<center, frame>` type for representing 
//! strongly typed spherical coordinates with various reference centers and frames.
//!
//! This module provides:
//! - The `SphericalCoord` struct, which represents a point in spherical coordinates
//!   with a specific reference center and frame.
//! - Methods for calculating distances, angular separations, and converting to string representations.
//!
//! # Coordinate Components
//! The `SphericalCoord` type uses:
//!
//! - **Polar (θ)** → `polar`: angle from a reference plane (e.g., latitude, declination, altitude), in degrees.
//! - **Azimuth (φ)** → `azimuth`: angle from a reference direction (e.g., longitude, right ascension, azimuth), in degrees.
//! - **Radial distance (r)** → `distance`: distance from the reference center, in AstronomicalUnits.
//!
//! # Example
//! ```rust
//! use siderust::units::Degrees;
//! use siderust::coordinates::spherical::direction::ICRS;
//!
//! let coord = ICRS::new(
//!     Degrees::new(120.0), Degrees::new(45.0)
//! );
//! println!("{}", coord.to_string());
//! ```

#[allow(clippy::all)]
mod spherical;
pub use spherical::*;

pub mod position;
pub use position::Position;

pub mod direction;
pub use direction::Direction;

mod equatorial;
mod ecliptic;
mod horizontal;
mod icrs;
mod ecef;
