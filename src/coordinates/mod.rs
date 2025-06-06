//! # Coordinates Module
//!
//! This module defines strongly typed spherical and cartesian coordinate systems used in astronomy.
//! The coordinate systems are implemented using Rust's type system with phantom types
//! to enforce compile-time safety. These phantom types represent the **frame** and **center** of
//! the coordinate system, ensuring that operations between incompatible coordinate systems are
//! disallowed unless explicitly converted.
//!
//! ## Key Concepts
//! - **SphericalCoord<Center, Frame>**: Represents a spherical coordinate system with a specific
//!   reference center (e.g., `Heliocentric`, `Geocentric`) and frame (e.g., `Ecliptic`, `Equatorial`).
//! - **CartesianCoord<Center, Frame>**: Represents a cartesian coordinate system with the same
//!   reference center and frame as `SphericalCoord`. It is mathematically equivalent to spherical
//!   coordinates but expressed in cartesian form (x, y, z).
//! - **Phantom Types**: The `Center` and `Frame` types are phantom types that provide strong typing
//!   without adding runtime overhead.
//! - **Type Safety**: By combining `Center` and `Frame`, the library ensures that only valid
//!   operations are allowed between compatible coordinate systems.
//! - **Conversions**: Conversion between `SphericalCoord` and `CartesianCoord`, as well as between
//!   different frames and centers, is supported using the `From`/`Into` traits. This allows seamless
//!   transformations between coordinate systems while maintaining type safety.
//!
//! ## Supported Coordinate Systems
//! - **Equatorial**: Based on the celestial equator, commonly used in star catalogs.
//! - **Ecliptic**: Aligned with the plane of Earth's orbit around the Sun.
//! - **Horizontal**: Local horizon-based coordinates for observers on Earth.
//! - **ICRS**: International Celestial Reference System, a standard celestial reference frame.
//! - **ECEF**: Earth-Centered Earth-Fixed, used for geographic coordinates.
//!
//! ## Example
//! ```rust
//! use siderust::coordinates::spherical::SphericalCoord;
//! use siderust::coordinates::cartesian::CartesianCoord;
//! use siderust::coordinates::centers::Heliocentric;
//! use siderust::coordinates::frames::Ecliptic;
//! use siderust::units::Degrees;
//!
//! // Create a spherical coordinate
//! let spherical = SphericalCoord::<Heliocentric, Ecliptic>::new(
//!     Degrees::new(45.0), Degrees::new(7.0), 1.0
//! );
//!
//! // Convert to cartesian coordinates
//! let cartesian: CartesianCoord<Heliocentric, Ecliptic> = (&spherical).into();
//!
//! // Convert back to spherical coordinates
//! let spherical_converted: SphericalCoord<Heliocentric, Ecliptic> = (&cartesian).into();
//!
//! println!("Spherical -> Cartesian -> Spherical: {:?}", spherical_converted);
//! ```
//!
//! ## Submodules
//! Each submodule provides specific implementations and utilities for the respective coordinate system:
//! - **Converters**: Functions and traits for converting between coordinate systems.
//! - **Cartesian**: Core implementation of cartesian coordinates.
//! - **Spherical**: Core implementation of spherical coordinates.
//! - **Frames**: Defines reference frames (e.g., `Ecliptic`, `Equatorial`, `ICRS`).
//! - **Centers**: Defines reference centers (e.g., `Heliocentric`, `Geocentric`).

pub mod transform;
pub mod cartesian;
pub mod spherical;
pub mod frames;
pub mod centers;
pub mod kinds;

//pub use cartesian::*;
//pub use spherical::*;
