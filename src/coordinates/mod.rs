//! Coordinate types and transformations.
//!
//! The `coordinates` module defines strongly typed Cartesian and spherical
//! representations. Phantom markers encode the reference centre, frame and
//! physical dimension so that mismatched coordinates are rejected at compile
//! time. Positions, directions and velocities therefore carry their semantics
//! directly in the type system.
//!
//! Conversions between frames and centres are provided through the
//! [`transform::Transform`] trait, and both representations implement
//! `From`/`Into` for mutual conversion.
//!
//! ## Example
//! ```rust
//! use siderust::coordinates::{cartesian, spherical};
//! use siderust::coordinates::centers::Heliocentric;
//! use siderust::coordinates::frames::Ecliptic;
//! use siderust::units::Degrees;
//!
//! let sph = spherical::Direction::<Heliocentric, Ecliptic>::new(
//!     Degrees::new(45.0),
//!     Degrees::new(7.0),
//! );
//! let cart: cartesian::Direction<Heliocentric, Ecliptic> = (&sph).into();
//! let sph2: spherical::Direction<Heliocentric, Ecliptic> = (&cart).into();
//! assert_eq!(sph, sph2);
//! ```
//!
//! ## Submodules
//! - [`transform`]: generic transformations between coordinate systems and frames.
//! - [`cartesian`]: Cartesian coordinate types and operations.
//! - [`spherical`]: spherical coordinate types and operations.
//! - [`frames`]: reference-frame marker types (e.g. `Ecliptic`, `Equatorial`, `ICRS`).
//! - [`centers`]: reference-centre marker types (e.g. `Heliocentric`, `Geocentric`).

pub mod transform;
pub mod cartesian;
pub mod spherical;
pub mod frames;
pub mod centers;

