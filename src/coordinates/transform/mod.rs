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
//! - **Cartesian and Spherical Coordinates**: The system supports both [`Vector`] and
//!   [`SphericalCoord`] types, parameterized by their reference center and frame. Transformations
//!   can occur between these types, as well as between different centers and frames.
//!
//! - **Reference Centers and Frames**: The transformation system is generic over both the center
//!   (e.g., Heliocentric, Geocentric) and the frame (e.g., Ecliptic, Equatorial), allowing for
//!   flexible and type-safe conversions.
//!
//! ## Design and Implementation
//!
//! - **Trait-based Extensibility**: By implementing the [`Transform`] trait for various combinations
//!   of coordinate types, centers, and frames, the system allows for seamless chaining and
//!   composition of transformations.
//!
//! - **Blanket Implementations**: The module provides blanket implementations for identity
//!   transformations (where the input and output types are the same), as well as for chaining
//!   multiple transformations (e.g., changing both center and frame in sequence).
//!
//! - **Conversion Chaining**: The [`From`] trait is implemented for converting between coordinate
//!   types with different centers and frames by chaining the appropriate [`Transform`] operations.
//!
//! - **Time Dependency**: Many transformations depend on the Julian Date (e.g., due to precession,
//!   nutation, or planetary positions). The trait method always receives a `JulianDate` parameter to
//!   support such cases, even if some transformations are time-independent.
//!
//! ## Usage Example
//!
//! ```rust
//! use siderust::coordinates::{cartesian::Position, frames::*, centers::*};
//! use siderust::coordinates::transform::{Transform, TransformFrame};
//! use siderust::units::AstronomicalUnit;
//! use siderust::astro::JulianDate;
//!
//! let cart_eq = Position::<Geocentric, Equatorial, AstronomicalUnit>::new(1.0, 2.0, 3.0);
//! let jd = JulianDate::J2000;
//! // Transform to Geocentric Ecliptic coordinates
//! let cart_geo_ecl: Position<Geocentric, Ecliptic, AstronomicalUnit> = cart_eq.to_frame();
//! // Transform to Heliocentric Ecliptic coordinates
//! let cart_helio_ecl: Position<Heliocentric, Ecliptic, AstronomicalUnit> = cart_geo_ecl.transform(jd);
//! ```
//!
//! ## Related Modules
//!
//! - [`centers`]: Transformations between reference centers (e.g., barycentric <-> heliocentric).
//! - [`frames`]: Transformations between reference frames (e.g., ecliptic <-> equatorial).
//! - [`to_cartesian`] and [`to_spherical`] — Conversions between Cartesian and Spherical forms.
//! - [`to_horizontal`] — Conversion to topocentric horizontal coordinates.
//!
//! ## Testing
//!
//! Each transformation submodule includes comprehensive tests to ensure correctness and
//! reversibility (where applicable), including edge cases and precision checks.
//!
//! ---
//! 
//! [`Transform`]: trait.Transform.html
//! [`Vector`]: ../struct.Vector.html
//! [`SphericalCoord`]: ../struct.SphericalCoord.html
//! [`centers`]: centers/index.html
//! [`frames`]: frames/index.html
//! [`to_cartesian`]: to_cartesian/index.html
//! [`to_spherical`]: to_spherical/index.html
//! [`to_horizontal`]: to_horizontal/index.html

mod centers;
mod frames;
mod to_cartesian;
mod to_spherical;
mod to_direction;
mod to_horizontal;

pub use frames::TransformFrame;
pub use centers::TransformCenter;

use crate::astro::JulianDate;
use crate::coordinates::{
    centers::ReferenceCenter,
    frames::MutableFrame,
    cartesian,
    cartesian::Vector,
    spherical::SphericalCoord
};
use crate::units::*;

pub trait Transform<Coord> {
    fn transform(&self, jd: crate::astro::JulianDate) -> Coord;
}

// Blanket identity transform for Vector<Center, Frame>
impl<C, F, U> Transform<cartesian::Position<C, F, U>> for cartesian::Position<C, F, U>
where
    C: ReferenceCenter,
    F: MutableFrame,
    U: LengthUnit
{
    fn transform(&self, _jd: crate::astro::JulianDate) -> cartesian::Position<C, F,U> {
        Vector::new(self.x(), self.y(), self.z())
    }
}

/// Blanket implementation to allow chaining two consecutive `From` operations.
///
/// This implementation allows converting a [`Vector`] in from one
/// reference center and frame (`C1`, `F1`) to another (`C2`, `F2`) by applying two 
/// transformations:
/// 1. Frame transformation (within the same center)
/// 2. Center transformation (within the new frame)
impl<C1, F1, C2, F2, U> From<&Vector<C1, F1, U>> for Vector<C2, F2, U>
where
    Vector<C1, F1, U>: TransformFrame<Vector<C1, F2, U>>, // transform frame
    Vector<C1, F2, U>: Transform<Vector<C2, F2, U>>, // transform center
    C1: ReferenceCenter,
    C2: ReferenceCenter,
    F1: MutableFrame,
    F2: MutableFrame,
    U: Unit,
{
    fn from(orig: &Vector<C1, F1, U>) -> Self {
        // Step 1: Transform to new frame, keeping the original center.
        let mid: Vector<C1, F2, U> = orig.to_frame();
        // Step 2: Transform to new center, now using the new frame.
        mid.transform(JulianDate::J2000)
    }
}

/// Blanket implementation for transforming [`SphericalCoord`],
/// involving frame and center changes. Internally uses Cartesian conversions.
///
/// The transformation follows these steps:
/// 1. Convert spherical coordinates to Cartesian.
/// 2. Apply frame transformation.
/// 3. Apply center transformation.
/// 4. Convert back to spherical coordinates.
impl<C1, F1, C2, F2, U> From<&SphericalCoord<C1, F1, U>> for SphericalCoord<C2, F2, U>
where
    Vector<C1, F1, U>: TransformFrame<Vector<C1, F2, U>>, // transform frame
    Vector<C1, F2, U>: Transform<Vector<C2, F2, U>>, // transform center
    Vector<C1, F1, U>: for<'a> From<&'a SphericalCoord<C1, F1, U>>, // to_cartesian
    SphericalCoord<C2, F2, U>: for<'a> From<&'a Vector<C2, F2, U>>, // to_spherical
    C1: ReferenceCenter,
    C2: ReferenceCenter,
    F1: MutableFrame,
    F2: MutableFrame,
    U: Unit,
{
    fn from(orig: &SphericalCoord<C1, F1, U>) -> Self {
        // Step 1: Convert spherical to Cartesian
        let cart: Vector<C1, F1, U> = orig.to_cartesian();
        // Step 2: Transform to new frame
        let cart_mid: Vector<C1, F2, U> = cart.to_frame();
        // Step 3: Transform to new center
        let cart_dest: Vector<C2, F2, U> = cart_mid.transform(JulianDate::J2000);
        // Step 4: Convert back to spherical
        cart_dest.to_spherical()
    }
}
