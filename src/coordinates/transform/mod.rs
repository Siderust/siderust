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
//! - **Cartesian and Spherical Coordinates**: The system supports both [`CartesianCoord`] and
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
//!   nutation, or planetary positions). The trait method always receives a `JulianDay` parameter to
//!   support such cases, even if some transformations are time-independent.
//!
//! ## Usage Example
//!
//! ```rust
//! use siderust::coordinates::{CartesianCoord, frames::*, centers::*};
//! use siderust::coordinates::transform::Transform;
//! use siderust::units::JulianDay;
//!
//! let cart_eq = CartesianCoord::<Geocentric, Equatorial>::new(1.0, 2.0, 3.0);
//! let jd = JulianDay::J2000;
//! // Transform to Geocentric Ecliptic coordinates
//! let cart_geo_ecl: CartesianCoord<Geocentric, Ecliptic> = cart_eq.transform(jd);
//! // Transform to Heliocentric Ecliptic coordinates
//! let cart_helio_ecl: CartesianCoord<Heliocentric, Ecliptic> = cart_geo_ecl.transform(jd);
//! ```
//!
//! ## Extending the System
//!
//! To add a new transformation, implement the [`Transform`] trait for the relevant source and
//! destination coordinate types. For example, to add a transformation from a new frame or center,
//! provide an implementation like:
//!
//! ```rust,ignore
//! impl Transform<DestinationType> for SourceType {
//!     fn transform(&self, jd: JulianDay) -> DestinationType {
//!         // ...conversion logic...
//!     }
//! }
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
//! [`CartesianCoord`]: ../struct.CartesianCoord.html
//! [`SphericalCoord`]: ../struct.SphericalCoord.html
//! [`centers`]: centers/index.html
//! [`frames`]: frames/index.html
//! [`to_cartesian`]: to_cartesian/index.html
//! [`to_spherical`]: to_spherical/index.html
//! [`to_horizontal`]: to_horizontal/index.html

pub trait Transform<Coord> {
    fn transform(&self, jd: crate::units::JulianDay) -> Coord;
}

mod centers;
mod frames;
mod to_cartesian;
mod to_spherical;
mod to_horizontal;

use crate::coordinates::cartesian;
use crate::coordinates::{
    centers::ReferenceCenter,
    frames::ReferenceFrame,
    kinds::Kind,
    cartesian::CartesianCoord,
    spherical::SphericalCoord
};
use crate::units::JulianDay;

// Blanket identity transform for CartesianCoord<Center, Frame>
impl<C, F> Transform<cartesian::Position<C, F>> for cartesian::Position<C, F>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
{
    fn transform(&self, _jd: crate::units::JulianDay) -> cartesian::Position<C, F> {
        CartesianCoord::new(self.x(), self.y(), self.z())
    }
}

/// Blanket implementation to allow chaining two consecutive `From` operations.
///
/// This implementation allows converting a [`CartesianCoord`] in from one
/// reference center and frame (`C1`, `F1`) to another (`C2`, `F2`) by applying two 
/// transformations:
/// 1. Frame transformation (within the same center)
/// 2. Center transformation (within the new frame)
impl<C1, F1, C2, F2, K> From<&CartesianCoord<C1, F1, K>> for CartesianCoord<C2, F2, K>
where
    CartesianCoord<C1, F1, K>: Transform<CartesianCoord<C1, F2, K>>, // transform frame
    CartesianCoord<C1, F2, K>: Transform<CartesianCoord<C2, F2, K>>, // transform center
    C1: ReferenceCenter,
    C2: ReferenceCenter,
    F1: ReferenceFrame,
    F2: ReferenceFrame,
    K: Kind,
{
    fn from(orig: &CartesianCoord<C1, F1, K>) -> Self {
        // Step 1: Transform to new frame, keeping the original center.
        let mid: CartesianCoord<C1, F2, K> = orig.transform(JulianDay::J2000);
        // Step 2: Transform to new center, now using the new frame.
        mid.transform(JulianDay::J2000)
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
impl<C1, F1, C2, F2, K> From<&SphericalCoord<C1, F1, K>> for SphericalCoord<C2, F2, K>
where
    CartesianCoord<C1, F1, K>: Transform<CartesianCoord<C1, F2, K>>, // transform frame
    CartesianCoord<C1, F2, K>: Transform<CartesianCoord<C2, F2, K>>, // transform center
    CartesianCoord<C1, F1, K>: for<'a> From<&'a SphericalCoord<C1, F1, K>>,
    SphericalCoord<C2, F2, K>: for<'a> From<&'a CartesianCoord<C2, F2, K>>,
    C1: ReferenceCenter,
    C2: ReferenceCenter,
    F1: ReferenceFrame,
    F2: ReferenceFrame,
    K: Kind,
{
    fn from(orig: &SphericalCoord<C1, F1, K>) -> Self {
        // Step 1: Convert spherical to Cartesian
        let cart: CartesianCoord<C1, F1, K> = orig.to_cartesian();
        // Step 2: Transform to new frame
        let cart_mid: CartesianCoord<C1, F2, K> = cart.transform(JulianDay::J2000);
        // Step 3: Transform to new center
        let cart_dest: CartesianCoord<C2, F2, K> = cart_mid.transform(JulianDay::J2000);
        // Step 4: Convert back to spherical
        cart_dest.to_spherical()
    }
}