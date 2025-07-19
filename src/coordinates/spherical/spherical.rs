//! # Spherical Coordinates
//!
//! This module defines the generic [`SphericalCoord<C, F>`] type for representing positions or directions
//! in spherical coordinates, parameterized by astronomical reference centers and frames for strong type safety.
//!
//! ## Overview
//!
//! - **Generic over Center and Frame:**
//!   - `C`: Reference center (e.g., `Heliocentric`, `Geocentric`, `Barycentric`).
//!   - `F`: Reference frame (e.g., `ICRS`, `Ecliptic`, `Equatorial`).
//! - **Type Safety:** Operations are only allowed between coordinates with matching type parameters.
//! - **Units:** Angles are stored as [`Degrees`]; distance is optional and typically in AstronomicalUnits or parsecs (see context).
//! - **Conversions:** Seamless conversion to and from [`Vector`] via `From`/`Into`.
//! - **Operations:** Compute Euclidean distance and angular separation between coordinates.
//!
//! ## Example
//! ```rust
//! use siderust::coordinates::spherical::Position;
//! use siderust::coordinates::centers::Heliocentric;
//! use siderust::coordinates::frames::Ecliptic;
//! use siderust::units::{Degrees, AstronomicalUnit};
//!
//! // Create a heliocentric ecliptic spherical position
//! let sph = Position::<Heliocentric, Ecliptic, AstronomicalUnit>::new(Degrees::new(45.0), Degrees::new(7.0), 1.0);
//! println!("θ: {}, φ: {}, r: {:?}", sph.polar, sph.azimuth, sph.distance);
//! ```
//!
//! ## Type Aliases
//! - [`Position<C, F>`]: Spherical position (with distance).
//! - [`Direction<C, F>`]: Spherical direction (distance is typically `None`).
//!
//! ## Methods
//! - [`new(polar, azimuth, distance)`]: Construct a new coordinate.
//! - [`distance_to(other)`]: Compute Euclidean distance to another coordinate.
//! - [`angular_separation(other)`]: Compute angular separation in degrees.
//!
//! ## Display
//! Implements `Display` for readable output including center, frame, angles, and distance.

use crate::units::*;

use crate::coordinates::{
    cartesian,
    frames::*,
    centers::*,
};
use std::marker::PhantomData;

/// Represents a point in spherical coordinates with a specific reference center and frame.
///
/// # Type Parameters
/// - `C`: The reference center (e.g., `Barycentric`, `Heliocentric`).
/// - `F`: The reference frame (e.g., `ICRS`, `Ecliptic`).
#[derive(Debug, Clone, Copy)]
pub struct SphericalCoord<
    C: ReferenceCenter,
    F: ReferenceFrame,
    U: Unit,
> {
    pub polar: Degrees,        // θ (polar/latitude/declination)
    pub azimuth: Degrees,      // φ (azimuth/longitude/right ascension)
    pub distance: Quantity<U>, // Optional distance (AstronomicalUnits, parsec, etc.)

    _center: PhantomData<C>,
    _frame: PhantomData<F>,
}

// filepath: src/coordinates/spherical/spherical.rs
impl<C, F, U> SphericalCoord<C, F, U>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
    U: Unit,
{
    pub const fn new_spherical_coord(polar: Degrees, azimuth: Degrees, distance: Quantity<U>) -> Self {
        Self {
            polar,
            azimuth,
            distance,
            _center: PhantomData,
            _frame: PhantomData,
        }
    }
}

impl<C, F, U> SphericalCoord<C, F, U>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
    U: LengthUnit,
    cartesian::Vector<C, F, U>: for<'a> From<&'a Self>,
{
    /// Calculates the Euclidean distance to another spherical coordinate.
    ///
    /// # Arguments
    /// - `other`: The other spherical coordinate.
    ///
    /// # Returns
    /// The distance to the other coordinate.
    pub fn distance_to(&self, other: &Self) -> Quantity<U>
    where
        U: std::cmp::PartialEq + std::fmt::Debug
    {
        self.to_cartesian()
            .distance_to(&other.to_cartesian())
    }

    /// Calculates the angular separation between this coordinate and another.
    ///
    /// # Arguments
    /// - `other`: The other spherical coordinate.
    ///
    /// # Returns
    /// The angular separation in degrees.
    pub fn angular_separation(&self, other: Self) -> Degrees {
        let az1 = self.azimuth.to::<Radian>();
        let po1 = self.polar.to::<Radian>();
        let az2 = other.azimuth.to::<Radian>();
        let po2 = other.polar.to::<Radian>();

        let x = (po1.cos() * po2.sin()) - (po1.sin() * po2.cos() * (az2 - az1).cos());
        let y = po2.cos() * (az2 - az1).sin();
        let z = (po1.sin() * po2.sin()) + (po1.cos() * po2.cos() * (az2 - az1).cos());

        let angle_rad = (x * x + y * y).sqrt().atan2(z);
        Radians::new(angle_rad).to::<Degree>()
    }
}

impl<C, F, U> std::fmt::Display for SphericalCoord<C, F, U>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
    U: LengthUnit,
    Quantity<U>: std::fmt::Display
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Center: {}, Frame: {}, θ: {:.6}, φ: {:.6}, r: {:.6}",
            C::center_name(),
            F::frame_name(),
            self.polar, self.azimuth, self.distance
        )
    }
}

impl<C, F> std::fmt::Display for SphericalCoord<C, F, Unitless>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Center: {}, Frame: {}, θ: {:.6}, φ: {:.6}",
            C::center_name(),
            F::frame_name(),
            self.polar, self.azimuth
        )
    }
}
