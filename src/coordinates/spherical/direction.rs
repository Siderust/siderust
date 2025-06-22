//! # Spherical Coordinates
//!
//! This module defines the generic [`SphericalCoord<C, F, K>`] type for representing positions or directions
//! in spherical coordinates, parameterized by astronomical reference centers and frames for strong type safety.
//!
//! ## Overview
//!
//! - **Generic over Center, Frame, and Kind:**
//!   - `C`: Reference center (e.g., `Heliocentric`, `Geocentric`, `Barycentric`).
//!   - `F`: Reference frame (e.g., `ICRS`, `Ecliptic`, `Equatorial`).
//!   - `K`: Kind marker (`PositionKind`, `DirectionKind`), enforcing semantic correctness.
//! - **Type Safety:** Operations are only allowed between coordinates with matching type parameters.
//! - **Units:** Angles are stored as [`Degrees`]; distance is optional and typically in AU or parsecs (see context).
//! - **Conversions:** Seamless conversion to and from [`Vector`] via `From`/`Into`.
//! - **Operations:** Compute Euclidean distance and angular separation between coordinates.
//!
//! ## Example
//! ```rust
//! use siderust::coordinates::spherical::Position;
//! use siderust::coordinates::centers::Heliocentric;
//! use siderust::coordinates::frames::Ecliptic;
//! use siderust::units::Degrees;
//!
//! // Create a heliocentric ecliptic spherical position
//! let sph = Position::<Heliocentric, Ecliptic>::new(Degrees::new(45.0), Degrees::new(7.0), 1.0);
//! println!("θ: {}, φ: {}, r: {:?}", sph.polar, sph.azimuth, sph.distance);
//! ```
//!
//! ## Type Aliases
//! - [`Position<C, F>`]: Spherical position (with distance).
//! - [`Direction<C, F>`]: Spherical direction (distance is typically `None`).
//!
//! ## Methods
//! - [`new(polar, azimuth, distance)`]: Construct a new coordinate.
//! - [`from_degrees(polar, azimuth, distance)`]: Construct from primitive values.
//! - [`distance_to(other)`]: Compute Euclidean distance to another coordinate.
//! - [`angular_separation(other)`]: Compute angular separation in degrees.
//!
//! ## Display
//! Implements `Display` for readable output including center, frame, angles, and distance.

use super::SphericalCoord;
use crate::coordinates::{
    centers::ReferenceCenter,
    frames::ReferenceFrame,
    kinds::DirectionKind
};
use crate::units::Unit;

pub type Direction<C, F, U=f64> = SphericalCoord<C, F, U, DirectionKind>;

impl<C, F, U> Direction<C, F, U>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
    U: Unit,
{
    pub fn position(&self, magnitude: U) -> super::Position<C, F, U> {
        Self::new_spherical_coord(self.polar, self.azimuth, Some(magnitude))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::units::Degrees;
    use crate::coordinates::centers::{Barycentric, Geocentric};
    use crate::coordinates::frames::ICRS;

    #[test]
    fn creates_valid_spherical_position() {
        let polar = Degrees::new(45.0);
        let azimuth = Degrees::new(90.0);

        let coord = Direction::<Barycentric, ICRS>::new(polar, azimuth);

        assert_eq!(coord.ra().as_f64(), 45.0);
        assert_eq!(coord.dec().as_f64(), 90.0);
        assert_eq!(coord.distance, None);
    }

    #[test]
    fn displays_coordinate_as_string_correctly() {
        let coord = Direction::<Geocentric, ICRS>::new(
            Degrees::new(30.0),
            Degrees::new(60.0),
        );

        let output = coord.to_string();

        assert!(output.contains("θ: 60"), "Missing polar angle");
        assert!(output.contains("φ: 30"), "Missing azimuth");
        assert!(output.contains("r: NaN"), "Distance should be NaN");
    }

    #[test]
    fn maintains_high_precision_on_values() {
        let polar = Degrees::new(90.654_321);
        let azimuth = Degrees::new(45.123_456);

        let coord = Direction::<Barycentric, ICRS>::new(polar, azimuth);

        assert!((coord.ra().as_f64() - 90.654_321).abs() < 1e-6);
        assert!((coord.dec().as_f64() - 45.123_456).abs() < 1e-6);
        assert_eq!(coord.distance, None);
    }
}
