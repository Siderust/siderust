//! # Center Transformations
//!
//! This module provides transformations between different astronomical reference centers
//! for **position** types only.
//!
//! ## Mathematical Foundations
//!
//! Center transformations are translations in affine space - they change the origin
//! from which positions are measured. Only affine objects (positions) can undergo
//! center transformations.
//!
//! **Directions and velocities are free vectors and cannot be center-transformed.**
//! Attempting to "change the center" of a direction is mathematically undefined.
//!
//! ## For Observer-Dependent Directions
//!
//! To compute the direction to a target as seen from an observer, use
//! [`line_of_sight`](crate::coordinates::cartesian::line_of_sight) with two positions:
//!
//! ```rust
//! use siderust::coordinates::cartesian::{line_of_sight, Position};
//! use siderust::coordinates::centers::Geocentric;
//! use siderust::coordinates::frames::Equatorial;
//! use qtty::*;
//!
//! let observer = Position::<Geocentric, Equatorial, AstronomicalUnit>::new(0.0, 0.0, 0.0);
//! let target = Position::<Geocentric, Equatorial, AstronomicalUnit>::new(1.0, 1.0, 1.0);
//!
//! let direction = line_of_sight(&observer, &target);
//! ```

pub mod position;

// Re-export extension traits for ergonomic imports
pub use position::{FromBodycentricExt, ToBodycentricExt, ToTopocentricExt};

use crate::astro::JulianDate;
use crate::coordinates::{cartesian::Position, centers::*, frames};
use qtty::LengthUnit;

/// Trait for transforming coordinates from one center to another.
///
/// This trait is only implemented for **position** types. Directions and velocities
/// are free vectors and cannot undergo center transformations.
///
/// # Type Parameters
///
/// - `Coord`: The target coordinate type after transformation.
pub trait TransformCenter<Coord> {
    /// Transform this coordinate to a different center.
    ///
    /// # Arguments
    ///
    /// - `jd`: The Julian Date at which to perform the transformation
    ///   (needed for time-dependent positions like Earth's location).
    fn to_center(&self, jd: crate::astro::JulianDate) -> Coord;
}

/// Identity transformation: a position in center C stays in center C.
impl<C, F, U> TransformCenter<Position<C, F, U>> for Position<C, F, U>
where
    C: ReferenceCenter,
    F: frames::ReferenceFrame,
    U: LengthUnit,
{
    fn to_center(&self, _jd: JulianDate) -> Position<C, F, U> {
        Position::<C, F, U>::from_vec3(self.center_params().clone(), *self.as_vec3())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::astro::JulianDate;
    use crate::bodies::solar_system::Earth;
    use crate::coordinates::cartesian;
    use crate::coordinates::transform::Transform;
    use crate::macros::assert_cartesian_eq;
    use qtty::AstronomicalUnit;

    const EPSILON: f64 = 1e-9;

    #[test]
    fn test_position_barycentric_to_geocentric() {
        let earth_bary = *Earth::vsop87e(JulianDate::J2000).get_position();
        let earth_geo: cartesian::position::Ecliptic<AstronomicalUnit, Geocentric> =
            earth_bary.transform(JulianDate::J2000);
        let expected_earth_geo =
            cartesian::position::Ecliptic::<AstronomicalUnit, Geocentric>::CENTER;
        assert_cartesian_eq!(
            &earth_geo,
            &expected_earth_geo,
            EPSILON,
            "Earth in Geocentric should be at origin. Current: {:?}",
            earth_geo
        );
    }

    #[test]
    fn test_position_heliocentric_to_geocentric() {
        let earth_helio = *Earth::vsop87a(JulianDate::J2000).get_position();
        let earth_geo: cartesian::position::Ecliptic<AstronomicalUnit, Geocentric> =
            earth_helio.transform(JulianDate::J2000);
        let expected = cartesian::position::Ecliptic::<AstronomicalUnit, Geocentric>::CENTER;
        assert_cartesian_eq!(&earth_geo, &expected, EPSILON);
    }
}
