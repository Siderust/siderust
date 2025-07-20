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
//! use siderust::coordinates::spherical::Direction;
//! use siderust::coordinates::centers::Heliocentric;
//! use siderust::coordinates::frames::Ecliptic;
//! use siderust::units::Degrees;
//!
//! // Create a heliocentric ecliptic spherical Direction
//! let sph = Direction::<Heliocentric, Ecliptic>::new(Degrees::new(45.0), Degrees::new(7.0));
//! println!("θ: {}, φ: {}", sph.polar, sph.azimuth);
//! ```
//!
//! ## Type Aliases
//! - [`Direction<C, F>`]: Spherical direction (distance is typically `None`).
//!
//! ## Type Aliases
//! - [`Direction<C, F>`]: Spherical direction (angles only, vector is unitary).
//! - [`Ecliptic`]: Heliocentric ecliptic direction (longitude, latitude).
//! - [`Equatorial`]: Geocentric equatorial direction (right ascension, declination).
//! - [`Horizontal`]: Topocentric horizontal direction (azimuth, altitude).
//! - [`ICRS`]: Barycentric ICRS direction (right ascension, declination).
//! - [`HCRS`]: Heliocentric ICRS direction (right ascension, declination).
//! - [`GCRS`]: Geocentric ICRS direction (right ascension, declination).
//! - [`Geographic`]: Geocentric ECEF

use super::SphericalCoord;
use crate::units::{Quantity, LengthUnit, Unitless, Degrees};
use crate::coordinates::{
    centers, frames,
    centers::ReferenceCenter,
    frames::ReferenceFrame,
};

pub type Direction<C, F> = SphericalCoord<C, F, Unitless>;
pub type Ecliptic   = Direction<centers::Heliocentric, frames::Ecliptic>;   // L (l), B (b)
pub type Equatorial = Direction<centers::Geocentric,   frames::Equatorial>; // Dec (δ), RA (α)
pub type Horizontal = Direction<centers::Topocentric,  frames::Horizontal>; // Alt (α), Az (θ)
pub type ICRS       = Direction<centers::Barycentric,  frames::ICRS>; // Dec (δ), RA (α)
pub type HCRS       = Direction<centers::Heliocentric, frames::ICRS>; // Dec (δ), RA (α)
pub type GCRS       = Direction<centers::Geocentric,   frames::ICRS>; // Dec (δ), RA (α)
pub type Geographic = Direction<centers::Geocentric,   frames::ECEF>;  //Latitude (φ),Longitude (λ)

impl<C, F> Direction<C, F>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
{

    pub const fn from_degrees(polar: f64, azimuth: f64) -> Self {
        Self::new_spherical_coord(Degrees::new(polar), Degrees::new(azimuth), Quantity::<Unitless>::new(1.0))
    }

    pub fn position<U: LengthUnit>(&self, magnitude: Quantity<U>) -> super::Position<C, F, U> {
        super::Position::new_spherical_coord(self.polar, self.azimuth, magnitude)
    }
}

impl<C, F> std::fmt::Display for Direction<C, F>
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


#[cfg(test)]
mod tests {
    use super::*;
    use crate::units::Degrees;

    #[test]
    fn creates_valid_spherical_position() {
        let polar = Degrees::new(45.0);
        let azimuth = Degrees::new(90.0);

        let coord = ICRS::new(polar, azimuth);

        assert_eq!(coord.ra().value(), 45.0);
        assert_eq!(coord.dec().value(), 90.0);
        assert_eq!(coord.distance.value(), 1.0);
    }

    #[test]
    fn displays_coordinate_as_string_correctly() {
        let coord = ICRS::new(
            Degrees::new(30.0),
            Degrees::new(60.0),
        );

        let output = coord.to_string();

        assert!(output.contains("θ: 60"), "Missing polar angle");
        assert!(output.contains("φ: 30"), "Missing azimuth");
    }

    #[test]
    fn maintains_high_precision_on_values() {
        let polar = Degrees::new(90.654_321);
        let azimuth = Degrees::new(45.123_456);

        let coord = ICRS::new(polar, azimuth);

        assert!((coord.ra().value() - 90.654_321).abs() < 1e-6);
        assert!((coord.dec().value() - 45.123_456).abs() < 1e-6);
    }
}
