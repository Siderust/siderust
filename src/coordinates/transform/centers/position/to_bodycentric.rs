//! # Bodycentric Position Transformations
//!
//! This module provides transformations to body-centered coordinate systems.
//! A body-centric coordinate system has its origin at an orbiting celestial body
//! (satellite, planet, moon, comet, asteroid, etc.).
//!
//! ## Overview
//!
//! The transformation computes the position relative to the body by:
//! 1. Computing the body's position at the given Julian date using Keplerian propagation
//! 2. Converting the body's position to match the source coordinate's center and frame
//! 3. Subtracting to get the relative position
//!
//! ## Usage
//!
//! Since body-centric coordinates require runtime parameters (the body's orbit),
//! transformations use an explicit `to_bodycentric()` method rather than the
//! blanket `Transform` trait.
//!
//! ```rust
//! use siderust::coordinates::centers::{Bodycentric, BodycentricParams, Geocentric};
//! use siderust::coordinates::cartesian::Position;
//! use siderust::coordinates::frames;
//! use siderust::astro::orbit::Orbit;
//! use siderust::astro::JulianDate;
//! use qtty::*;
//!
//! // Define satellite orbit
//! let satellite_orbit = Orbit::new(
//!     0.0000426 * AU,       // ~6378 km in AU
//!     0.001,
//!     Degrees::new(51.6),
//!     Degrees::new(0.0),
//!     Degrees::new(0.0),
//!     Degrees::new(0.0),
//!     JulianDate::J2000,
//! );
//!
//! let sat_params = BodycentricParams::geocentric(satellite_orbit);
//!
//! // Transform a geocentric position to satellite-centric
//! let target_geo: Position<Geocentric, frames::Equatorial, AstronomicalUnit> =
//!     Position::new(0.00257, 0.0, 0.0); // Moon distance in AU
//!
//! let target_from_sat = target_geo.to_bodycentric(sat_params, JulianDate::J2000);
//! ```

use crate::astro::JulianDate;
use crate::coordinates::cartesian::position::{Ecliptic, Position};
use crate::coordinates::centers::{
    Barycentric, Bodycentric, BodycentricParams, Geocentric, Heliocentric, OrbitReferenceCenter,
};
use crate::coordinates::frames::MutableFrame;
use crate::coordinates::transform::centers::TransformCenter;
use crate::coordinates::transform::TransformFrame;
use qtty::{AstronomicalUnits, LengthUnit, Quantity};

// =============================================================================
// Geocentric → Bodycentric
// =============================================================================

impl<F: MutableFrame, U: LengthUnit> Position<Geocentric, F, U>
where
    Quantity<U>: From<AstronomicalUnits>,
    Position<Geocentric, F, U>: TransformFrame<Ecliptic<U, Geocentric>>,
    Ecliptic<U, Geocentric>: TransformFrame<Position<Geocentric, F, U>>,
{
    /// Transform to body-centric coordinates.
    ///
    /// This computes the position relative to an orbiting body (satellite, moon, etc.)
    /// at the given Julian date.
    ///
    /// # Arguments
    ///
    /// - `body_params`: The orbital parameters of the body to use as the new center.
    /// - `jd`: The Julian date at which to compute the body's position.
    ///
    /// # Returns
    ///
    /// A position in the body-centric coordinate system.
    ///
    /// # Example
    ///
    /// ```rust
    /// use siderust::coordinates::centers::{BodycentricParams, Geocentric};
    /// use siderust::coordinates::cartesian::Position;
    /// use siderust::coordinates::frames;
    /// use siderust::astro::orbit::Orbit;
    /// use siderust::astro::JulianDate;
    /// use qtty::*;
    ///
    /// let satellite_orbit = Orbit::new(
    ///     0.0000426 * AU,
    ///     0.001,
    ///     Degrees::new(51.6),
    ///     Degrees::new(0.0),
    ///     Degrees::new(0.0),
    ///     Degrees::new(0.0),
    ///     JulianDate::J2000,
    /// );
    ///
    /// let target: Position<Geocentric, frames::Ecliptic, AstronomicalUnit> =
    ///     Position::new(0.00257, 0.0, 0.0);
    ///
    /// let sat_params = BodycentricParams::geocentric(satellite_orbit);
    /// let relative = target.to_bodycentric(sat_params, JulianDate::J2000);
    /// ```
    pub fn to_bodycentric(
        &self,
        body_params: BodycentricParams,
        jd: JulianDate,
    ) -> Position<Bodycentric, F, U> {
        // Get the body's position in ecliptic coordinates at the given time
        let body_ecliptic_au = body_params.orbit.kepler_position(jd);

        // Convert the body's position to our center and frame
        let body_in_our_center: Position<Geocentric, F, U> = match body_params.orbit_center {
            OrbitReferenceCenter::Geocentric => {
                // Body orbits Earth - position is already geocentric
                let body_ecl: Ecliptic<U, Geocentric> = Ecliptic::new(
                    body_ecliptic_au.x(),
                    body_ecliptic_au.y(),
                    body_ecliptic_au.z(),
                );
                body_ecl.to_frame()
            }
            OrbitReferenceCenter::Heliocentric => {
                // Body orbits Sun - need to convert heliocentric to geocentric
                // First create heliocentric ecliptic, then transform to geocentric
                let body_helio_ecl: Ecliptic<U, Heliocentric> = Ecliptic::new(
                    body_ecliptic_au.x(),
                    body_ecliptic_au.y(),
                    body_ecliptic_au.z(),
                );
                let body_geo_ecl: Ecliptic<U, Geocentric> = body_helio_ecl.to_center(jd);
                body_geo_ecl.to_frame()
            }
            OrbitReferenceCenter::Barycentric => {
                // Body orbits barycenter - transform to geocentric
                let body_bary_ecl: Ecliptic<U, Barycentric> = Ecliptic::new(
                    body_ecliptic_au.x(),
                    body_ecliptic_au.y(),
                    body_ecliptic_au.z(),
                );
                let body_geo_ecl: Ecliptic<U, Geocentric> = body_bary_ecl.to_center(jd);
                body_geo_ecl.to_frame()
            }
        };

        // Compute relative position: target - body
        let relative_vec = self.as_vec3() - body_in_our_center.as_vec3();

        Position::<Bodycentric, F, U>::from_vec3(body_params, relative_vec)
    }
}

// =============================================================================
// Heliocentric → Bodycentric
// =============================================================================

impl<F: MutableFrame, U: LengthUnit> Position<Heliocentric, F, U>
where
    Quantity<U>: From<AstronomicalUnits>,
    Position<Heliocentric, F, U>: TransformFrame<Ecliptic<U, Heliocentric>>,
    Ecliptic<U, Heliocentric>: TransformFrame<Position<Heliocentric, F, U>>,
    Ecliptic<U, Geocentric>: TransformCenter<Ecliptic<U, Heliocentric>>,
    Ecliptic<U, Barycentric>: TransformCenter<Ecliptic<U, Heliocentric>>,
{
    /// Transform to body-centric coordinates.
    ///
    /// This computes the position relative to an orbiting body (planet, comet, asteroid, etc.)
    /// at the given Julian date.
    pub fn to_bodycentric(
        &self,
        body_params: BodycentricParams,
        jd: JulianDate,
    ) -> Position<Bodycentric, F, U> {
        // Get the body's position in ecliptic coordinates at the given time
        let body_ecliptic_au = body_params.orbit.kepler_position(jd);

        // Convert the body's position to our center and frame
        let body_in_our_center: Position<Heliocentric, F, U> = match body_params.orbit_center {
            OrbitReferenceCenter::Heliocentric => {
                // Body orbits Sun - position is already heliocentric
                let body_ecl: Ecliptic<U, Heliocentric> = Ecliptic::new(
                    body_ecliptic_au.x(),
                    body_ecliptic_au.y(),
                    body_ecliptic_au.z(),
                );
                body_ecl.to_frame()
            }
            OrbitReferenceCenter::Geocentric => {
                // Body orbits Earth - need to convert geocentric to heliocentric
                let body_geo_ecl: Ecliptic<U, Geocentric> = Ecliptic::new(
                    body_ecliptic_au.x(),
                    body_ecliptic_au.y(),
                    body_ecliptic_au.z(),
                );
                let body_helio_ecl: Ecliptic<U, Heliocentric> = body_geo_ecl.to_center(jd);
                body_helio_ecl.to_frame()
            }
            OrbitReferenceCenter::Barycentric => {
                // Body orbits barycenter - transform to heliocentric
                let body_bary_ecl: Ecliptic<U, Barycentric> = Ecliptic::new(
                    body_ecliptic_au.x(),
                    body_ecliptic_au.y(),
                    body_ecliptic_au.z(),
                );
                let body_helio_ecl: Ecliptic<U, Heliocentric> = body_bary_ecl.to_center(jd);
                body_helio_ecl.to_frame()
            }
        };

        // Compute relative position: target - body
        let relative_vec = self.as_vec3() - body_in_our_center.as_vec3();

        Position::<Bodycentric, F, U>::from_vec3(body_params, relative_vec)
    }
}

// =============================================================================
// Barycentric → Bodycentric
// =============================================================================

impl<F: MutableFrame, U: LengthUnit> Position<Barycentric, F, U>
where
    Quantity<U>: From<AstronomicalUnits>,
    Position<Barycentric, F, U>: TransformFrame<Ecliptic<U, Barycentric>>,
    Ecliptic<U, Barycentric>: TransformFrame<Position<Barycentric, F, U>>,
    Ecliptic<U, Geocentric>: TransformCenter<Ecliptic<U, Barycentric>>,
    Ecliptic<U, Heliocentric>: TransformCenter<Ecliptic<U, Barycentric>>,
{
    /// Transform to body-centric coordinates.
    ///
    /// This computes the position relative to an orbiting body at the given Julian date.
    pub fn to_bodycentric(
        &self,
        body_params: BodycentricParams,
        jd: JulianDate,
    ) -> Position<Bodycentric, F, U> {
        // Get the body's position in ecliptic coordinates at the given time
        let body_ecliptic_au = body_params.orbit.kepler_position(jd);

        // Convert the body's position to our center and frame
        let body_in_our_center: Position<Barycentric, F, U> = match body_params.orbit_center {
            OrbitReferenceCenter::Barycentric => {
                // Body orbits barycenter - position is already barycentric
                let body_ecl: Ecliptic<U, Barycentric> = Ecliptic::new(
                    body_ecliptic_au.x(),
                    body_ecliptic_au.y(),
                    body_ecliptic_au.z(),
                );
                body_ecl.to_frame()
            }
            OrbitReferenceCenter::Heliocentric => {
                // Body orbits Sun - need to convert heliocentric to barycentric
                let body_helio_ecl: Ecliptic<U, Heliocentric> = Ecliptic::new(
                    body_ecliptic_au.x(),
                    body_ecliptic_au.y(),
                    body_ecliptic_au.z(),
                );
                let body_bary_ecl: Ecliptic<U, Barycentric> = body_helio_ecl.to_center(jd);
                body_bary_ecl.to_frame()
            }
            OrbitReferenceCenter::Geocentric => {
                // Body orbits Earth - need to convert geocentric to barycentric
                let body_geo_ecl: Ecliptic<U, Geocentric> = Ecliptic::new(
                    body_ecliptic_au.x(),
                    body_ecliptic_au.y(),
                    body_ecliptic_au.z(),
                );
                let body_bary_ecl: Ecliptic<U, Barycentric> = body_geo_ecl.to_center(jd);
                body_bary_ecl.to_frame()
            }
        };

        // Compute relative position: target - body
        let relative_vec = self.as_vec3() - body_in_our_center.as_vec3();

        Position::<Bodycentric, F, U>::from_vec3(body_params, relative_vec)
    }
}

// =============================================================================
// Bodycentric → Other Centers
// =============================================================================

impl<F: MutableFrame, U: LengthUnit> Position<Bodycentric, F, U>
where
    Quantity<U>: From<AstronomicalUnits>,
    Position<Geocentric, F, U>: TransformFrame<Ecliptic<U, Geocentric>>,
    Ecliptic<U, Geocentric>: TransformFrame<Position<Geocentric, F, U>>,
{
    /// Transform from body-centric to geocentric coordinates.
    ///
    /// This adds back the body's position to recover the geocentric position.
    pub fn to_geocentric(&self, jd: JulianDate) -> Position<Geocentric, F, U> {
        let body_params = *self.center_params();

        // Get the body's position in ecliptic coordinates at the given time
        let body_ecliptic_au = body_params.orbit.kepler_position(jd);

        // Convert the body's position to geocentric
        let body_geo: Position<Geocentric, F, U> = match body_params.orbit_center {
            OrbitReferenceCenter::Geocentric => {
                let body_ecl: Ecliptic<U, Geocentric> = Ecliptic::new(
                    body_ecliptic_au.x(),
                    body_ecliptic_au.y(),
                    body_ecliptic_au.z(),
                );
                body_ecl.to_frame()
            }
            OrbitReferenceCenter::Heliocentric => {
                let body_helio_ecl: Ecliptic<U, Heliocentric> = Ecliptic::new(
                    body_ecliptic_au.x(),
                    body_ecliptic_au.y(),
                    body_ecliptic_au.z(),
                );
                let body_geo_ecl: Ecliptic<U, Geocentric> = body_helio_ecl.to_center(jd);
                body_geo_ecl.to_frame()
            }
            OrbitReferenceCenter::Barycentric => {
                let body_bary_ecl: Ecliptic<U, Barycentric> = Ecliptic::new(
                    body_ecliptic_au.x(),
                    body_ecliptic_au.y(),
                    body_ecliptic_au.z(),
                );
                let body_geo_ecl: Ecliptic<U, Geocentric> = body_bary_ecl.to_center(jd);
                body_geo_ecl.to_frame()
            }
        };

        // Add back the body's position: relative + body = geocentric
        Position::<Geocentric, F, U>::from_vec3_origin(self.as_vec3() + body_geo.as_vec3())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::astro::orbit::Orbit;
    use crate::bodies::solar_system::Earth;
    use crate::coordinates::frames;
    use crate::macros::assert_cartesian_eq;
    use qtty::{AstronomicalUnit, Au, Degrees, AU};

    #[test]
    fn test_geocentric_to_bodycentric_geocentric_orbit() {
        // Create a simple geocentric orbit (satellite at 0.0001 AU from Earth)
        let satellite_orbit = Orbit::new(
            0.0001 * AU,
            0.0,
            Degrees::new(0.0),
            Degrees::new(0.0),
            Degrees::new(0.0),
            Degrees::new(0.0),
            JulianDate::J2000,
        );

        let sat_params = BodycentricParams::geocentric(satellite_orbit);

        // Target at 0.001 AU from Earth
        let target: Position<Geocentric, frames::Ecliptic, AstronomicalUnit> =
            Position::new(0.001, 0.0, 0.0);

        let result = target.to_bodycentric(sat_params, JulianDate::J2000);

        // The satellite is at 0.0001 AU, target at 0.001 AU
        // Relative position should be ~0.0009 AU
        assert!(result.x().value() > 0.0);
        assert!(result.x().value() < 0.001);
    }

    #[test]
    fn test_heliocentric_to_bodycentric_mars_view() {
        // Mars-like orbit
        let mars_orbit = Orbit::new(
            1.524 * AU,
            0.0934,
            Degrees::new(1.85),
            Degrees::new(49.56),
            Degrees::new(286.5),
            Degrees::new(19.41),
            JulianDate::J2000,
        );

        let mars_params = BodycentricParams::heliocentric(mars_orbit);

        // Earth's position (heliocentric)
        let earth_helio = *Earth::vsop87a(JulianDate::J2000).get_position();

        let earth_from_mars: Position<Bodycentric, frames::Ecliptic, Au> =
            earth_helio.to_bodycentric(mars_params, JulianDate::J2000);

        // Mars is further from the Sun than Earth, so the relative position should exist
        assert!(!earth_from_mars.x().value().is_nan());
        assert!(!earth_from_mars.y().value().is_nan());
        assert!(!earth_from_mars.z().value().is_nan());
    }

    #[test]
    fn test_bodycentric_roundtrip() {
        // Create a geocentric orbit
        let satellite_orbit = Orbit::new(
            0.0001 * AU,
            0.0,
            Degrees::new(0.0),
            Degrees::new(0.0),
            Degrees::new(0.0),
            Degrees::new(0.0),
            JulianDate::J2000,
        );

        let sat_params = BodycentricParams::geocentric(satellite_orbit);

        // Original geocentric position
        let original: Position<Geocentric, frames::Ecliptic, AstronomicalUnit> =
            Position::new(0.001, 0.002, 0.003);

        // Transform to bodycentric and back
        let bodycentric = original.to_bodycentric(sat_params, JulianDate::J2000);
        let recovered = bodycentric.to_geocentric(JulianDate::J2000);

        // Should recover the original position
        assert_cartesian_eq!(
            &original,
            &recovered,
            1e-12,
            "Roundtrip should preserve position"
        );
    }

    #[test]
    fn test_body_at_origin() {
        // If we ask for the body's own position in body-centric coords, it should be at origin
        let orbit = Orbit::new(
            0.0001 * AU,
            0.0,
            Degrees::new(0.0),
            Degrees::new(0.0),
            Degrees::new(0.0),
            Degrees::new(0.0),
            JulianDate::J2000,
        );

        let params = BodycentricParams::geocentric(orbit);

        // Get the body's geocentric position
        let body_geo_ecl = orbit.kepler_position(JulianDate::J2000);
        let body_geo: Position<Geocentric, frames::Ecliptic, AstronomicalUnit> = Position::new(
            body_geo_ecl.x().value(),
            body_geo_ecl.y().value(),
            body_geo_ecl.z().value(),
        );

        // Transform to body-centric
        let body_from_body = body_geo.to_bodycentric(params, JulianDate::J2000);

        // Should be at origin
        assert!(
            body_from_body.distance().value().abs() < 1e-10,
            "Body's own position in body-centric should be at origin, got distance {}",
            body_from_body.distance().value()
        );
    }
}
