// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Bodycentric Position Transformations
//!
//! This module provides [`TransformCenter`] impls for bodycentric coordinates.
//! A body-centric coordinate system has its origin at an orbiting celestial body
//! (satellite, planet, moon, comet, asteroid, etc.).
//!
//! ## Usage
//!
//! Forward (any standard center → Bodycentric):
//! ```rust
//! use siderust::coordinates::transform::TransformCenter;
//! use siderust::coordinates::centers::{Bodycentric, BodycentricParams, Geocentric};
//! use siderust::coordinates::cartesian::Position;
//! use siderust::coordinates::frames;
//! use siderust::astro::orbit::Orbit;
//! use siderust::time::JulianDate;
//! use qtty::*;
//!
//! let satellite_orbit = Orbit::new(
//!     0.0000426 * AU, 0.001,
//!     Degrees::new(51.6), Degrees::new(0.0), Degrees::new(0.0), Degrees::new(0.0),
//!     JulianDate::J2000,
//! );
//! let sat_params = BodycentricParams::geocentric(satellite_orbit);
//!
//! let target_geo: Position<Geocentric, frames::EquatorialMeanJ2000, AstronomicalUnit> =
//!     Position::new(0.00257, 0.0, 0.0);
//!
//! let target_from_sat: Position<Bodycentric, frames::EquatorialMeanJ2000, AstronomicalUnit> = target_geo.to_center((sat_params, JulianDate::J2000));
//! ```
//!
//! Reverse (Bodycentric → Geocentric):
//! ```rust,ignore
//! let geo: Position<Geocentric, _, _> = bodycentric_pos.to_center(jd);
//! ```

use crate::coordinates::cartesian::position::{EclipticMeanJ2000, Position};
use crate::coordinates::centers::{
    Barycentric, Bodycentric, BodycentricParams, Geocentric, Heliocentric, OrbitReferenceCenter,
};
use crate::coordinates::frames::MutableFrame;
use crate::coordinates::transform::centers::TransformCenter;
use crate::coordinates::transform::context::AstroContext;
use crate::coordinates::transform::TransformFrame;
use crate::time::JulianDate;
use qtty::{AstronomicalUnits, LengthUnit, Quantity};

// =============================================================================
// Geocentric → Bodycentric
// =============================================================================

impl<F: MutableFrame, U: LengthUnit> TransformCenter<Bodycentric, F, U>
    for Position<Geocentric, F, U>
where
    Quantity<U>: From<AstronomicalUnits>,
    Position<Geocentric, F, U>: TransformFrame<EclipticMeanJ2000<U, Geocentric>>,
    EclipticMeanJ2000<U, Geocentric>: TransformFrame<Position<Geocentric, F, U>>,
{
    fn to_center_with(
        &self,
        body_params: BodycentricParams,
        jd: JulianDate,
        _ctx: &AstroContext,
    ) -> Position<Bodycentric, F, U> {
        let body_ecliptic_au = body_params.orbit.kepler_position(jd);

        let body_in_our_center: Position<Geocentric, F, U> = match body_params.orbit_center {
            OrbitReferenceCenter::Geocentric => {
                let body_ecl: EclipticMeanJ2000<U, Geocentric> = EclipticMeanJ2000::new(
                    body_ecliptic_au.x(),
                    body_ecliptic_au.y(),
                    body_ecliptic_au.z(),
                );
                body_ecl.to_frame()
            }
            OrbitReferenceCenter::Heliocentric => {
                let body_helio_ecl: EclipticMeanJ2000<U, Heliocentric> = EclipticMeanJ2000::new(
                    body_ecliptic_au.x(),
                    body_ecliptic_au.y(),
                    body_ecliptic_au.z(),
                );
                let body_geo_ecl: EclipticMeanJ2000<U, Geocentric> = body_helio_ecl.to_center(jd);
                body_geo_ecl.to_frame()
            }
            OrbitReferenceCenter::Barycentric => {
                let body_bary_ecl: EclipticMeanJ2000<U, Barycentric> = EclipticMeanJ2000::new(
                    body_ecliptic_au.x(),
                    body_ecliptic_au.y(),
                    body_ecliptic_au.z(),
                );
                let body_geo_ecl: EclipticMeanJ2000<U, Geocentric> = body_bary_ecl.to_center(jd);
                body_geo_ecl.to_frame()
            }
        };

        let relative_vec = self.as_vec3() - body_in_our_center.as_vec3();
        Position::<Bodycentric, F, U>::from_vec3(body_params, relative_vec)
    }
}

// =============================================================================
// Heliocentric → Bodycentric
// =============================================================================

impl<F: MutableFrame, U: LengthUnit> TransformCenter<Bodycentric, F, U>
    for Position<Heliocentric, F, U>
where
    Quantity<U>: From<AstronomicalUnits>,
    Position<Heliocentric, F, U>: TransformFrame<EclipticMeanJ2000<U, Heliocentric>>,
    EclipticMeanJ2000<U, Heliocentric>: TransformFrame<Position<Heliocentric, F, U>>,
{
    fn to_center_with(
        &self,
        body_params: BodycentricParams,
        jd: JulianDate,
        _ctx: &AstroContext,
    ) -> Position<Bodycentric, F, U> {
        let body_ecliptic_au = body_params.orbit.kepler_position(jd);

        let body_in_our_center: Position<Heliocentric, F, U> = match body_params.orbit_center {
            OrbitReferenceCenter::Heliocentric => {
                let body_ecl: EclipticMeanJ2000<U, Heliocentric> = EclipticMeanJ2000::new(
                    body_ecliptic_au.x(),
                    body_ecliptic_au.y(),
                    body_ecliptic_au.z(),
                );
                body_ecl.to_frame()
            }
            OrbitReferenceCenter::Geocentric => {
                let body_geo_ecl: EclipticMeanJ2000<U, Geocentric> = EclipticMeanJ2000::new(
                    body_ecliptic_au.x(),
                    body_ecliptic_au.y(),
                    body_ecliptic_au.z(),
                );
                let body_helio_ecl: EclipticMeanJ2000<U, Heliocentric> = body_geo_ecl.to_center(jd);
                body_helio_ecl.to_frame()
            }
            OrbitReferenceCenter::Barycentric => {
                let body_bary_ecl: EclipticMeanJ2000<U, Barycentric> = EclipticMeanJ2000::new(
                    body_ecliptic_au.x(),
                    body_ecliptic_au.y(),
                    body_ecliptic_au.z(),
                );
                let body_helio_ecl: EclipticMeanJ2000<U, Heliocentric> =
                    body_bary_ecl.to_center(jd);
                body_helio_ecl.to_frame()
            }
        };

        let relative_vec = self.as_vec3() - body_in_our_center.as_vec3();
        Position::<Bodycentric, F, U>::from_vec3(body_params, relative_vec)
    }
}

// =============================================================================
// Barycentric → Bodycentric
// =============================================================================

impl<F: MutableFrame, U: LengthUnit> TransformCenter<Bodycentric, F, U>
    for Position<Barycentric, F, U>
where
    Quantity<U>: From<AstronomicalUnits>,
    Position<Barycentric, F, U>: TransformFrame<EclipticMeanJ2000<U, Barycentric>>,
    EclipticMeanJ2000<U, Barycentric>: TransformFrame<Position<Barycentric, F, U>>,
{
    fn to_center_with(
        &self,
        body_params: BodycentricParams,
        jd: JulianDate,
        _ctx: &AstroContext,
    ) -> Position<Bodycentric, F, U> {
        let body_ecliptic_au = body_params.orbit.kepler_position(jd);

        let body_in_our_center: Position<Barycentric, F, U> = match body_params.orbit_center {
            OrbitReferenceCenter::Barycentric => {
                let body_ecl: EclipticMeanJ2000<U, Barycentric> = EclipticMeanJ2000::new(
                    body_ecliptic_au.x(),
                    body_ecliptic_au.y(),
                    body_ecliptic_au.z(),
                );
                body_ecl.to_frame()
            }
            OrbitReferenceCenter::Heliocentric => {
                let body_helio_ecl: EclipticMeanJ2000<U, Heliocentric> = EclipticMeanJ2000::new(
                    body_ecliptic_au.x(),
                    body_ecliptic_au.y(),
                    body_ecliptic_au.z(),
                );
                let body_bary_ecl: EclipticMeanJ2000<U, Barycentric> = body_helio_ecl.to_center(jd);
                body_bary_ecl.to_frame()
            }
            OrbitReferenceCenter::Geocentric => {
                let body_geo_ecl: EclipticMeanJ2000<U, Geocentric> = EclipticMeanJ2000::new(
                    body_ecliptic_au.x(),
                    body_ecliptic_au.y(),
                    body_ecliptic_au.z(),
                );
                let body_bary_ecl: EclipticMeanJ2000<U, Barycentric> = body_geo_ecl.to_center(jd);
                body_bary_ecl.to_frame()
            }
        };

        let relative_vec = self.as_vec3() - body_in_our_center.as_vec3();
        Position::<Bodycentric, F, U>::from_vec3(body_params, relative_vec)
    }
}

// =============================================================================
// Bodycentric → Geocentric
// =============================================================================

impl<F: MutableFrame, U: LengthUnit> TransformCenter<Geocentric, F, U>
    for Position<Bodycentric, F, U>
where
    Quantity<U>: From<AstronomicalUnits>,
    Position<Geocentric, F, U>: TransformFrame<EclipticMeanJ2000<U, Geocentric>>,
    EclipticMeanJ2000<U, Geocentric>: TransformFrame<Position<Geocentric, F, U>>,
{
    fn to_center_with(
        &self,
        _params: (),
        jd: JulianDate,
        _ctx: &AstroContext,
    ) -> Position<Geocentric, F, U> {
        let body_params = *self.center_params();
        let body_ecliptic_au = body_params.orbit.kepler_position(jd);

        let body_geo: Position<Geocentric, F, U> = match body_params.orbit_center {
            OrbitReferenceCenter::Geocentric => {
                let body_ecl: EclipticMeanJ2000<U, Geocentric> = EclipticMeanJ2000::new(
                    body_ecliptic_au.x(),
                    body_ecliptic_au.y(),
                    body_ecliptic_au.z(),
                );
                body_ecl.to_frame()
            }
            OrbitReferenceCenter::Heliocentric => {
                let body_helio_ecl: EclipticMeanJ2000<U, Heliocentric> = EclipticMeanJ2000::new(
                    body_ecliptic_au.x(),
                    body_ecliptic_au.y(),
                    body_ecliptic_au.z(),
                );
                let body_geo_ecl: EclipticMeanJ2000<U, Geocentric> = body_helio_ecl.to_center(jd);
                body_geo_ecl.to_frame()
            }
            OrbitReferenceCenter::Barycentric => {
                let body_bary_ecl: EclipticMeanJ2000<U, Barycentric> = EclipticMeanJ2000::new(
                    body_ecliptic_au.x(),
                    body_ecliptic_au.y(),
                    body_ecliptic_au.z(),
                );
                let body_geo_ecl: EclipticMeanJ2000<U, Geocentric> = body_bary_ecl.to_center(jd);
                body_geo_ecl.to_frame()
            }
        };

        Position::<Geocentric, F, U>::from_vec3_origin(self.as_vec3() + body_geo.as_vec3())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::astro::orbit::Orbit;
    use crate::bodies::solar_system::Earth;
    use crate::coordinates::frames;
    use crate::coordinates::transform::TransformCenter;
    use crate::macros::assert_cartesian_eq;
    use qtty::{AstronomicalUnit, Au, Degrees, AU};

    #[test]
    fn test_geocentric_to_bodycentric_geocentric_orbit() {
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

        let target: Position<Geocentric, frames::EclipticMeanJ2000, AstronomicalUnit> =
            Position::new(0.001, 0.0, 0.0);

        let result: Position<Bodycentric, frames::EclipticMeanJ2000, AstronomicalUnit> =
            target.to_center((sat_params, JulianDate::J2000));

        assert!(result.x() > 0.0);
        assert!(result.x() < 0.001);
    }

    #[test]
    fn test_heliocentric_to_bodycentric_mars_view() {
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

        let earth_helio = Earth::vsop87a(JulianDate::J2000);

        let earth_from_mars: Position<Bodycentric, frames::EclipticMeanJ2000, Au> =
            earth_helio.to_center((mars_params, JulianDate::J2000));

        assert!(earth_from_mars.x().is_finite());
        assert!(earth_from_mars.y().is_finite());
        assert!(earth_from_mars.z().is_finite());
    }

    #[test]
    fn test_bodycentric_roundtrip() {
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

        let original: Position<Geocentric, frames::EclipticMeanJ2000, AstronomicalUnit> =
            Position::new(0.001, 0.002, 0.003);

        let bodycentric: Position<Bodycentric, frames::EclipticMeanJ2000, AstronomicalUnit> =
            original.to_center((sat_params, JulianDate::J2000));
        let recovered: Position<Geocentric, frames::EclipticMeanJ2000, AstronomicalUnit> =
            bodycentric.to_center(JulianDate::J2000);

        assert_cartesian_eq!(
            &original,
            &recovered,
            1e-12,
            "Roundtrip should preserve position"
        );
    }

    #[test]
    fn test_body_at_origin() {
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

        let body_geo_ecl = orbit.kepler_position(JulianDate::J2000);
        let body_geo: Position<Geocentric, frames::EclipticMeanJ2000, AstronomicalUnit> =
            Position::from_vec3_origin(*body_geo_ecl.as_vec3());

        let body_from_body: Position<Bodycentric, frames::EclipticMeanJ2000, AstronomicalUnit> =
            body_geo.to_center((params, JulianDate::J2000));

        assert!(
            body_from_body.distance().abs() < 1e-10,
            "Body's own position in body-centric should be at origin, got distance {}",
            body_from_body.distance()
        );
    }
}
