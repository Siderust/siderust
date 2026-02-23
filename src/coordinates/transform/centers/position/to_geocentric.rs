// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

use crate::calculus::ephemeris::Ephemeris;
use crate::coordinates::transform::centers::TransformCenter;
use crate::coordinates::transform::context::DefaultEphemeris;
use crate::coordinates::{
    cartesian::position::{EclipticMeanJ2000, Position},
    centers::{Barycentric, Geocentric, Heliocentric},
    frames::{self, MutableFrame},
    transform::Transform,
};
use crate::time::JulianDate;
use qtty::{AstronomicalUnits, LengthUnit, Quantity};

// =============================================================================
// Barycentric → Geocentric (pure translation, no aberration)
// =============================================================================
//
// Center transforms are pure geometry: they only translate positions from one
// origin to another. Aberration is an observation-model effect that depends on
// observer velocity and must be applied explicitly via the observation module.

impl<F: MutableFrame, U: LengthUnit> TransformCenter<Position<Geocentric, F, U>>
    for Position<Barycentric, F, U>
where
    Quantity<U>: From<AstronomicalUnits>,
    (): crate::coordinates::transform::FrameRotationProvider<frames::EclipticMeanJ2000, F>,
{
    fn to_center(&self, jd: JulianDate) -> Position<Geocentric, F, U> {
        // Get Earth's position in barycentric ecliptic coordinates
        let earth_bary_ecl_au = DefaultEphemeris::earth_barycentric(jd);
        let earth_ecl = EclipticMeanJ2000::<U, Barycentric>::new(
            earth_bary_ecl_au.x(),
            earth_bary_ecl_au.y(),
            earth_bary_ecl_au.z(),
        );

        // Transform Earth to the target frame
        let earth: Position<Barycentric, F, U> = earth_ecl.transform(jd);

        // Pure translation: geocentric = barycentric - earth_position
        Position::<Geocentric, F, U>::from_vec3_origin(self.as_vec3() - earth.as_vec3())
    }
}

// =============================================================================
// Heliocentric → Geocentric (pure translation, no aberration)
// =============================================================================

impl<F: MutableFrame, U: LengthUnit> TransformCenter<Position<Geocentric, F, U>>
    for Position<Heliocentric, F, U>
where
    Quantity<U>: From<AstronomicalUnits>,
    (): crate::coordinates::transform::FrameRotationProvider<frames::EclipticMeanJ2000, F>,
{
    fn to_center(&self, jd: JulianDate) -> Position<Geocentric, F, U> {
        // Get Earth's position in heliocentric ecliptic coordinates
        let earth_helio_ecl_au = DefaultEphemeris::earth_heliocentric(jd);
        let earth_ecl = EclipticMeanJ2000::<U>::new(
            earth_helio_ecl_au.x(),
            earth_helio_ecl_au.y(),
            earth_helio_ecl_au.z(),
        );

        // Transform Earth to the target frame
        let earth: Position<Heliocentric, F, U> = earth_ecl.transform(jd);

        // Pure translation: geocentric = heliocentric - earth_position
        Position::<Geocentric, F, U>::from_vec3_origin(self.as_vec3() - earth.as_vec3())
    }
}

#[cfg(test)]
mod tests {
    use crate::calculus::ephemeris::Ephemeris;
    use crate::coordinates::transform::context::DefaultEphemeris;
    use crate::coordinates::{cartesian, centers::*, spherical, transform::Transform};
    use crate::macros::assert_cartesian_eq;
    use crate::time::JulianDate;
    use qtty::*;

    const EPSILON: f64 = 1e-9; // Precision tolerance for floating-point comparisons

    #[test] // Barycentric -> Geocentric
    fn test_bary_to_geo() {
        let earth_bary = DefaultEphemeris::earth_barycentric(JulianDate::J2000);
        let earth_geo: cartesian::position::EclipticMeanJ2000<Au, Geocentric> =
            earth_bary.transform(JulianDate::J2000);
        let expected_earth_geo = cartesian::position::EclipticMeanJ2000::<Au, Geocentric>::CENTER;
        assert_cartesian_eq!(
            &earth_geo,
            &expected_earth_geo,
            EPSILON,
            "Earth in Geocentric shall be (0,0,0). Current Value {:?}",
            earth_geo
        );
    }

    #[test] // Heliocentric -> Geocentric
    fn test_helio_to_geo() {
        let earth_helio = DefaultEphemeris::earth_heliocentric(JulianDate::J2000);
        let earth_geo: cartesian::position::EclipticMeanJ2000<Au, Geocentric> =
            earth_helio.transform(JulianDate::J2000);
        let expected_earth_geo = cartesian::position::EclipticMeanJ2000::<Au, Geocentric>::CENTER;
        assert_cartesian_eq!(&earth_geo, &expected_earth_geo, EPSILON);
    }

    #[test]
    /// Test that center transforms are pure translations (no aberration).
    ///
    /// This test verifies that ICRS → GCRS is a geometric translation only.
    /// Aberration must be applied separately using the observation module.
    fn test_icrs_to_gcrs_is_pure_translation() {
        const AU_PER_PC: f64 = 206_264.806_f64;
        const SIRIUS_PARALLAX: f64 = 0.37921_f64; // arcsec (Hipparcos van Leeuwen 2007)
        let sirius_distance_au = (1.0 / SIRIUS_PARALLAX) * AU_PER_PC;

        // Sirius catalog position (astrometric)
        let sirius_barycentric_spherical = spherical::position::ICRS::<Au>::new(
            Degrees::new(101.287_155_33), // RA
            Degrees::new(-16.716_115_86), // Dec
            sirius_distance_au,
        );

        let jd = JulianDate::new(2460792.157638889);

        let sirius_barycentric_cartesian =
            cartesian::position::ICRS::<Au>::from_spherical(&sirius_barycentric_spherical);
        let sirius_geocentric_cartesian: cartesian::position::GCRS<Au> =
            Transform::transform(&sirius_barycentric_cartesian, jd);
        let sirius_geocentric_spherical =
            spherical::position::GCRS::<Au>::from_cartesian(&sirius_geocentric_cartesian);

        // Distance should be preserved (approximately - slight change due to Earth's offset)
        assert!(
            (sirius_geocentric_spherical.distance - AstronomicalUnits::new(sirius_distance_au))
                .abs()
                < AstronomicalUnits::new(100.0),
            "Distance should be approximately preserved"
        );

        // Coordinates should be close to the original (small parallax at stellar distances)
        // The shift should be very small for distant stars
        let delta_ra = (sirius_geocentric_spherical.azimuth - Degrees::new(101.287_155_33)).abs();
        let delta_dec = (sirius_geocentric_spherical.polar - Degrees::new(-16.716_115_86)).abs();

        // For stars at ~500,000 AU, Earth's ~1 AU offset causes ~1/500000 radian ≈ 0.4 arcsec
        // change in direction, or about 0.0001 degrees
        assert!(
            delta_ra < Degrees::new(0.001) && delta_dec < Degrees::new(0.001),
            "Astrometric position should change only slightly due to parallax: dRA={}, dDec={}",
            delta_ra,
            delta_dec
        );
    }
}
