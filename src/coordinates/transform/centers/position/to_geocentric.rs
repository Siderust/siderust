// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

use crate::coordinates::transform::centers::TransformCenter;
use crate::coordinates::transform::context::AstroContext;
use crate::coordinates::transform::providers::CenterShiftProvider;
use crate::coordinates::{
    cartesian::Position,
    centers::{Barycentric, Geocentric, Heliocentric},
    frames::{self, MutableFrame},
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
//
// This implementation delegates to CenterShiftProvider, which sources positions
// from the ephemeris trait. Using default AstroContext provides VSOP87 positions.

impl<F: MutableFrame, U: LengthUnit> TransformCenter<Position<Geocentric, F, U>>
    for Position<Barycentric, F, U>
where
    Quantity<U>: From<AstronomicalUnits>,
    (): crate::coordinates::transform::FrameRotationProvider<frames::Ecliptic, F>,
    (): CenterShiftProvider<Barycentric, Geocentric, F>,
{
    fn to_center(&self, jd: JulianDate) -> Position<Geocentric, F, U> {
        let ctx = AstroContext::default();
        let [sx, sy, sz] = <() as CenterShiftProvider<Barycentric, Geocentric, F>>::shift(jd, &ctx);

        // The shift is in AU; convert to unit U
        let x = self.x().value() + sx;
        let y = self.y().value() + sy;
        let z = self.z().value() + sz;

        Position::<Geocentric, F, U>::new(
            Quantity::<U>::new(x),
            Quantity::<U>::new(y),
            Quantity::<U>::new(z),
        )
    }
}

// =============================================================================
// Heliocentric → Geocentric (pure translation, no aberration)
// =============================================================================

impl<F: MutableFrame, U: LengthUnit> TransformCenter<Position<Geocentric, F, U>>
    for Position<Heliocentric, F, U>
where
    Quantity<U>: From<AstronomicalUnits>,
    (): crate::coordinates::transform::FrameRotationProvider<frames::Ecliptic, F>,
    (): CenterShiftProvider<Heliocentric, Geocentric, F>,
{
    fn to_center(&self, jd: JulianDate) -> Position<Geocentric, F, U> {
        let ctx = AstroContext::default();
        let [sx, sy, sz] =
            <() as CenterShiftProvider<Heliocentric, Geocentric, F>>::shift(jd, &ctx);

        // Apply the shift
        let x = self.x().value() + sx;
        let y = self.y().value() + sy;
        let z = self.z().value() + sz;

        Position::<Geocentric, F, U>::new(
            Quantity::<U>::new(x),
            Quantity::<U>::new(y),
            Quantity::<U>::new(z),
        )
    }
}

#[cfg(test)]
mod tests {
    use crate::bodies::solar_system::Earth;
    use crate::coordinates::{cartesian, centers::*, spherical, transform::Transform};
    use crate::macros::assert_cartesian_eq;
    use crate::time::JulianDate;
    use qtty::*;

    const EPSILON: f64 = 1e-9; // Precision tolerance for floating-point comparisons

    #[test] // Barycentric -> Geocentric
    fn test_bary_to_geo() {
        let earth_bary = *Earth::vsop87e(JulianDate::J2000).get_position();
        let earth_geo: cartesian::position::Ecliptic<Au, Geocentric> =
            earth_bary.transform(JulianDate::J2000);
        let expected_earth_geo = cartesian::position::Ecliptic::<Au, Geocentric>::CENTER;
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
        let earth_helio = *Earth::vsop87a(JulianDate::J2000).get_position();
        let earth_geo: cartesian::position::Ecliptic<Au, Geocentric> =
            earth_helio.transform(JulianDate::J2000);
        let expected_earth_geo = cartesian::position::Ecliptic::<Au, Geocentric>::CENTER;
        // Note: The tolerance here is slightly relaxed because heliocentric transforms
        // now go through barycentric as a hub (vsop87e), while vsop87a is the original
        // heliocentric series. The difference (~5e-9 AU ≈ 0.7m) is well within VSOP87 accuracy.
        assert_cartesian_eq!(&earth_geo, &expected_earth_geo, 1e-8);
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
