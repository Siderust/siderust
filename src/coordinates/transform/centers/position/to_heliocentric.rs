use crate::astro::JulianDate;
use crate::bodies::solar_system::{Earth, Sun};
use crate::coordinates::transform::centers::TransformCenter;
use crate::coordinates::{
    cartesian::position::{Ecliptic, Position},
    centers::{Barycentric, Geocentric, Heliocentric},
    frames::{self, MutableFrame},
    transform::Transform,
};
use qtty::{AstronomicalUnits, LengthUnit, Quantity};

// =============================================================================
// Geocentric → Heliocentric (pure translation, no aberration)
// =============================================================================
//
// Center transforms are pure geometry: they only translate positions from one
// origin to another. Aberration is an observation-model effect that depends on
// observer velocity and must be applied explicitly via the observation module.

impl<F: MutableFrame, U: LengthUnit> TransformCenter<Position<Heliocentric, F, U>>
    for Position<Geocentric, F, U>
where
    Quantity<U>: From<AstronomicalUnits> + PartialEq + std::fmt::Debug,
    (): crate::coordinates::transform::FrameRotationProvider<frames::Ecliptic, F>,
{
    fn to_center(&self, jd: JulianDate) -> Position<Heliocentric, F, U> {
        // Get Earth's position in heliocentric ecliptic coordinates
        let earth_helio_ecl_au = *Earth::vsop87a(jd).get_position();
        let earth_helio_ecl = Ecliptic::<U, Heliocentric>::new(
            earth_helio_ecl_au.x(),
            earth_helio_ecl_au.y(),
            earth_helio_ecl_au.z(),
        );

        // Transform Earth to the target frame
        let earth: Position<Heliocentric, F, U> = earth_helio_ecl.transform(jd);

        // Pure translation: heliocentric = geocentric + earth_position
        Position::<Heliocentric, F, U>::from_vec3_origin(self.as_vec3() + earth.as_vec3())
    }
}

impl<F: MutableFrame, U: LengthUnit> TransformCenter<Position<Heliocentric, F, U>>
    for Position<Barycentric, F, U>
where
    Quantity<U>: From<AstronomicalUnits> + PartialEq + std::fmt::Debug,
    (): crate::coordinates::transform::FrameRotationProvider<frames::Ecliptic, F>,
{
    fn to_center(&self, jd: JulianDate) -> Position<Heliocentric, F, U> {
        // Barycentric to Heliocentric
        let sun_bary_ecl_au = *Sun::vsop87e(jd).get_position();
        let sun_bary_ecl = Ecliptic::<U, Barycentric>::new(
            sun_bary_ecl_au.x(),
            sun_bary_ecl_au.y(),
            sun_bary_ecl_au.z(),
        );

        let sun: Position<Barycentric, F, U> = sun_bary_ecl.transform(jd);
        Position::from_vec3_origin(self.as_vec3() - sun.as_vec3())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bodies::solar_system::Sun;
    use crate::coordinates::centers::*;
    use crate::coordinates::transform::Transform;
    use crate::macros::assert_cartesian_eq;
    use qtty::Au;

    const EPSILON: f64 = 1e-9; // Precision tolerance for floating-point comparisons

    #[test] // Barycentric -> Heliocentric
    fn test_bary() {
        let sun_bary = *Sun::vsop87e(JulianDate::J2000).get_position();
        let sun_helio: Ecliptic<Au, Heliocentric> = sun_bary.transform(JulianDate::J2000);
        let expected_sun_helio = Ecliptic::<Au>::CENTER;
        assert_cartesian_eq!(
            &sun_helio,
            &expected_sun_helio,
            EPSILON,
            "Sun in Heliocentric shall be (0,0,0). Current Value {:?}",
            sun_helio
        );
    }

    #[test] // Geocentric -> Heliocentric
    fn test_geo() {
        let sun = *Sun::vsop87e(JulianDate::J2000).get_position();
        let sun_geo: Ecliptic<Au, Geocentric> = sun.transform(JulianDate::J2000);
        let sun_helio: Ecliptic<Au> = sun_geo.transform(JulianDate::J2000);
        let expected_sun_helio = Ecliptic::<Au>::CENTER;
        assert_cartesian_eq!(
            &sun_helio,
            &expected_sun_helio,
            1e-8,
            "Sun in Heliocentric shall be (0,0,0). Current Value {:?}",
            sun_helio
        );
    }
}
