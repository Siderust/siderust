use crate::astro::aberration::remove_aberration;
use crate::bodies::solar_system::{Earth, Sun};
use crate::coordinates::{
    centers::{Barycentric, Geocentric, Heliocentric},
    frames::MutableFrame,
    cartesian::position::{Position, Ecliptic, Equatorial},
};
use crate::units::{Quantity, AstronomicalUnits, LengthUnit};
use crate::astro::JulianDate;
use crate::coordinates::transform::centers::TransformCenter;
use crate::coordinates::transform::TransformFrame;

impl<F: MutableFrame, U: LengthUnit> TransformCenter<Position<Heliocentric, F, U>>
    for Position<Geocentric, F, U>
where
    Quantity<U>: From<AstronomicalUnits> + PartialEq + std::fmt::Debug,
    Ecliptic<U>: TransformFrame<Equatorial<U, Heliocentric>>, // Required by VSOP
    Position<Geocentric, F, U>: TransformFrame<Equatorial<U, Geocentric>>, // Required by Aberration
    Equatorial<U, Heliocentric>: TransformFrame<Position<Heliocentric, F, U>>, // Required by Aberration
{
    fn to_center(&self, jd: JulianDate) -> Position<Heliocentric, F, U> {
        //geocentric_to_heliocentric(self, jd)
        let earth_helio_ecl_au = *Earth::vsop87a(jd).get_position();
        let earth_helio_ecl = Ecliptic::<U, Heliocentric>::new(
            earth_helio_ecl_au.x(),
            earth_helio_ecl_au.y(),
            earth_helio_ecl_au.z(),
        );

        let earth_helio_equ: Equatorial::<U, Heliocentric> = earth_helio_ecl.to_frame(); // (Helio-Ecl) -> (Helio-Equ)
        let target_geo_equ:  Equatorial::<U, Geocentric> = self.to_frame(); // (Geo-F) -> (Geo-Equ)
        let target_geo_equ_no_aberration = remove_aberration(target_geo_equ, jd);

        let helio_equ = Equatorial::<U, Heliocentric>::from_vec3(
            target_geo_equ_no_aberration.as_vec3() + earth_helio_equ.as_vec3(),
        ); // Geocentric -> Heliocentric
        helio_equ.to_frame() // Equatorial -> F
    }
}

impl<F: MutableFrame, U: LengthUnit> TransformCenter<Position<Heliocentric, F, U>>
    for Position<Barycentric, F, U>
where
    Quantity<U>: From<AstronomicalUnits> + PartialEq + std::fmt::Debug,
    Ecliptic<U, Barycentric>: TransformFrame<Position<Barycentric, F, U>>,
{
    fn to_center(&self, jd: JulianDate) -> Position<Heliocentric, F, U> {
        // Barycentric to Heliocentric
        let sun_bary_ecl_au = *Sun::vsop87e(jd).get_position();
        let sun_bary_ecl = Ecliptic::<U, Barycentric>::new(
            sun_bary_ecl_au.x(),
            sun_bary_ecl_au.y(),
            sun_bary_ecl_au.z(),
        );

        let sun: Position::<Barycentric, F, U> = sun_bary_ecl.to_frame();
        Position::from_vec3(self.as_vec3() - sun.as_vec3())
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::bodies::solar_system::Sun;
    use crate::coordinates::centers::*;
    use crate::macros::assert_cartesian_eq;
    use crate::units::Au;
    use crate::coordinates::transform::Transform;

    const EPSILON: f64 = 1e-9; // Precision tolerance for floating-point comparisons

    #[test] // Barycentric -> Heliocentric
    fn test_bary() {
        let sun_bary = Sun::vsop87e(JulianDate::J2000).get_position().clone();
        let sun_helio: Ecliptic::<Au, Heliocentric> = (&sun_bary).transform(JulianDate::J2000);
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
        let sun = Sun::vsop87e(JulianDate::J2000).get_position().clone();
        let sun_geo: Ecliptic::<Au, Geocentric> = sun.transform(JulianDate::J2000);
        let sun_helio: Ecliptic::<Au> = sun_geo.transform(JulianDate::J2000);
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
