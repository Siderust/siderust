use crate::astro::aberration::remove_aberration;
use crate::astro::JulianDate;
use crate::bodies::solar_system::{Earth, Sun};
use crate::coordinates::transform::centers::TransformCenter;
use crate::coordinates::transform::TransformFrame;
use crate::coordinates::{
    cartesian::position::{Ecliptic, Equatorial},
    cartesian::Position,
    centers::{Barycentric, Geocentric, Heliocentric},
    frames::MutableFrame,
};
use qtty::{AstronomicalUnits, LengthUnit, Quantity};

// Heliocentric To Barycentric
impl<F: MutableFrame, U: LengthUnit> TransformCenter<Position<Barycentric, F, U>>
    for Position<Heliocentric, F, U>
where
    Quantity<U>: From<AstronomicalUnits>,
    Ecliptic<U, Barycentric>: TransformFrame<Position<Barycentric, F, U>>,
{
    fn to_center(&self, jd: JulianDate) -> Position<Barycentric, F, U> {
        let sun_bary_ecl_au = *Sun::vsop87e(jd).get_position();

        // VSOP87 gives the Sun's position in AstronomicalUnits
        let sun_bary_ecl = Ecliptic::<U, Barycentric>::new(
            sun_bary_ecl_au.x(),
            sun_bary_ecl_au.y(),
            sun_bary_ecl_au.z(),
        );

        let sun_bary_f: Position<Barycentric, F, U> = sun_bary_ecl.to_frame(); // (Bary-Ecl) -> (Bary-F)
        Position::from_vec3(self.as_vec3() + sun_bary_f.as_vec3())
    }
}

// Geocentric To Barycentric
impl<F: MutableFrame, U: LengthUnit> TransformCenter<Position<Barycentric, F, U>>
    for Position<Geocentric, F, U>
where
    Quantity<U>: From<AstronomicalUnits>,
    Ecliptic<U, Barycentric>: TransformFrame<Equatorial<U, Barycentric>>, // Required by VSOP
    Position<Geocentric, F, U>: TransformFrame<Equatorial<U, Geocentric>>, // Required by Aberration
    Equatorial<U, Barycentric>: TransformFrame<Position<Barycentric, F, U>>, // Required by Aberration
{
    fn to_center(&self, jd: JulianDate) -> Position<Barycentric, F, U> {
        let earth_bary_ecl_au = *Earth::vsop87e(jd).get_position();

        // VSOP87 gives the Earth's position in AstronomicalUnits
        let earth_bary_ecl = Ecliptic::<U, Barycentric>::new(
            earth_bary_ecl_au.x(),
            earth_bary_ecl_au.y(),
            earth_bary_ecl_au.z(),
        );

        let earth_bary_equ: Equatorial<U, Barycentric> = earth_bary_ecl.to_frame(); // (Bary-Ecl) -> (Bary-Equ)
        let target_geo_equ: Equatorial<U, Geocentric> = self.to_frame(); // (Geo-F) -> (Geo-Equ)
        let target_geo_equ_no_aberration = remove_aberration(target_geo_equ, jd);

        let bary_equ = Equatorial::<U, Barycentric>::from_vec3(
            target_geo_equ_no_aberration.as_vec3() + earth_bary_equ.as_vec3(),
        ); // Geocentric -> Barycentric
        bary_equ.to_frame() // Equatorial -> F
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bodies::solar_system::{Earth, Sun};
    use crate::coordinates::cartesian::position::Ecliptic;
    use crate::coordinates::centers::*;
    use crate::coordinates::transform::Transform;
    use crate::macros::assert_cartesian_eq;
    use qtty::Au;

    const EPSILON: f64 = 1e-9; // Precision tolerance for floating-point comparisons

    #[test] // Heliocentric -> Barycentric
    fn test_helio() {
        let sun_helio = Ecliptic::<Au>::CENTER;
        let sun_bary: Ecliptic<Au, Barycentric> = sun_helio.transform(JulianDate::J2000);
        let expected_sun_bary = *Sun::vsop87e(JulianDate::J2000).get_position();
        assert_cartesian_eq!(sun_bary, expected_sun_bary, EPSILON);
    }

    #[test] // Geocentric -> Barycentric
    fn test_geo() {
        let earth_geo = Ecliptic::<Au, Geocentric>::CENTER;
        let earth_bary: Ecliptic<Au, Barycentric> = earth_geo.transform(JulianDate::J2000);
        let expected_earth_bary = *Earth::vsop87e(JulianDate::J2000).get_position();
        assert_cartesian_eq!(&earth_bary, &expected_earth_bary, EPSILON);
    }
}
