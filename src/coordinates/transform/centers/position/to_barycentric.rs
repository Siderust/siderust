use crate::units::{LengthUnit, AstronomicalUnits, Quantity};
use crate::astro::JulianDate;
use crate::bodies::solar_system::{Sun, Earth};
use crate::coordinates::{
    frames::MutableFrame,
    centers::{Barycentric, Heliocentric, Geocentric},
    cartesian::Position,
    cartesian::position::{Ecliptic, Equatorial}
};
use crate::coordinates::transform::Transform;
use crate::coordinates::transform::centers::TransformCenter;
use crate::astro::aberration::remove_aberration;

// Heliocentric To Barycentric
impl<F: MutableFrame, U: LengthUnit> TransformCenter<Position<Barycentric, F, U>> for Position<Heliocentric, F, U>
where
    Quantity<U>: From<AstronomicalUnits>,
    for<'a> Position<Barycentric, F, U>: From<&'a Ecliptic<U, Barycentric>>,
{
    fn to_center(&self, jd: JulianDate) -> Position<Barycentric, F, U> {
        let sun_bary_ecl_au = Sun::vsop87e(jd).get_position().clone();

        // VSOP87 gives the Sun's position in AstronomicalUnits
        let sun_bary_ecl = Ecliptic::<U, Barycentric>::new(
            sun_bary_ecl_au.x(),
            sun_bary_ecl_au.y(),
            sun_bary_ecl_au.z(),
        );

        let sun_bary_f = Position::<Barycentric, F, U>::from(&sun_bary_ecl);
        Position::from_vec3(self.as_vec3() + sun_bary_f.as_vec3())
    }
}

// Geocentric To Barycentric
impl<F: MutableFrame, U: LengthUnit> TransformCenter<Position<Barycentric, F, U>> for Position<Geocentric, F, U>
where
    Quantity<U>: From<AstronomicalUnits>,
    for<'a> Equatorial<U, Barycentric>: From<&'a Ecliptic<U, Barycentric>>, // Required by VSOP
    for<'a> Equatorial<U>: From<&'a Position<Geocentric, F, U>>, // Required by Aberration
    for<'a> Position<Barycentric, F, U>: From<&'a Equatorial<U, Barycentric>>, // Required by Aberration
{
    fn to_center(&self, jd: JulianDate) -> Position<Barycentric, F, U> {
        let earth_bary_ecl_au = Earth::vsop87e(jd).get_position().clone();

        // VSOP87 gives the Earth's position in AstronomicalUnits
        let earth_bary_ecl = Ecliptic::<U, Barycentric>::new(
            earth_bary_ecl_au.x(),
            earth_bary_ecl_au.y(),
            earth_bary_ecl_au.z(),
        );

        let earth_bary_equ = Equatorial::<U, Barycentric>::from(&earth_bary_ecl); // (Bary-Ecl) -> (Bary-Equ)
        let target_geo_equ  = Equatorial::<U, Geocentric>::from(self); // (Geo-F) -> (Geo-Equ)
        let target_geo_equ_no_aberration = remove_aberration(target_geo_equ, jd);

        let bary_equ = Equatorial::<U, Barycentric>::from_vec3(target_geo_equ_no_aberration.as_vec3() + earth_bary_equ.as_vec3()); // Geocentric -> Barycentric
        Position::<Barycentric, F, U>::from(&bary_equ) // Equatorial -> F
    }
}

impl<F: MutableFrame, U: LengthUnit> Transform<Position<Barycentric, F, U>> for Position<Heliocentric, F, U>
where
    Quantity<U>: From<AstronomicalUnits>,
    for<'a> Position<Barycentric, F, U>: From<&'a Ecliptic<U, Barycentric>>,
{
    fn transform(&self, jd: JulianDate) -> Position<Barycentric, F, U> {
        self.to_center(jd)
    }
}

impl<F: MutableFrame, U: LengthUnit> Transform<Position<Barycentric, F, U>> for Position<Geocentric, F, U>
where
    Quantity<U>: From<AstronomicalUnits>,
    for<'a> Equatorial<U, Barycentric>: From<&'a Ecliptic<U, Barycentric>>, // Required by VSOP
    for<'a> Equatorial<U>: From<&'a Position<Geocentric, F, U>>, // Required by Aberration
    for<'a> Position<Barycentric, F, U>: From<&'a Equatorial<U, Barycentric>>, // Required by Aberration
{
    fn transform(&self, jd: JulianDate) -> Position<Barycentric, F, U> {
        self.to_center(jd)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coordinates::centers::*;
    use crate::bodies::solar_system::{Sun, Earth};
    use crate::macros::assert_cartesian_eq;
    use crate::units::Au;
    use crate::coordinates::cartesian::position::Ecliptic;

    const EPSILON: f64 = 1e-9; // Precision tolerance for floating-point comparisons

    #[test] // Heliocentric -> Barycentric
    fn test_helio() {
        let sun_helio = Ecliptic::<Au>::CENTER;
        let sun_bary = Ecliptic::<Au, Barycentric>::from(&sun_helio);
        let expected_sun_bary = Sun::vsop87e(JulianDate::J2000).get_position().clone();
        assert_cartesian_eq!(sun_bary, expected_sun_bary, EPSILON);
    }

    #[test] // Geocentric -> Barycentric
    fn test_geo() {
        let earth_geo = Ecliptic::<Au, Geocentric>::CENTER;
        let earth_bary = Ecliptic::<Au, Barycentric>::from(&earth_geo);
        let expected_earth_bary = Earth::vsop87e(JulianDate::J2000).get_position().clone();
        assert_cartesian_eq!(&earth_bary, &expected_earth_bary, EPSILON);
    }

}
