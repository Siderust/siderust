use crate::units::{JulianDay, Unit};
use crate::bodies::solar_system::{Sun, Earth};
use crate::coordinates::{
    frames::{MutableFrame, Ecliptic, Equatorial},
    centers::{Barycentric, Heliocentric, Geocentric},
    cartesian::Position,
};
use crate::coordinates::transform::Transform;
use crate::astro::aberration::remove_aberration;

pub fn heliocentric_to_barycentric<F: MutableFrame, U: Unit>(
    helio_f: &Position<Heliocentric, F, U>,
    jd: JulianDay
) -> Position<Barycentric, F, U>
where
    for<'a> Position<Barycentric, F, U>: From<&'a Position<Barycentric, Ecliptic, U>>,
{
    let sun_bary_ecl = Sun::vsop87e(jd).get_position();
    let sun_bary_f = Position::<Barycentric, F>::from(sun_bary_ecl);
    Position::from_vec3(helio_f.as_vec3() + sun_bary_f.as_vec3())
}

pub fn geocentric_to_barycentric<F: MutableFrame, U: Unit>(
    geo: &Position<Geocentric, F, U>,
    jd: JulianDay
) -> Position<Barycentric, F, U>
where
    for<'a> Position<Barycentric, Equatorial, U>: From<&'a Position<Barycentric, Ecliptic, U>>, // Required by VSOP
    for<'a> Position<Geocentric, Equatorial, U>: From<&'a Position<Geocentric, F, U>>,   // Required by Aberration
    for<'a> Position<Barycentric, F, U>: From<&'a Position<Barycentric, Equatorial, U>>, // Required by Aberration
{
    let earth_bary_ecl = Earth::vsop87e(jd).get_position().clone();
    let earth_bary_equ = Position::<Barycentric, Equatorial>::from(&earth_bary_ecl); // (Bary-Ecl) -> (Bary-Equ)

    let target_geo_equ  = Position::<Geocentric, Equatorial>::from(geo); // (Geo-F) -> (Geo-Equ)
    let target_geo_equ_no_aberration = remove_aberration(target_geo_equ, jd);

    let bary_equ = Position::<Barycentric, Equatorial, U>::from_vec3(target_geo_equ_no_aberration.as_vec3() + earth_bary_equ.as_vec3()); // Geocentric -> Barycentric
    Position::<Barycentric, F, U>::from(&bary_equ) // Equatorial -> F
}

impl<F: MutableFrame, U: Unit> Transform<Position<Barycentric, F, U>> for Position<Heliocentric, F, U>
where
    for<'a> Position<Barycentric, F, U>: From<&'a Position<Barycentric, Ecliptic, U>>,
{
    fn transform(&self, jd: JulianDay) -> Position<Barycentric, F, U> {
        heliocentric_to_barycentric(self, jd)
    }
}

impl<F: MutableFrame, U: Unit> Transform<Position<Barycentric, F, U>> for Position<Geocentric, F, U>
where
    for<'a> Position<Barycentric, Equatorial, U>: From<&'a Position<Barycentric, Ecliptic, U>>, // Required by VSOP
    for<'a> Position<Geocentric, Equatorial, U>: From<&'a Position<Geocentric, F, U>>, // Required by Aberration
    for<'a> Position<Barycentric, F, U>: From<&'a Position<Barycentric, Equatorial, U>>, // Required by Aberration
{
    fn transform(&self, jd: JulianDay) -> Position<Barycentric, F, U> {
        geocentric_to_barycentric(self, jd)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coordinates::frames::*;
    use crate::coordinates::centers::*;
    use crate::bodies::solar_system::{Sun, Earth};
    use crate::macros::assert_cartesian_eq;

    const EPSILON: f64 = 1e-9; // Precision tolerance for floating-point comparisons

    #[test] // Heliocentric -> Barycentric
    fn test_helio() {
        let sun_helio = Position::<Heliocentric, Ecliptic>::new(0.0, 0.0, 0.0);
        let sun_bary = Position::<Barycentric, Ecliptic>::from(&sun_helio);
        let expected_sun_bary = Sun::vsop87e(JulianDay::J2000).get_position().clone();
        assert_cartesian_eq!(&sun_bary, &expected_sun_bary, EPSILON);
    }

    #[test] // Geocentric -> Barycentric
    fn test_geo() {
        let earth_geo = Position::<Geocentric, Ecliptic>::new(0.0, 0.0, 0.0);
        let earth_bary = Position::<Barycentric, Ecliptic>::from(&earth_geo);
        let expected_earth_bary = Earth::vsop87e(JulianDay::J2000).get_position().clone();
        assert_cartesian_eq!(&earth_bary, &expected_earth_bary, EPSILON);
    }

}
