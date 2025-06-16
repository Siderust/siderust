use crate::units::JulianDay;
use crate::bodies::solar_system::{Sun, Earth};
use crate::coordinates::{
    frames::{ReferenceFrame, Ecliptic, Equatorial},
    centers::{Barycentric, Heliocentric, Geocentric},
    cartesian::Position,
};
use crate::coordinates::transform::Transform;
use crate::astro::aberration::remove_aberration;

pub fn heliocentric_to_barycentric<F: ReferenceFrame>(
    helio: &Position<Heliocentric, F>,
    jd: JulianDay
) -> Position<Barycentric, F>
where
    for<'a> Position<Barycentric, F>: From<&'a Position<Barycentric, Ecliptic>>,
{
    let sun = Position::<Barycentric, F>::from(
        Sun::vsop87e(jd).get_position()
    );
    Position::new(
        helio.x() + sun.x(),
        helio.y() + sun.y(),
        helio.z() + sun.z(),
    )
}

pub fn geocentric_to_barycentric<F: ReferenceFrame>(
    geo: &Position<Geocentric, F>,
    jd: JulianDay
) -> Position<Barycentric, F>
where
    for<'a> Position<Barycentric, Equatorial>: From<&'a Position<Barycentric, Ecliptic>>, // Required by VSOP
    for<'a> Position<Geocentric, Equatorial>: From<&'a Position<Geocentric, F>>,   // Required by Aberration
    for<'a> Position<Barycentric, F>: From<&'a Position<Barycentric, Equatorial>>, // Required by Aberration
{
    let earth_bary_ecl = Earth::vsop87e(jd).get_position().clone();
    let earth_bary_equ = Position::<Barycentric, Equatorial>::from(&earth_bary_ecl); // (Bary-Ecl) -> (Bary-Equ)

    let target_geo_equ  = Position::<Geocentric, Equatorial>::from(geo); // (Geo-F) -> (Geo-Equ)
    let target_geo_equ_no_aberration = remove_aberration(target_geo_equ, jd);

    let bary_equ = Position::<Barycentric, Equatorial>::from_vec3(target_geo_equ_no_aberration.as_vec3() + earth_bary_equ.as_vec3()); // Geocentric -> Barycentric
    Position::<Barycentric, F>::from(&bary_equ) // Equatorial -> F
}


impl<F: ReferenceFrame> Transform<Position<Barycentric, F>> for Position<Heliocentric, F>
where
    for<'a> Position<Barycentric, F>: From<&'a Position<Barycentric, Ecliptic>>,
{
    fn transform(&self, jd: JulianDay) -> Position<Barycentric, F> {
        heliocentric_to_barycentric(self, jd)
    }
}

impl<F: ReferenceFrame>  Transform<Position<Barycentric, F>> for Position<Geocentric, F>
where
    for<'a> Position<Barycentric, Equatorial>: From<&'a Position<Barycentric, Ecliptic>>, // Required by VSOP
    for<'a> Position<Geocentric, Equatorial>: From<&'a Position<Geocentric, F>>, // Required by Aberration
    for<'a> Position<Barycentric, F>: From<&'a Position<Barycentric, Equatorial>>, // Required by Aberration
{
    fn transform( &self, jd: JulianDay) -> Position<Barycentric, F> {
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
