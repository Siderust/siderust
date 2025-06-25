use crate::units::{JulianDay, Unit};
use crate::bodies::solar_system::{Sun, Earth};
use crate::coordinates::{
    cartesian::Position,
    frames::{MutableFrame, Ecliptic, Equatorial},
    centers::{Heliocentric, Barycentric, Geocentric}
};
use crate::coordinates::transform::Transform;
use crate::astro::aberration::remove_aberration;

pub fn barycentric_to_heliocentric<F: MutableFrame, U: Unit>(
    bary: &Position<Barycentric, F, U>,
    jd: JulianDay
) -> Position<Heliocentric, F, U>
where
    for<'a> Position<Barycentric, F, U>: From<&'a Position<Barycentric, Ecliptic, U>>,
{
    // VSOP87 gives the Sun's position in AstronomicalUnits, so we need to convert to U
    let sun_bary_ecl_au = Sun::vsop87e(jd).get_position().clone();
    let x: U = sun_bary_ecl_au.x().into();
    let y: U = sun_bary_ecl_au.y().into();
    let z: U = sun_bary_ecl_au.z().into();
    let sun_bary_ecl = Position::<Barycentric, Ecliptic, U>::new(x, y, z);

    let sun = Position::<Barycentric, F, U>::from(&sun_bary_ecl);
    Position::from_vec3(bary.as_vec3() - sun.as_vec3())
}

pub fn geocentric_to_heliocentric<F: MutableFrame, U: Unit>(
    geo: &Position<Geocentric, F, U>,
    jd: JulianDay
) -> Position<Heliocentric, F, U>
where
    for<'a> Position<Heliocentric, Equatorial, U>: From<&'a Position<Heliocentric, Ecliptic, U>>, // Required by VSOP
    for<'a> Position<Geocentric, Equatorial, U>: From<&'a Position<Geocentric, F, U>>,   // Required by Aberration
    for<'a> Position<Heliocentric, F, U>: From<&'a Position<Heliocentric, Equatorial, U>>, // Required by Aberration
{
    // VSOP87 gives the Earth's position in AstronomicalUnits, so we need to convert to U
    let earth_helio_ecl_au = Earth::vsop87a(jd).get_position().clone();
    let x: U = earth_helio_ecl_au.x().into();
    let y: U = earth_helio_ecl_au.y().into();
    let z: U = earth_helio_ecl_au.z().into();
    let earth_helio_ecl = Position::<Heliocentric, Ecliptic, U>::new(x, y, z);

    let earth_helio_equ = Position::<Heliocentric, Equatorial, U>::from(&earth_helio_ecl); // (Helio-Ecl) -> (Helio-Equ)

    let target_geo_equ  = Position::<Geocentric, Equatorial, U>::from(geo); // (Geo-F) -> (Geo-Equ)
    let target_geo_equ_no_aberration = remove_aberration(target_geo_equ, jd);

    let helio_equ = Position::<Heliocentric, Equatorial, U>::from_vec3(target_geo_equ_no_aberration.as_vec3() + earth_helio_equ.as_vec3()); // Geocentric -> Heliocentric
    Position::<Heliocentric, F, U>::from(&helio_equ) // Equatorial -> F
}

impl<F: MutableFrame, U: Unit> Transform<Position<Heliocentric, F, U>> for Position<Geocentric, F, U>
where
    for<'a> Position<Heliocentric, Equatorial, U>: From<&'a Position<Heliocentric, Ecliptic, U>>, // Required by VSOP
    for<'a> Position<Geocentric, Equatorial, U>: From<&'a Position<Geocentric, F, U>>,   // Required by Aberration
    for<'a> Position<Heliocentric, F, U>: From<&'a Position<Heliocentric, Equatorial, U>>, // Required by Aberration
{
    fn transform(
        &self,
        jd: JulianDay,
    ) -> Position<Heliocentric, F, U> {
        geocentric_to_heliocentric(self, jd)
    }
}

impl<F: MutableFrame, U: Unit> Transform<Position<Heliocentric, F, U>> for Position<Barycentric, F, U>
where
    for<'a> Position<Barycentric, F, U>: From<&'a Position<Barycentric, Ecliptic, U>>,
{
    fn transform(
        &self,
        jd: JulianDay,
    ) -> Position<Heliocentric, F, U> {
        barycentric_to_heliocentric(self, jd)
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::coordinates::frames::*;
    use crate::coordinates::centers::*;
    use crate::bodies::solar_system::Sun;

    const EPSILON: f64 = 1e-9; // Precision tolerance for floating-point comparisons

    /// Helper function to compare floating-point values within a small tolerance
    fn approx_eq(a: f64, b: f64, epsilon: f64) -> bool {
        (a - b).abs() < epsilon
    }

    fn coords_approx_eq(a: &Position<impl ReferenceCenter, impl MutableFrame>,
                        b: &Position<impl ReferenceCenter, impl MutableFrame>,
                        epsilon: f64) -> bool {
        approx_eq(a.x(), b.x(), epsilon) &&
        approx_eq(a.y(), b.y(), epsilon) &&
        approx_eq(a.z(), b.z(), epsilon)
    }

    #[test] // Barycentric -> Heliocentric
    fn test_bary() {
        let sun_bary = Sun::vsop87e(JulianDay::J2000).get_position().clone();
        let sun_helio = Position::<Heliocentric, Ecliptic>::from(&sun_bary);
        let expected_sun_helio = Position::<Heliocentric, Ecliptic>::new(0.0, 0.0, 0.0);
        assert!(coords_approx_eq(&sun_helio, &expected_sun_helio, EPSILON), 
                "Sun in Heliocentric shall be (0,0,0). Current Value {:?}", sun_helio);
    }

    #[test] // Geocentric -> Heliocentric
    fn test_geo() {
        let sun_geo = Position::<Geocentric, Ecliptic>::from(&Sun::vsop87e(JulianDay::J2000).get_position().clone()); // Sun in Geocentric
        let sun_helio = Position::<Heliocentric, Ecliptic>::from(&sun_geo);
        let expected_sun_helio = Position::<Heliocentric, Ecliptic>::new(0.0, 0.0, 0.0);
        assert!(coords_approx_eq(&sun_helio, &expected_sun_helio, 1e-8),
                "Sun in Heliocentric shall be (0,0,0). Current Value {:?}", sun_helio);
    }
}
