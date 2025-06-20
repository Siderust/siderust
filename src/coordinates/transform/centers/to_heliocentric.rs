use crate::units::JulianDay;
use crate::bodies::solar_system::{Sun, Earth};
use crate::coordinates::{
    cartesian::Position,
    frames::{MutableFrame, Ecliptic, Equatorial},
    centers::{Heliocentric, Barycentric, Geocentric}
};
use crate::coordinates::transform::Transform;
use crate::astro::aberration::remove_aberration;

pub fn barycentric_to_heliocentric<F: MutableFrame>(
    bary: &Position<Barycentric, F>,
    jd: JulianDay
) -> Position<Heliocentric, F>
where
    for<'a> Position<Barycentric, F>: From<&'a Position<Barycentric, Ecliptic>>,
{
    let sun = Position::<Barycentric, F>::from(
        Sun::vsop87e(jd).get_position()
    );
    Position::new(
        bary.x() - sun.x(),
        bary.y() - sun.y(),
        bary.z() - sun.z(),
    )
}

pub fn geocentric_to_heliocentric<F: MutableFrame>(
    geo: &Position<Geocentric, F>,
    jd: JulianDay
) -> Position<Heliocentric, F>
where
    for<'a> Position<Heliocentric, Equatorial>: From<&'a Position<Heliocentric, Ecliptic>>, // Required by VSOP
    for<'a> Position<Geocentric, Equatorial>: From<&'a Position<Geocentric, F>>,   // Required by Aberration
    for<'a> Position<Heliocentric, F>: From<&'a Position<Heliocentric, Equatorial>>, // Required by Aberration
{
    let earth_helio_ecl = Earth::vsop87a(jd).get_position().clone();
    let earth_helio_equ = Position::<Heliocentric, Equatorial>::from(&earth_helio_ecl); // (Helio-Ecl) -> (Helio-Equ)

    let target_geo_equ  = Position::<Geocentric, Equatorial>::from(geo); // (Geo-F) -> (Geo-Equ)
    let target_geo_equ_no_aberration = remove_aberration(target_geo_equ, jd);

    let helio_equ = Position::<Heliocentric, Equatorial>::from_vec3(target_geo_equ_no_aberration.as_vec3() + earth_helio_equ.as_vec3()); // Geocentric -> Barycentric
    Position::<Heliocentric, F>::from(&helio_equ) // Equatorial -> F
}

impl<F: MutableFrame> Transform<Position<Heliocentric, F>> for Position<Geocentric, F>
where
    for<'a> Position<Heliocentric, Equatorial>: From<&'a Position<Heliocentric, Ecliptic>>, // Required by VSOP
    for<'a> Position<Geocentric, Equatorial>: From<&'a Position<Geocentric, F>>,   // Required by Aberration
    for<'a> Position<Heliocentric, F>: From<&'a Position<Heliocentric, Equatorial>>, // Required by Aberration
{
    fn transform(
        &self,
        jd: JulianDay,
    ) -> Position<Heliocentric, F> {
        geocentric_to_heliocentric(self, jd)
    }
}

impl<F: MutableFrame> Transform<Position<Heliocentric, F>> for Position<Barycentric, F>
where
    for<'a> Position<Barycentric, F>: From<&'a Position<Barycentric, Ecliptic>>,
{
    fn transform(
        &self,
        jd: JulianDay,
    ) -> Position<Heliocentric, F> {
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
