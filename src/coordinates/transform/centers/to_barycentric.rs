use crate::units::JulianDay;
use crate::bodies::solar_system::{Sun, Earth};
use crate::coordinates::{
    frames::{ReferenceFrame, Ecliptic},
    centers::{Barycentric, Heliocentric, Geocentric},
    CartesianCoord
};
use crate::coordinates::transform::Transform;

pub fn heliocentric_to_barycentric<F: ReferenceFrame>(
    helio: &CartesianCoord<Heliocentric, F>,
    jd: JulianDay
) -> CartesianCoord<Barycentric, F>
where
    for<'a> CartesianCoord<Barycentric, F>: From<&'a CartesianCoord<Barycentric, Ecliptic>>,
{
    let sun = CartesianCoord::<Barycentric, F>::from(
        Sun::vsop87e(jd).get_position()
    );
    CartesianCoord::new(
        helio.x() + sun.x(),
        helio.y() + sun.y(),
        helio.z() + sun.z(),
    )
}

pub fn geocentric_to_barycentric<F: ReferenceFrame>(
    geo: &CartesianCoord<Geocentric, F>,
    jd: JulianDay
) -> CartesianCoord<Barycentric, F>
where
    for<'a> CartesianCoord<Barycentric, F>: From<&'a CartesianCoord<Barycentric, Ecliptic>>,
{
    let earth = CartesianCoord::<Barycentric, F>::from(
        Earth::vsop87e(jd).get_position()
    );
    CartesianCoord::new(
        geo.x() + earth.x(),
        geo.y() + earth.y(),
        geo.z() + earth.z(),
    )
}


impl<F: ReferenceFrame> Transform<CartesianCoord<Barycentric, F>> for CartesianCoord<Heliocentric, F>
where
    for<'a> CartesianCoord<Barycentric, F>: From<&'a CartesianCoord<Barycentric, Ecliptic>>,
{
    fn transform(
        &self,
        jd: JulianDay,
    ) -> CartesianCoord<Barycentric, F> {
        heliocentric_to_barycentric(self, jd)
    }
}

impl<F: ReferenceFrame>  Transform<CartesianCoord<Barycentric, F>> for CartesianCoord<Geocentric, F>
where
    for<'a> CartesianCoord<Barycentric, F>: From<&'a CartesianCoord<Barycentric, Ecliptic>>,
{

    fn transform(
        &self,
        jd: JulianDay,
    ) -> CartesianCoord<Barycentric, F> {
        geocentric_to_barycentric(self, jd)
    }
}



#[cfg(test)]
mod tests {
    use super::*;
    use crate::coordinates::*;
    use crate::coordinates::frames::*;
    use crate::coordinates::centers::*;
    use crate::bodies::solar_system::{Sun, Earth};

    const EPSILON: f64 = 1e-9; // Precision tolerance for floating-point comparisons

    /// Helper function to compare floating-point values within a small tolerance
    fn approx_eq(a: f64, b: f64, epsilon: f64) -> bool {
        (a - b).abs() < epsilon
    }

    fn coords_approx_eq(a: &CartesianCoord<impl ReferenceCenter, impl ReferenceFrame>,
                        b: &CartesianCoord<impl ReferenceCenter, impl ReferenceFrame>,
                        epsilon: f64) -> bool {
        approx_eq(a.x(), b.x(), epsilon) &&
        approx_eq(a.y(), b.y(), epsilon) &&
        approx_eq(a.z(), b.z(), epsilon)
    }

    #[test] // Heliocentric -> Barycentric
    fn test_helio() {
        let sun_helio = CartesianCoord::<Heliocentric, Ecliptic>::new(0.0, 0.0, 0.0);
        let sun_bary = CartesianCoord::<Barycentric, Ecliptic>::from(&sun_helio);
        let expected_sun_bary = Sun::vsop87e(JulianDay::J2000).get_position().clone();
        assert!(coords_approx_eq(&sun_bary, &expected_sun_bary, EPSILON), 
                "Sun in Barycentric shall be {:?}. Current Value {:?}", expected_sun_bary, sun_bary);
    }

    #[test] // Geocentric -> Barycentric
    fn test_geo() {
        let earth_geo = CartesianCoord::<Geocentric, Ecliptic>::new(0.0, 0.0, 0.0);
        let earth_bary = CartesianCoord::<Barycentric, Ecliptic>::from(&earth_geo);
        let expected_earth_bary = Earth::vsop87e(JulianDay::J2000).get_position().clone();
        assert!(coords_approx_eq(&earth_bary, &expected_earth_bary, EPSILON), 
                "Earth in Geocentric shall be {:?}. Current Value {:?}", expected_earth_bary, earth_bary);
    }

}
