use crate::units::JulianDay;
use crate::bodies::solar_system::{Sun, Earth};
use crate::coordinates::{
    CartesianCoord,
    frames::{ReferenceFrame, Ecliptic},
    centers::{Heliocentric, Barycentric, Geocentric}
};
use crate::coordinates::transform::Transform;

pub fn barycentric_to_heliocentric<F: ReferenceFrame>(
    bary: &CartesianCoord<Barycentric, F>,
    jd: JulianDay
) -> CartesianCoord<Heliocentric, F>
where
    for<'a> CartesianCoord<Barycentric, F>: From<&'a CartesianCoord<Barycentric, Ecliptic>>,
{
    let sun = CartesianCoord::<Barycentric, F>::from(
        Sun::vsop87e(jd).get_position()
    );
    CartesianCoord::new(
        bary.x() - sun.x(),
        bary.y() - sun.y(),
        bary.z() - sun.z(),
    )
}

pub fn geocentric_to_heliocentric<F: ReferenceFrame>(
    geo: &CartesianCoord<Geocentric, F>,
    jd: JulianDay
) -> CartesianCoord<Heliocentric, F>
where
    for<'a> CartesianCoord<Heliocentric, F>: From<&'a CartesianCoord<Heliocentric, Ecliptic>>,
{
    let earth = CartesianCoord::<Heliocentric, F>::from(
        Earth::vsop87a(jd).get_position()
    );
    CartesianCoord::new(
        geo.x() + earth.x(),
        geo.y() + earth.y(),
        geo.z() + earth.z(),
    )
}

impl<F: ReferenceFrame> Transform<CartesianCoord<Heliocentric, F>> for CartesianCoord<Geocentric, F>
where
    for<'a> CartesianCoord<Heliocentric, F>: From<&'a CartesianCoord<Heliocentric, Ecliptic>>,
{
    fn transform(
        &self,
        jd: JulianDay,
    ) -> CartesianCoord<Heliocentric, F> {
        geocentric_to_heliocentric(self, jd)
    }
}

impl<F: ReferenceFrame> Transform<CartesianCoord<Heliocentric, F>> for CartesianCoord<Barycentric, F>
where
    for<'a> CartesianCoord<Barycentric, F>: From<&'a CartesianCoord<Barycentric, Ecliptic>>,
{
    fn transform(
        &self,
        jd: JulianDay,
    ) -> CartesianCoord<Heliocentric, F> {
        barycentric_to_heliocentric(self, jd)
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::coordinates::*;
    use crate::coordinates::frames::*;
    use crate::coordinates::centers::*;
    use crate::bodies::solar_system::Sun;

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

    #[test] // Barycentric -> Heliocentric
    fn test_bary() {
        let sun_bary = Sun::vsop87e(JulianDay::J2000).get_position().clone();
        let sun_helio = CartesianCoord::<Heliocentric, Ecliptic>::from(&sun_bary);
        let expected_sun_helio = CartesianCoord::<Heliocentric, Ecliptic>::new(0.0, 0.0, 0.0);
        assert!(coords_approx_eq(&sun_helio, &expected_sun_helio, EPSILON), 
                "Sun in Heliocentric shall be (0,0,0). Current Value {:?}", sun_helio);
    }

    #[test] // Geocentric -> Heliocentric
    fn test_geo() {
        let sun_geo = CartesianCoord::<Geocentric, Ecliptic>::from(&Sun::vsop87e(JulianDay::J2000).get_position().clone()); // Sun in Geocentric
        let sun_helio = CartesianCoord::<Heliocentric, Ecliptic>::from(&sun_geo);
        let expected_sun_helio = CartesianCoord::<Heliocentric, Ecliptic>::new(0.0, 0.0, 0.0);
        assert!(coords_approx_eq(&sun_helio, &expected_sun_helio, 5e-7),  // WARNING:  TOLERANCE TOO HIGH
                "Sun in Heliocentric shall be (0,0,0). Current Value {:?}", sun_helio);
    }
}
