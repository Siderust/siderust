use crate::bodies::solar_system::Earth;
use crate::units::JulianDay;
use crate::coordinates::{
    CartesianCoord,
    frames::{ReferenceFrame, Ecliptic},
    centers::{Heliocentric, Barycentric, Geocentric}
};
use crate::coordinates::transform::Transform;
use crate::astro::aberration::ecl_aberration;

pub fn barycentric_to_geocentric<F: ReferenceFrame>(
    bary: &CartesianCoord<Barycentric, F>,
    jd: JulianDay
) -> CartesianCoord<Geocentric, F>
where
    for<'a> CartesianCoord<Barycentric, Ecliptic>: From<&'a CartesianCoord<Barycentric, F>>, // Required by VSOP
    for<'a> CartesianCoord<Geocentric,  F>: From<&'a CartesianCoord<Geocentric, Ecliptic>>,  // Required by Aberration
{
    let earth_ecl = Earth::vsop87e(jd).get_position().clone();
    let bary_ecl  = CartesianCoord::<Barycentric, Ecliptic>::from(bary);
    let geo_ecl   = CartesianCoord::<Geocentric, Ecliptic>::from_vec3(bary_ecl.as_vec3() - earth_ecl.as_vec3());

    CartesianCoord::<Geocentric, F>::from(
        &ecl_aberration(geo_ecl, jd)
    )
}

pub fn heliocentric_to_geocentric<F: ReferenceFrame>(
    helio: &CartesianCoord<Heliocentric, F>,
    jd: JulianDay
) -> CartesianCoord<Geocentric, F>
where
    for<'a> CartesianCoord<Heliocentric, F>: From<&'a CartesianCoord<Heliocentric, Ecliptic>>,
{
    let earth = CartesianCoord::<Heliocentric, F>::from(
        Earth::vsop87a(jd).get_position()
    );
    CartesianCoord::from_vec3(helio.as_vec3() - earth.as_vec3())
}

impl<F: ReferenceFrame> Transform<CartesianCoord<Geocentric, F>> for CartesianCoord<Barycentric, F>
where
    for<'a> CartesianCoord<Barycentric, Ecliptic>: From<&'a CartesianCoord<Barycentric, F>>, // Required by VSOP
    for<'a> CartesianCoord<Geocentric,  F>: From<&'a CartesianCoord<Geocentric, Ecliptic>>,  // Required by Aberration
{
    fn transform(
        &self,
        jd: JulianDay,
    ) -> CartesianCoord<Geocentric, F> {
        barycentric_to_geocentric(self, jd)
    }
}

impl<F: ReferenceFrame> Transform<CartesianCoord<Geocentric, F>> for CartesianCoord<Heliocentric, F>
where
    for<'a> CartesianCoord<Heliocentric, F>: From<&'a CartesianCoord<Heliocentric, Ecliptic>>,
{
    fn transform(
        &self,
        jd: JulianDay,
    ) -> CartesianCoord<Geocentric, F> {
        heliocentric_to_geocentric(self, jd)
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::coordinates::*;
    use crate::coordinates::frames::*;
    use crate::coordinates::centers::*;
    use crate::bodies::solar_system::Earth;

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

    #[test] // Barycentric -> Geocentric
    fn test_helio() {
        let earth_bary = Earth::vsop87e(JulianDay::J2000).get_position().clone();
        let earth_geo = CartesianCoord::<Geocentric, Ecliptic>::from(&earth_bary);
        let expected_earth_geo = CartesianCoord::<Geocentric, Ecliptic>::new(0.0, 0.0, 0.0);
        assert!(coords_approx_eq(&earth_geo, &expected_earth_geo, EPSILON), 
                "Earth in Geocentric shall be (0,0,0). Current Value {:?}", earth_geo);
    }

    #[test] // Heliocentric -> Geocentric
    fn test_eart_bary() {
        let earth_helio = Earth::vsop87a(JulianDay::J2000).get_position().clone();
        let earth_geo = CartesianCoord::<Geocentric, Ecliptic>::from(&earth_helio);
        let expected_earth_geo = CartesianCoord::<Geocentric, Ecliptic>::new(0.0, 0.0, 0.0);
        assert!(coords_approx_eq(&earth_geo, &expected_earth_geo, EPSILON), 
                "Earth in Geocentric shall be (0,0,0). Current Value {:?}", earth_geo);
    }
}
