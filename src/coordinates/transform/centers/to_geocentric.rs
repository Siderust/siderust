use crate::bodies::solar_system::Earth;
use crate::units::JulianDay;
use crate::coordinates::{
    CartesianCoord,
    frames::{ReferenceFrame, Ecliptic, Equatorial},
    centers::{Heliocentric, Barycentric, Geocentric}
};
use crate::coordinates::transform::Transform;
use crate::astro::aberration::apply_aberration;

pub fn barycentric_to_geocentric<F: ReferenceFrame>(
    bary: &CartesianCoord<Barycentric, F>,
    jd: JulianDay
) -> CartesianCoord<Geocentric, F>
where
    for<'a> CartesianCoord<Barycentric, Equatorial>: From<&'a CartesianCoord<Barycentric, Ecliptic>>, // Required by VSOP
    for<'a> CartesianCoord<Barycentric, Equatorial>: From<&'a CartesianCoord<Barycentric, F>>, // Required by Aberration
    for<'a> CartesianCoord<Geocentric, F>: From<&'a CartesianCoord<Geocentric, Equatorial>>, // Required by Aberration
{
    let earth_ecl = Earth::vsop87e(jd).get_position().clone(); // Barycentric Ecliptic Earth
    let earth_equ = CartesianCoord::<Barycentric, Equatorial>::from(&earth_ecl); // (Bary-Ecl) -> (Bary-Equ)
    let bary_equ  = CartesianCoord::<Barycentric, Equatorial>::from(bary);       // (Bary-F)   -> (Bary-Equ)
    let geo_equ    = CartesianCoord::<Geocentric, Equatorial>::from_vec3(bary_equ.as_vec3() - earth_equ.as_vec3()); // Barycentric -> Geocentric

    let gcrs = apply_aberration(geo_equ, jd);
    CartesianCoord::<Geocentric, F>::from(&gcrs) // Equatorial -> F
}

pub fn heliocentric_to_geocentric<F: ReferenceFrame>(
    helio: &CartesianCoord<Heliocentric, F>,
    jd: JulianDay
) -> CartesianCoord<Geocentric, F>
where
    for<'a> CartesianCoord<Heliocentric, Equatorial>: From<&'a CartesianCoord<Heliocentric, Ecliptic>>, // Required by VSOP
    for<'a> CartesianCoord<Heliocentric, Equatorial>: From<&'a CartesianCoord<Heliocentric, F>>, // Required by Aberration
    for<'a> CartesianCoord<Geocentric, F>: From<&'a CartesianCoord<Geocentric, Equatorial>>, // Required by Aberration
{
    let earth_ecl = Earth::vsop87a(jd).get_position().clone(); // Heliocentric Ecliptic Earth
    let earth_equ = CartesianCoord::<Heliocentric, Equatorial>::from(&earth_ecl); // (Helio-Ecl) -> (Helio-Equ)
    let helio_equ = CartesianCoord::<Heliocentric, Equatorial>::from(helio); // (Helio-F)   -> (Helio-Equ)
    let geo_equ     = CartesianCoord::<Geocentric, Equatorial>::from_vec3(helio_equ.as_vec3() - earth_equ.as_vec3()); // Heliocentric -> Geocentric

    let gcrs = apply_aberration(geo_equ, jd);
    CartesianCoord::<Geocentric, F>::from(&gcrs) // Equatorial -> F
}

impl<F: ReferenceFrame> Transform<CartesianCoord<Geocentric, F>> for CartesianCoord<Barycentric, F>
where
    for<'a> CartesianCoord<Barycentric, Equatorial>: From<&'a CartesianCoord<Barycentric, Ecliptic>>, // Required by VSOP
    for<'a> CartesianCoord<Barycentric, Equatorial>: From<&'a CartesianCoord<Barycentric, F>>, // Required by Aberration
    for<'a> CartesianCoord<Geocentric, F>: From<&'a CartesianCoord<Geocentric, Equatorial>>, // Required by Aberration
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
    for<'a> CartesianCoord<Heliocentric, Equatorial>: From<&'a CartesianCoord<Heliocentric, Ecliptic>>, // Required by VSOP
    for<'a> CartesianCoord<Heliocentric, Equatorial>: From<&'a CartesianCoord<Heliocentric, F>>, // Required by Aberration
    for<'a> CartesianCoord<Geocentric, F>: From<&'a CartesianCoord<Geocentric, Equatorial>>, // Required by Aberration
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
    use crate::units::Degrees;

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

    fn sph_coords_approx_eq(a: &SphericalCoord<impl ReferenceCenter, impl ReferenceFrame>,
                            b: &SphericalCoord<impl ReferenceCenter, impl ReferenceFrame>,
                            epsilon: f64) -> bool {
        approx_eq(a.polar.as_f64(), b.polar.as_f64(), epsilon) &&
        approx_eq(a.azimuth.as_f64(), b.azimuth.as_f64(), epsilon) &&
        approx_eq(a.radial_distance, b.radial_distance, epsilon)
    }

    #[test] // Barycentric -> Geocentric
    fn test_bary_to_geo() {
        let earth_bary = Earth::vsop87e(JulianDay::J2000).get_position().clone();
        let earth_geo = CartesianCoord::<Geocentric, Ecliptic>::from(&earth_bary);
        let expected_earth_geo = CartesianCoord::<Geocentric, Ecliptic>::new(0.0, 0.0, 0.0);
        assert!(coords_approx_eq(&earth_geo, &expected_earth_geo, EPSILON), 
                "Earth in Geocentric shall be (0,0,0). Current Value {:?}", earth_geo);
    }

    #[test] // Heliocentric -> Geocentric
    fn test_helio_to_geo() {
        let earth_helio = Earth::vsop87a(JulianDay::J2000).get_position().clone();
        let earth_geo = CartesianCoord::<Geocentric, Ecliptic>::from(&earth_helio);
        let expected_earth_geo = CartesianCoord::<Geocentric, Ecliptic>::new(0.0, 0.0, 0.0);
        assert!(coords_approx_eq(&earth_geo, &expected_earth_geo, EPSILON), 
                "Earth in Geocentric shall be (0,0,0). Current Value {:?}", earth_geo);
    }

    #[test] // ICRS -> GCRS
    fn test_icrs_to_gcrs() {
        let sirius = SphericalCoord::<Barycentric, frames::ICRS>::new(Degrees::new(101.28715533), Degrees::new(-16.71611586), 1.0);
        let sirius_geo = barycentric_to_geocentric(&sirius.to_cartesian(), JulianDay::new(2460792.157638889)).to_spherical();
        let expected_geo = SphericalCoord::<Geocentric, frames::ICRS>::new(Degrees::new(101.2846608), Degrees::new(-16.71925194), 1.0);
        assert!(sph_coords_approx_eq(&sirius_geo, &expected_geo, EPSILON), 
                "Sirius in Geocentric shall be {}. Current Value {}", expected_geo, sirius_geo);
    }
}
