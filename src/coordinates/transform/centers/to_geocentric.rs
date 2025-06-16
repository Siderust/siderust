use crate::bodies::solar_system::Earth;
use crate::units::JulianDay;
use crate::coordinates::{
    cartesian::Position,
    frames::{ReferenceFrame, Ecliptic, Equatorial},
    centers::{Heliocentric, Barycentric, Geocentric}
};
use crate::coordinates::transform::Transform;
use crate::astro::aberration::apply_aberration;

pub fn barycentric_to_geocentric<F: ReferenceFrame>(
    bary: &Position<Barycentric, F>,
    jd: JulianDay
) -> Position<Geocentric, F>
where
    for<'a> Position<Barycentric, Equatorial>: From<&'a Position<Barycentric, Ecliptic>>, // Required by VSOP
    for<'a> Position<Barycentric, Equatorial>: From<&'a Position<Barycentric, F>>, // Required by Aberration
    for<'a> Position<Geocentric, F>: From<&'a Position<Geocentric, Equatorial>>, // Required by Aberration
{
    let earth_ecl = Earth::vsop87e(jd).get_position().clone(); // Barycentric Ecliptic Earth
    let earth_equ = Position::<Barycentric, Equatorial>::from(&earth_ecl); // (Bary-Ecl) -> (Bary-Equ)
    let bary_equ  = Position::<Barycentric, Equatorial>::from(bary);       // (Bary-F)   -> (Bary-Equ)
    let geo_equ    = Position::<Geocentric, Equatorial>::from_vec3(bary_equ.as_vec3() - earth_equ.as_vec3()); // Barycentric -> Geocentric

    let gcrs = apply_aberration(geo_equ, jd);
    Position::<Geocentric, F>::from(&gcrs) // Equatorial -> F
}

pub fn heliocentric_to_geocentric<F: ReferenceFrame>(
    helio: &Position<Heliocentric, F>,
    jd: JulianDay
) -> Position<Geocentric, F>
where
    for<'a> Position<Heliocentric, Equatorial>: From<&'a Position<Heliocentric, Ecliptic>>, // Required by VSOP
    for<'a> Position<Heliocentric, Equatorial>: From<&'a Position<Heliocentric, F>>, // Required by Aberration
    for<'a> Position<Geocentric, F>: From<&'a Position<Geocentric, Equatorial>>, // Required by Aberration
{
    let earth_ecl = Earth::vsop87a(jd).get_position().clone(); // Heliocentric Ecliptic Earth
    let earth_equ = Position::<Heliocentric, Equatorial>::from(&earth_ecl); // (Helio-Ecl) -> (Helio-Equ)
    let helio_equ = Position::<Heliocentric, Equatorial>::from(helio); // (Helio-F)   -> (Helio-Equ)
    let geo_equ     = Position::<Geocentric, Equatorial>::from_vec3(helio_equ.as_vec3() - earth_equ.as_vec3()); // Heliocentric -> Geocentric

    let gcrs = apply_aberration(geo_equ, jd);
    Position::<Geocentric, F>::from(&gcrs) // Equatorial -> F
}

impl<F: ReferenceFrame> Transform<Position<Geocentric, F>> for Position<Barycentric, F>
where
    for<'a> Position<Barycentric, Equatorial>: From<&'a Position<Barycentric, Ecliptic>>, // Required by VSOP
    for<'a> Position<Barycentric, Equatorial>: From<&'a Position<Barycentric, F>>, // Required by Aberration
    for<'a> Position<Geocentric, F>: From<&'a Position<Geocentric, Equatorial>>, // Required by Aberration
{
    fn transform(
        &self,
        jd: JulianDay,
    ) -> Position<Geocentric, F> {
        barycentric_to_geocentric(self, jd)
    }
}

impl<F: ReferenceFrame> Transform<Position<Geocentric, F>> for Position<Heliocentric, F>
where
    for<'a> Position<Heliocentric, Equatorial>: From<&'a Position<Heliocentric, Ecliptic>>, // Required by VSOP
    for<'a> Position<Heliocentric, Equatorial>: From<&'a Position<Heliocentric, F>>, // Required by Aberration
    for<'a> Position<Geocentric, F>: From<&'a Position<Geocentric, Equatorial>>, // Required by Aberration
{
    fn transform(
        &self,
        jd: JulianDay,
    ) -> Position<Geocentric, F> {
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
    use crate::macros::{assert_cartesian_eq, assert_spherical_eq};

    const EPSILON: f64 = 1e-9; // Precision tolerance for floating-point comparisons

    #[test] // Barycentric -> Geocentric
    fn test_bary_to_geo() {
        let earth_bary = Earth::vsop87e(JulianDay::J2000).get_position().clone();
        let earth_geo = Position::<Geocentric, Ecliptic>::from(&earth_bary);
        let expected_earth_geo = Position::<Geocentric, Ecliptic>::new(0.0, 0.0, 0.0);
        assert_cartesian_eq!(&earth_geo, &expected_earth_geo, EPSILON, "Earth in Geocentric shall be (0,0,0). Current Value {:?}", earth_geo);
    }

    #[test] // Heliocentric -> Geocentric
    fn test_helio_to_geo() {
        let earth_helio = Earth::vsop87a(JulianDay::J2000).get_position().clone();
        let earth_geo = Position::<Geocentric, Ecliptic>::from(&earth_helio);
        let expected_earth_geo = Position::<Geocentric, Ecliptic>::new(0.0, 0.0, 0.0);
        assert_cartesian_eq!(&earth_geo, &expected_earth_geo, EPSILON, "Earth in Geocentric shall be (0,0,0). Current Value {:?}", earth_geo);
    }

    #[test] // ICRS -> GCRS
    fn test_icrs_to_gcrs() {
        const AU_PER_PC: f64 = 206_264.806_f64;
        const SIRIUS_PARALLAX: f64 = 0.37921_f64;          // arcsec  (Hipparcos van Leeuwen 2007)
        let sirius_distance_au = (1.0 / SIRIUS_PARALLAX) * AU_PER_PC;

        let sirius_barycentric_spherical = spherical::Position::<Barycentric, frames::ICRS>::new(
            Degrees::new(101.287_155_33),
            Degrees::new(-16.716_115_86),
            sirius_distance_au
        );
        let expected_sirius_coordinates = spherical::Position::<Geocentric, frames::ICRS>::new(
            Degrees::new(101.2846608),
            Degrees::new(-16.71925194),
            543933.225421
        );

        let sirius_barycentric_cartesian = Position::<Barycentric, frames::ICRS>::from(&sirius_barycentric_spherical);
        let sirius_geocentric_cartesian = barycentric_to_geocentric(&sirius_barycentric_cartesian, JulianDay::new(2460792.157638889));
        let sirius_geocentric_spherical = spherical::Position::<Geocentric, frames::ICRS>::from(&sirius_geocentric_cartesian);
        assert_spherical_eq!(sirius_geocentric_spherical, expected_sirius_coordinates, 2e-4, "Sirius in Geocentric shall be {}. Current Value {}", expected_sirius_coordinates, sirius_geocentric_spherical);
    }

}
