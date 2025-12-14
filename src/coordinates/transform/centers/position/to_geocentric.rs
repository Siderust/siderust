use crate::astro::aberration::apply_aberration;
use crate::astro::JulianDate;
use crate::bodies::solar_system::Earth;
use crate::coordinates::transform::centers::TransformCenter;
use crate::coordinates::{
    cartesian::position::{Ecliptic, Equatorial, Position},
    centers::{Barycentric, Geocentric, Heliocentric},
    frames::MutableFrame,
    transform::TransformFrame,
};
use qtty::{AstronomicalUnits, LengthUnit, Quantity};

// Barycentric To Geocentric
impl<F: MutableFrame, U: LengthUnit> TransformCenter<Position<Geocentric, F, U>>
    for Position<Barycentric, F, U>
where
    Quantity<U>: From<AstronomicalUnits>,
    Position<Barycentric, F, U>: TransformFrame<Equatorial<U, Barycentric>>, // Required by Aberration
    Equatorial<U, Barycentric>: TransformFrame<Position<Barycentric, F, U>>, // Required by Aberration
    Equatorial<U, Geocentric>: TransformFrame<Position<Geocentric, F, U>>,
{
    fn to_center(&self, jd: JulianDate) -> Position<Geocentric, F, U> {
        let earth_bary_ecl_au = *Earth::vsop87e(jd).get_position();
        let earth_ecl = Ecliptic::<U, Barycentric>::new(
            earth_bary_ecl_au.x(),
            earth_bary_ecl_au.y(),
            earth_bary_ecl_au.z(),
        );

        let earth_equ: Equatorial<U, Barycentric> = earth_ecl.to_frame(); // (Bary-Ecl) -> (Bary-Equ)
        let bary_equ: Equatorial<U, Barycentric> = self.to_frame(); // (Bary-Any) -> (Bary-Equ)
        let geo_equ = Equatorial::<U>::from_vec3(bary_equ.as_vec3() - earth_equ.as_vec3()); // Barycentric -> Geocentric

        let gcrs = apply_aberration(geo_equ, jd);
        gcrs.to_frame() // Equatorial -> F
    }
}

// Heliocentric To Geocentric
impl<F: MutableFrame, U: LengthUnit> TransformCenter<Position<Geocentric, F, U>>
    for Position<Heliocentric, F, U>
where
    Quantity<U>: From<AstronomicalUnits>,
    Position<Heliocentric, F, U>: TransformFrame<Equatorial<U, Heliocentric>>, // Required by Aberration
    Equatorial<U, Heliocentric>: TransformFrame<Position<Heliocentric, F, U>>, // Required by Aberration
    Equatorial<U, Geocentric>: TransformFrame<Position<Geocentric, F, U>>,
{
    fn to_center(&self, jd: JulianDate) -> Position<Geocentric, F, U> {
        let earth_helio_ecl_au = *Earth::vsop87a(jd).get_position();
        let earth_ecl = Ecliptic::<U>::new(
            earth_helio_ecl_au.x(),
            earth_helio_ecl_au.y(),
            earth_helio_ecl_au.z(),
        );

        let earth_equ: Equatorial<U, Heliocentric> = earth_ecl.to_frame(); // (Helio-Ecl) -> (Helio-Equ)
        let helio_equ: Equatorial<U, Heliocentric> = self.to_frame(); // (Helio-any) -> (Helio-Equ)
        let geo_equ = Equatorial::<U>::from_vec3(helio_equ.as_vec3() - earth_equ.as_vec3()); // Heliocentric -> Geocentric

        let gcrs = apply_aberration(geo_equ, jd);
        gcrs.to_frame() // Equatorial -> any
    }
}

#[cfg(test)]
mod tests {
    use crate::astro::JulianDate;
    use crate::bodies::solar_system::Earth;
    use crate::coordinates::{cartesian, centers::*, spherical, transform::Transform};
    use crate::macros::{assert_cartesian_eq, assert_spherical_eq};
    use qtty::*;

    const EPSILON: f64 = 1e-9; // Precision tolerance for floating-point comparisons

    #[test] // Barycentric -> Geocentric
    fn test_bary_to_geo() {
        let earth_bary = *Earth::vsop87e(JulianDate::J2000).get_position();
        let earth_geo: cartesian::position::Ecliptic<Au, Geocentric> =
            earth_bary.transform(JulianDate::J2000);
        let expected_earth_geo = cartesian::position::Ecliptic::<Au, Geocentric>::CENTER;
        assert_cartesian_eq!(
            &earth_geo,
            &expected_earth_geo,
            EPSILON,
            "Earth in Geocentric shall be (0,0,0). Current Value {:?}",
            earth_geo
        );
    }

    #[test] // Heliocentric -> Geocentric
    fn test_helio_to_geo() {
        let earth_helio = *Earth::vsop87a(JulianDate::J2000).get_position();
        let earth_geo: cartesian::position::Ecliptic<Au, Geocentric> =
            earth_helio.transform(JulianDate::J2000);
        let expected_earth_geo = cartesian::position::Ecliptic::<Au, Geocentric>::CENTER;
        assert_cartesian_eq!(&earth_geo, &expected_earth_geo, EPSILON);
    }

    #[test] // ICRS -> GCRS
    fn test_icrs_to_gcrs() {
        const AU_PER_PC: f64 = 206_264.806_f64;
        const SIRIUS_PARALLAX: f64 = 0.37921_f64; // arcsec  (Hipparcos van Leeuwen 2007)
        let sirius_distance_au = (1.0 / SIRIUS_PARALLAX) * AU_PER_PC;

        let sirius_barycentric_spherical = spherical::position::ICRS::<Au>::new(
            Degrees::new(101.287_155_33),
            Degrees::new(-16.716_115_86),
            sirius_distance_au,
        );
        let expected_sirius_coordinates = spherical::position::GCRS::<Au>::new(
            Degrees::new(101.2846608),
            Degrees::new(-16.71925194),
            543933.225421,
        );

        let sirius_barycentric_cartesian =
            cartesian::position::ICRS::<Au>::from(&sirius_barycentric_spherical);
        let sirius_geocentric_cartesian: cartesian::position::GCRS<Au> = Transform::transform(
            &sirius_barycentric_cartesian,
            JulianDate::new(2460792.157638889),
        );
        let sirius_geocentric_spherical =
            spherical::position::GCRS::<Au>::from(&sirius_geocentric_cartesian);
        assert_spherical_eq!(
            sirius_geocentric_spherical,
            expected_sirius_coordinates,
            2e-4,
            "Sirius in Geocentric shall be {}. Current Value {}",
            expected_sirius_coordinates,
            sirius_geocentric_spherical
        );
    }
}
