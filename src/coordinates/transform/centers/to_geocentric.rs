use crate::astro::aberration::{apply_aberration, apply_aberration_to_direction};
use crate::bodies::solar_system::Earth;
use crate::coordinates::transform::Transform;
use crate::coordinates::{
    cartesian::Direction, cartesian::Position, centers::*, frames::*, spherical,
};
use crate::units::{AstronomicalUnits, JulianDay, LengthUnit, Quantity};

// Barycentric To Geocentric
impl<F: MutableFrame, U: LengthUnit> Transform<Position<Geocentric, F, U>>
    for Position<Barycentric, F, U>
where
    Quantity<U>: From<AstronomicalUnits> + PartialEq + std::fmt::Debug,
    for<'a> Position<Barycentric, Equatorial, U>: From<&'a Position<Barycentric, Ecliptic, U>>, // Required by VSOP
    for<'a> Position<Barycentric, Equatorial, U>: From<&'a Position<Barycentric, F, U>>, // Required by Aberration
    for<'a> Position<Geocentric, F, U>: From<&'a Position<Geocentric, Equatorial, U>>, // Required by Aberration
{
    fn transform(&self, jd: JulianDay) -> Position<Geocentric, F, U> {
        let earth_bary_ecl_au = Earth::vsop87e(jd).get_position().clone();
        let earth_ecl = Position::<Barycentric, Ecliptic, U>::new(
            earth_bary_ecl_au.x(),
            earth_bary_ecl_au.y(),
            earth_bary_ecl_au.z(),
        );

        let earth_equ = Position::<Barycentric, Equatorial, U>::from(&earth_ecl); // (Bary-Ecl) -> (Bary-Equ)
        let bary_equ = Position::<Barycentric, Equatorial, U>::from(self); // (Bary-F)   -> (Bary-Equ)
        let geo_equ = Position::<Geocentric, Equatorial, U>::from_vec3(
            bary_equ.as_vec3() - earth_equ.as_vec3(),
        ); // Barycentric -> Geocentric

        let gcrs = apply_aberration(geo_equ, jd);
        Position::<Geocentric, F, U>::from(&gcrs) // Equatorial -> F
    }
}

// Heliocentric To Geocentric
impl<F: MutableFrame, U: LengthUnit> Transform<Position<Geocentric, F, U>>
    for Position<Heliocentric, F, U>
where
    Quantity<U>: From<AstronomicalUnits> + PartialEq + std::fmt::Debug,
    for<'a> Position<Heliocentric, Equatorial, U>: From<&'a Position<Heliocentric, Ecliptic, U>>, // Required by VSOP
    for<'a> Position<Heliocentric, Equatorial, U>: From<&'a Position<Heliocentric, F, U>>, // Required by Aberration
    for<'a> Position<Geocentric, F, U>: From<&'a Position<Geocentric, Equatorial, U>>, // Required by Aberration
{
    fn transform(&self, jd: JulianDay) -> Position<Geocentric, F, U> {
        let earth_helio_ecl_au = Earth::vsop87a(jd).get_position().clone();
        let earth_ecl = Position::<Heliocentric, Ecliptic, U>::new(
            earth_helio_ecl_au.x(),
            earth_helio_ecl_au.y(),
            earth_helio_ecl_au.z(),
        );

        let earth_equ = Position::<Heliocentric, Equatorial, U>::from(&earth_ecl); // (Helio-Ecl) -> (Helio-Equ)
        let helio_equ = Position::<Heliocentric, Equatorial, U>::from(self); // (Helio-F)   -> (Helio-Equ)
        let geo_equ = Position::<Geocentric, Equatorial, U>::from_vec3(
            helio_equ.as_vec3() - earth_equ.as_vec3(),
        ); // Heliocentric -> Geocentric

        let gcrs = apply_aberration(geo_equ, jd);
        Position::<Geocentric, F, U>::from(&gcrs) // Equatorial -> F
    }
}

// ------------- If we transform TO a Geocentric Direction, we only need to apply aberration ------------------
impl<C, F> Transform<Direction<Geocentric, F>> for Direction<C, F>
where
    C: ReferenceCenter + NonGeocentric,
    F: ReferenceFrame + MutableFrame,
    Direction<Geocentric, F>: Transform<Direction<Geocentric, Equatorial>>, // Required by Aberration
    Direction<Geocentric, Equatorial>: Transform<Direction<Geocentric, F>>,
{
    #[inline]
    fn transform(&self, jd: JulianDay) -> Direction<Geocentric, F> {
        // 1. Convert to Geocentric Equatorial coordinates
        let geocentric = Direction::<Geocentric, F>::from_vec3(self.as_vec3());
        // 2. Transform to Geocentric Equatorial
        let equatorial: Direction<Geocentric, Equatorial> = geocentric.transform(jd);
        // 3. Apply aberration
        let aberrated = apply_aberration_to_direction(
            Direction::<Geocentric, Equatorial>::from_vec3(equatorial.as_vec3()),
            jd,
        );
        // 4. Recover target Frame
        Transform::<Direction<Geocentric, F>>::transform(&aberrated, jd)
    }
}

impl<C, F> Transform<spherical::Direction<Geocentric, F>> for spherical::Direction<C, F>
where
    C: ReferenceCenter + NonGeocentric,
    F: ReferenceFrame + MutableFrame,
    Direction<Geocentric, F>: Transform<Direction<Geocentric, Equatorial>>,
    Direction<Geocentric, Equatorial>: Transform<Direction<Geocentric, F>>,
{
    #[inline]
    fn transform(&self, jd: JulianDay) -> spherical::Direction<Geocentric, F> {
        let cart = self.to_cartesian();
        let cart_transformed = Transform::<Direction<Geocentric, F>>::transform(&cart, jd);
        cart_transformed.to_spherical()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bodies::solar_system::Earth;
    use crate::coordinates::*;
    use crate::macros::{assert_cartesian_eq, assert_spherical_eq};
    use crate::units::{Degrees, AstronomicalUnit};

    const EPSILON: f64 = 1e-9; // Precision tolerance for floating-point comparisons

    #[test] // Barycentric -> Geocentric
    fn test_bary_to_geo() {
        let earth_bary = Earth::vsop87e(JulianDay::J2000).get_position().clone();
        let earth_geo = Position::<Geocentric, Ecliptic, AstronomicalUnit>::from(&earth_bary);
        let expected_earth_geo = Position::<Geocentric, Ecliptic, AstronomicalUnit>::CENTER;
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
        let earth_helio: cartesian::Vector<Heliocentric, Ecliptic, crate::units::AstronomicalUnit> =
            Earth::vsop87a(JulianDay::J2000).get_position().clone();
        let earth_geo = Position::<Geocentric, Ecliptic, AstronomicalUnit>::from(&earth_helio);
        let expected_earth_geo = Position::<Geocentric, Ecliptic, AstronomicalUnit>::CENTER;
        assert_cartesian_eq!(&earth_geo, &expected_earth_geo, EPSILON);
    }

    #[test] // ICRS -> GCRS
    fn test_icrs_to_gcrs() {
        const AU_PER_PC: f64 = 206_264.806_f64;
        const SIRIUS_PARALLAX: f64 = 0.37921_f64; // arcsec  (Hipparcos van Leeuwen 2007)
        let sirius_distance_au = (1.0 / SIRIUS_PARALLAX) * AU_PER_PC;

        let sirius_barycentric_spherical = spherical::position::ICRS::<AstronomicalUnit>::new(
            Degrees::new(101.287_155_33),
            Degrees::new(-16.716_115_86),
            sirius_distance_au,
        );
        let expected_sirius_coordinates = spherical::position::GCRS::<AstronomicalUnit>::new(
            Degrees::new(101.2846608),
            Degrees::new(-16.71925194),
            543933.225421,
        );

        let sirius_barycentric_cartesian =
            cartesian::position::ICRS::<AstronomicalUnit>::from(&sirius_barycentric_spherical);
        let sirius_geocentric_cartesian: cartesian::position::GCRS<AstronomicalUnit> = Transform::transform(
            &sirius_barycentric_cartesian,
            JulianDay::new(2460792.157638889),
        );
        let sirius_geocentric_spherical =
            spherical::Position::<Geocentric, frames::ICRS, AstronomicalUnit>::from(&sirius_geocentric_cartesian);
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
