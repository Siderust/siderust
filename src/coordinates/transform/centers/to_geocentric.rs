use crate::bodies::solar_system::Earth;
use crate::units::{JulianDay, Unit};
use crate::coordinates::{
    cartesian::Position, cartesian::Direction,
    frames::*, centers::*,
    spherical
};
use crate::coordinates::transform::Transform;
use crate::astro::aberration::{apply_aberration, apply_aberration_to_direction};

pub fn barycentric_to_geocentric<F: MutableFrame, U: Unit>(
    bary: &Position<Barycentric, F, U>,
    jd: JulianDay
) -> Position<Geocentric, F, U>
where
    for<'a> Position<Barycentric, Equatorial, U>: From<&'a Position<Barycentric, Ecliptic, U>>, // Required by VSOP
    for<'a> Position<Barycentric, Equatorial, U>: From<&'a Position<Barycentric, F, U>>, // Required by Aberration
    for<'a> Position<Geocentric, F, U>: From<&'a Position<Geocentric, Equatorial, U>>, // Required by Aberration
{
    // VSOP87 gives the Earth's position in AstronomicalUnits, so we need to convert to U
    let earth_bary_ecl_au = Earth::vsop87e(jd).get_position().clone();
    let x: U = earth_bary_ecl_au.x().into();
    let y: U = earth_bary_ecl_au.y().into();
    let z: U = earth_bary_ecl_au.z().into();
    let earth_ecl = Position::<Barycentric, Ecliptic, U>::new(x, y, z);

    let earth_equ = Position::<Barycentric, Equatorial, U>::from(&earth_ecl); // (Bary-Ecl) -> (Bary-Equ)
    let bary_equ  = Position::<Barycentric, Equatorial, U>::from(bary);       // (Bary-F)   -> (Bary-Equ)
    let geo_equ   = Position::<Geocentric, Equatorial, U>::from_vec3(bary_equ.as_vec3() - earth_equ.as_vec3()); // Barycentric -> Geocentric

    let gcrs = apply_aberration(geo_equ, jd);
    Position::<Geocentric, F, U>::from(&gcrs) // Equatorial -> F
}

pub fn heliocentric_to_geocentric<F: MutableFrame, U: Unit>(
    helio: &Position<Heliocentric, F, U>,
    jd: JulianDay
) -> Position<Geocentric, F, U>
where
    for<'a> Position<Heliocentric, Equatorial, U>: From<&'a Position<Heliocentric, Ecliptic, U>>, // Required by VSOP
    for<'a> Position<Heliocentric, Equatorial, U>: From<&'a Position<Heliocentric, F, U>>, // Required by Aberration
    for<'a> Position<Geocentric, F, U>: From<&'a Position<Geocentric, Equatorial, U>>, // Required by Aberration
{
    // VSOP87 gives the Earth's position in AstronomicalUnits, so we need to convert to U
    let earth_helio_ecl_au = Earth::vsop87a(jd).get_position().clone();
    let x: U = earth_helio_ecl_au.x().into();
    let y: U = earth_helio_ecl_au.y().into();
    let z: U = earth_helio_ecl_au.z().into();
    let earth_ecl = Position::<Barycentric, Ecliptic, U>::new(x, y, z);

    let earth_equ = Position::<Heliocentric, Equatorial, U>::from(&earth_ecl); // (Helio-Ecl) -> (Helio-Equ)
    let helio_equ = Position::<Heliocentric, Equatorial, U>::from(helio); // (Helio-F)   -> (Helio-Equ)
    let geo_equ   = Position::<Geocentric, Equatorial, U>::from_vec3(helio_equ.as_vec3() - earth_equ.as_vec3()); // Heliocentric -> Geocentric

    let gcrs = apply_aberration(geo_equ, jd);
    Position::<Geocentric, F, U>::from(&gcrs) // Equatorial -> F
}

impl<F: MutableFrame, U: Unit> Transform<Position<Geocentric, F, U>> for Position<Barycentric, F, U>
where
    for<'a> Position<Barycentric, Equatorial, U>: From<&'a Position<Barycentric, Ecliptic, U>>, // Required by VSOP
    for<'a> Position<Barycentric, Equatorial, U>: From<&'a Position<Barycentric, F, U>>, // Required by Aberration
    for<'a> Position<Geocentric, F, U>: From<&'a Position<Geocentric, Equatorial, U>>, // Required by Aberration
{
    fn transform(
        &self,
        jd: JulianDay,
    ) -> Position<Geocentric, F, U> {
        barycentric_to_geocentric(self, jd)
    }
}

impl<F: MutableFrame, U: Unit> Transform<Position<Geocentric, F, U>> for Position<Heliocentric, F, U>
where
    for<'a> Position<Heliocentric, Equatorial, U>: From<&'a Position<Heliocentric, Ecliptic, U>>, // Required by VSOP
    for<'a> Position<Heliocentric, Equatorial, U>: From<&'a Position<Heliocentric, F, U>>, // Required by Aberration
    for<'a> Position<Geocentric, F, U>: From<&'a Position<Geocentric, Equatorial, U>>, // Required by Aberration
{
    fn transform(
        &self,
        jd: JulianDay,
    ) -> Position<Geocentric, F, U> {
        heliocentric_to_geocentric(self, jd)
    }
}

// ------------- If we transform TO a Geocentric Direction, we only need to apply aberration ------------------
impl<C, F, U> Transform<Direction<Geocentric, F, U>> for Direction<C, F, U>
where
    C: ReferenceCenter + NonGeocentric,
    F: ReferenceFrame + MutableFrame,
    U: Unit,
    Direction<Geocentric, F, U>: Transform<Direction<Geocentric, Equatorial, U>>, // Required by Aberration
    Direction<Geocentric, Equatorial, U>: Transform<Direction<Geocentric, F, U>>,
{
    #[inline]
    fn transform(&self, jd: JulianDay) -> Direction<Geocentric, F, U> {
        // 1. Convert to Geocentric Equatorial coordinates
        let geocentric = Direction::<Geocentric, F, U>::from_vec3(
            self.as_vec3()
        );
        // 2. Transform to Geocentric Equatorial
        let equatorial: Direction<Geocentric, Equatorial, U> = geocentric.transform(jd);
        // 3. Apply aberration
        let aberrated = apply_aberration_to_direction(
            Direction::<Geocentric, Equatorial, U>::from_vec3(equatorial.as_vec3()), jd
        );
        // 4. Recover target Frame
        Transform::<Direction<Geocentric, F, U>>::transform(&aberrated, jd)
    }
}

impl<C, F, U> Transform<spherical::Direction<Geocentric, F, U>> for spherical::Direction<C, F, U>
where
    C: ReferenceCenter + NonGeocentric,
    F: ReferenceFrame + MutableFrame,
    U: Unit,
    Direction<Geocentric, F, U>: Transform<Direction<Geocentric, Equatorial, U>>,
    Direction<Geocentric, Equatorial, U>: Transform<Direction<Geocentric, F, U>>,
{
    #[inline]
    fn transform(&self, jd: JulianDay) -> spherical::Direction<Geocentric, F, U> {
        let cart = self.to_cartesian();
        let cart_transformed = Transform::<Direction<Geocentric, F, U>>::transform(&cart, jd);
        cart_transformed.to_spherical()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coordinates::*;
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
