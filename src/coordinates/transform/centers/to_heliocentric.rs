use crate::astro::aberration::remove_aberration;
use crate::bodies::solar_system::{Earth, Sun};
use crate::coordinates::transform::Transform;
use crate::coordinates::{
    cartesian::Position,
    centers::{Barycentric, Geocentric, Heliocentric},
    frames::{Ecliptic, Equatorial, MutableFrame},
};
use crate::units::{Quantity, AstronomicalUnits, JulianDay, LengthUnit};

impl<F: MutableFrame, U: LengthUnit> Transform<Position<Heliocentric, F, U>>
    for Position<Geocentric, F, U>
where
    Quantity<U>: From<AstronomicalUnits> + PartialEq + std::fmt::Debug,
    for<'a> Position<Heliocentric, Equatorial, U>: From<&'a Position<Heliocentric, Ecliptic, U>>, // Required by VSOP
    for<'a> Position<Geocentric, Equatorial, U>: From<&'a Position<Geocentric, F, U>>, // Required by Aberration
    for<'a> Position<Heliocentric, F, U>: From<&'a Position<Heliocentric, Equatorial, U>>, // Required by Aberration
{
    fn transform(&self, jd: JulianDay) -> Position<Heliocentric, F, U> {
        //geocentric_to_heliocentric(self, jd)
        // VSOP87 gives the Earth's position in AstronomicalUnits, so we need to convert to U
        let earth_helio_ecl_au = Earth::vsop87a(jd).get_position().clone();
        let earth_helio_ecl = Position::<Heliocentric, Ecliptic, U>::new(
            earth_helio_ecl_au.x().into(),
            earth_helio_ecl_au.y().into(),
            earth_helio_ecl_au.z().into(),
        );

        let earth_helio_equ = Position::<Heliocentric, Equatorial, U>::from(&earth_helio_ecl); // (Helio-Ecl) -> (Helio-Equ)

        let target_geo_equ = Position::<Geocentric, Equatorial, U>::from(self); // (Geo-F) -> (Geo-Equ)
        let target_geo_equ_no_aberration = remove_aberration(target_geo_equ, jd);

        let helio_equ = Position::<Heliocentric, Equatorial, U>::from_vec3(
            target_geo_equ_no_aberration.as_vec3() + earth_helio_equ.as_vec3(),
        ); // Geocentric -> Heliocentric
        Position::<Heliocentric, F, U>::from(&helio_equ) // Equatorial -> F
    }
}

impl<F: MutableFrame, U: LengthUnit> Transform<Position<Heliocentric, F, U>>
    for Position<Barycentric, F, U>
where
    Quantity<U>: From<AstronomicalUnits> + PartialEq + std::fmt::Debug,
    for<'a> Position<Barycentric, F, U>: From<&'a Position<Barycentric, Ecliptic, U>>,
{
    fn transform(&self, jd: JulianDay) -> Position<Heliocentric, F, U> {
        // Barycentric to Heliocentric
        // VSOP87 gives the Sun's position in AstronomicalUnits, so we need to convert to U
        let sun_bary_ecl_au = Sun::vsop87e(jd).get_position().clone();
        let sun_bary_ecl = Position::<Barycentric, Ecliptic, U>::new(
            sun_bary_ecl_au.x().into(),
            sun_bary_ecl_au.y().into(),
            sun_bary_ecl_au.z().into(),
        );

        let sun = Position::<Barycentric, F, U>::from(&sun_bary_ecl);
        Position::from_vec3(self.as_vec3() - sun.as_vec3())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bodies::solar_system::Sun;
    use crate::coordinates::centers::*;
    use crate::coordinates::frames::*;
    use crate::macros::assert_cartesian_eq;
    use crate::units::AU;

    const EPSILON: f64 = 1e-9; // Precision tolerance for floating-point comparisons

    #[test] // Barycentric -> Heliocentric
    fn test_bary() {
        let sun_bary = Sun::vsop87e(JulianDay::J2000).get_position().clone();
        let sun_helio = Position::<Heliocentric, Ecliptic, AU>::from(&sun_bary);
        let expected_sun_helio =
            Position::<Heliocentric, Ecliptic, AU>::new(0.0 * AU, 0.0 * AU, 0.0 * AU);
        assert_cartesian_eq!(
            &sun_helio,
            &expected_sun_helio,
            EPSILON,
            "Sun in Heliocentric shall be (0,0,0). Current Value {:?}",
            sun_helio
        );
    }

    #[test] // Geocentric -> Heliocentric
    fn test_geo() {
        let sun_geo = Position::<Geocentric, Ecliptic, AU>::from(
            &Sun::vsop87e(JulianDay::J2000).get_position().clone(),
        ); // Sun in Geocentric
        let sun_helio = Position::<Heliocentric, Ecliptic, AU>::from(&sun_geo);
        let expected_sun_helio =
            Position::<Heliocentric, Ecliptic, AU>::new(0.0 * AU, 0.0 * AU, 0.0 * AU);
        assert_cartesian_eq!(
            &sun_helio,
            &expected_sun_helio,
            1e-8,
            "Sun in Heliocentric shall be (0,0,0). Current Value {:?}",
            sun_helio
        );
    }
}
