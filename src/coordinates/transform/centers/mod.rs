pub mod to_barycentric;
pub mod to_heliocentric;
pub mod to_geocentric;

use crate::units::JulianDay;
use crate::coordinates::{
    frames::*, centers::*,
    cartesian, spherical,
};
use crate::coordinates::transform::Transform;
use crate::astro::aberration::{
    apply_aberration_to_direction,
    remove_aberration_from_direction,
};

impl<C: ReferenceCenter> Transform<cartesian::Direction<Geocentric, Equatorial>> for cartesian::Direction<C, Equatorial>
where
    C: NonGeocentric
{
    #[inline]
    fn transform(&self, jd: JulianDay) -> cartesian::Direction<Geocentric, Equatorial> {
        apply_aberration_to_direction(
            cartesian::Direction::from_vec3(self.as_vec3()), jd
        )
    }
}

impl<C: ReferenceCenter> Transform<cartesian::Direction<C, Equatorial>> for cartesian::Direction<Geocentric, Equatorial>
where
    C: NonGeocentric
{
    #[inline]
    fn transform(&self, jd: JulianDay) -> cartesian::Direction<C, Equatorial> {
        cartesian::Direction::<C, Equatorial>::from_vec3(
            remove_aberration_from_direction(*self, jd).as_vec3()
        )
    }
}

impl<C: ReferenceCenter> Transform<spherical::Direction<Geocentric, Equatorial>> for spherical::Direction<C, Equatorial>
where
    C: NonGeocentric
{
    #[inline]
    fn transform(&self, jd: JulianDay) -> spherical::Direction<Geocentric, Equatorial> {
        let cart = self.to_cartesian();
        let catr_trasnformed: cartesian::Direction<Geocentric, Equatorial> = cart.transform(jd);
        catr_trasnformed.to_spherical()
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::coordinates::centers::{Geocentric, Heliocentric, Barycentric};
    use crate::coordinates::spherical::Direction;
    use crate::units::Degrees;

    #[test]
    fn test_from_heliocentric_to_geocentric() {
        let dir = Direction::<Heliocentric, Equatorial>::new(Degrees::new(120.0), Degrees::new(-30.0));
        let jd = crate::units::JulianDay::J2000;
        let transformed: Direction<Geocentric, Equatorial> = dir.transform(jd);
        assert!(!transformed.polar.as_f64().is_nan(), "Polar should not be NaN");
        assert!(!transformed.azimuth.as_f64().is_nan(), "Azimuth should not be NaN");
    }

        #[test]
    fn test_from_barycentric_to_geocentric() {
        let dir = Direction::<Barycentric, Equatorial>::new(Degrees::new(120.0), Degrees::new(-30.0));
        let jd = crate::units::JulianDay::J2000;
        let transformed: Direction<Geocentric, Equatorial> = dir.transform(jd);
        assert!(!transformed.polar.as_f64().is_nan(), "Polar should not be NaN");
        assert!(!transformed.azimuth.as_f64().is_nan(), "Azimuth should not be NaN");
    }
/*

    #[test]
    fn test_from_geocentric_to_heliocentric() {
        let dir = Direction::<Geocentric, Equatorial>::new(Degrees::new(10.0), Degrees::new(20.0));
        let jd = crate::units::JulianDay::J2000;
        let transformed: Direction<Heliocentric, Equatorial> = dir.transform(jd);
        // Should be close in direction, but not identical due to aberration, so just check type and not NaN
        assert!(!transformed.polar.as_f64().is_nan(), "Polar should not be NaN");
        assert!(!transformed.azimuth.as_f64().is_nan(), "Azimuth should not be NaN");
    }


    #[test]
    fn test_barycentric_to_heliocentric() {
        let dir = Direction::<Barycentric, ICRS>::new(Degrees::new(200.0), Degrees::new(45.0));
        let jd = crate::units::JulianDay::J2000;
        let transformed: Direction<Heliocentric, ICRS> = dir.transform(jd);
        assert!(!transformed.polar.as_f64().is_nan(), "Polar should not be NaN");
        assert!(!transformed.azimuth.as_f64().is_nan(), "Azimuth should not be NaN");
    }

    #[test]
    fn test_heliocentric_to_barycentric() {
        let dir = Direction::<Heliocentric, ICRS>::new(Degrees::new(100.0), Degrees::new(-10.0));
        let jd = crate::units::JulianDay::J2000;
        let transformed: Direction<Barycentric, ICRS> = dir.transform(jd);
        assert!(!transformed.polar.as_f64().is_nan(), "Polar should not be NaN");
        assert!(!transformed.azimuth.as_f64().is_nan(), "Azimuth should not be NaN");
    }*/
}
