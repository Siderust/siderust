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

impl<C: ReferenceCenter, F: ReferenceFrame> Transform<cartesian::Direction<Geocentric, F>> for cartesian::Direction<C, F>
where
    C: NonGeocentric,
    F: MutableFrame,
    cartesian::Direction<Geocentric, F>: Transform<cartesian::Direction<Geocentric, Equatorial>>, // Required by Aberration
    cartesian::Direction<Geocentric, Equatorial>: Transform<cartesian::Direction<Geocentric, F>>,
{
    #[inline]
    fn transform(&self, jd: JulianDay) -> cartesian::Direction<Geocentric, F> {
        let geocentric = cartesian::Direction::<Geocentric, F>::from_vec3(
            self.as_vec3()
        );
        let equatorial: cartesian::Direction<Geocentric, Equatorial> = geocentric.transform(jd);
        let aberrated = apply_aberration_to_direction(
            cartesian::Direction::from_vec3(equatorial.as_vec3()), jd
        );
        aberrated.transform(jd)
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

impl<C1: ReferenceCenter, C2: ReferenceCenter> Transform<cartesian::Direction<C1, Equatorial>> for cartesian::Direction<C2, Equatorial>
where
    C1: NonGeocentric,
    C2: NonGeocentric
{
    #[inline]
    fn transform(&self, _jd: JulianDay) -> cartesian::Direction<C1, Equatorial> {
        cartesian::Direction::from_vec3(self.as_vec3())
    }
}

impl Transform<cartesian::Direction<Geocentric, Equatorial>>
    for cartesian::Direction<Geocentric, Equatorial>
{
    #[inline]
    fn transform(&self, _jd: JulianDay) -> cartesian::Direction<Geocentric, Equatorial> {
        cartesian::Direction::from_vec3(self.as_vec3())
    }
}

impl<C: ReferenceCenter, F: ReferenceFrame> Transform<spherical::Direction<Geocentric, F>> for spherical::Direction<C, F>
where
    C: NonGeocentric,
    F: MutableFrame,
    cartesian::Direction<Geocentric, F>: Transform<cartesian::Direction<Geocentric, Equatorial>>, // Required by Aberration
    cartesian::Direction<Geocentric, Equatorial>: Transform<cartesian::Direction<Geocentric, F>>,
{
    #[inline]
    fn transform(&self, jd: JulianDay) -> spherical::Direction<Geocentric, F> {
        let cart = self.to_cartesian();
        let cart_transformed = Transform::<cartesian::Direction<Geocentric, F>>::transform(&cart, jd);
        cart_transformed.to_spherical()
    }
}

impl<C: ReferenceCenter> Transform<spherical::Direction<C, Equatorial>> for spherical::Direction<Geocentric, Equatorial>
where
    C: NonGeocentric,
{
    #[inline]
    fn transform(&self, jd: JulianDay) -> spherical::Direction<C, Equatorial> {
        let cart = self.to_cartesian();
        let catr_trasnformed: cartesian::Direction<C, Equatorial> = cart.transform(jd);
        catr_trasnformed.to_spherical()
    }
}

impl<C1: ReferenceCenter, C2: ReferenceCenter> Transform<spherical::Direction<C1, Equatorial>> for spherical::Direction<C2, Equatorial>
where
    C1: NonGeocentric,
    C2: NonGeocentric
{
    #[inline]
    fn transform(&self, jd: JulianDay) -> spherical::Direction<C1, Equatorial> {
        let cart = self.to_cartesian();
        let catr_trasnformed: cartesian::Direction<C1, Equatorial> = cart.transform(jd);
        catr_trasnformed.to_spherical()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
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

    #[test]
    fn test_from_geocentric_to_heliocentric() {
        let dir = Direction::<Geocentric, Equatorial>::new(Degrees::new(10.0), Degrees::new(20.0));
        let jd = crate::units::JulianDay::J2000;
        let transformed: Direction<Heliocentric, Equatorial> = dir.transform(jd);
        // Should be close in direction, but not identical due to aberration, so just check type and not NaN
        assert!(!transformed.polar.as_f64().is_nan(), "Polar should not be NaN");
        assert!(!transformed.azimuth.as_f64().is_nan(), "Azimuth should not be NaN");
    }

/*
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
    }
    */
}
