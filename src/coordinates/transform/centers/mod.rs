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

// ------------- If we transform TO a Geocentric ReferenceCenter, we need to apply aberration ------------------
impl<C: ReferenceCenter, F: ReferenceFrame> Transform<cartesian::Direction<Geocentric, F>> for cartesian::Direction<C, F>
where
    C: NonGeocentric,
    F: MutableFrame,
    cartesian::Direction<Geocentric, F>: Transform<cartesian::Direction<Geocentric, Equatorial>>, // Required by Aberration
    cartesian::Direction<Geocentric, Equatorial>: Transform<cartesian::Direction<Geocentric, F>>,
{
    #[inline]
    fn transform(&self, jd: JulianDay) -> cartesian::Direction<Geocentric, F> {
        // 1. Convert to Geocentric Equatorial coordinates
        let geocentric = cartesian::Direction::<Geocentric, F>::from_vec3(
            self.as_vec3()
        );
        // 2. Transform to Geocentric Equatorial
        let equatorial: cartesian::Direction<Geocentric, Equatorial> = geocentric.transform(jd);
        // 3. Apply aberration
        let aberrated = apply_aberration_to_direction(
            cartesian::Direction::from_vec3(equatorial.as_vec3()), jd
        );
        // 4. Recover target Frame
        Transform::<cartesian::Direction<Geocentric, F>>::transform(&aberrated, jd)
    }
}

// ------------- If we transform FROM a Geocentric ReferenceCenter, we need to remove aberration ------------------
impl<C, F> Transform<cartesian::Direction<C, F>> for cartesian::Direction<Geocentric, F>
where
    C: ReferenceCenter + NonGeocentric,
    F: MutableFrame,
    cartesian::Direction<Geocentric, F>: Transform<cartesian::Direction<Geocentric, Equatorial>>, // Required by Aberration
    cartesian::Direction<Geocentric, Equatorial>: Transform<cartesian::Direction<Geocentric, F>>,
{
    #[inline]
    fn transform(&self, jd: JulianDay) -> cartesian::Direction<C, F> {
        // 1. Transform to Equatorial (already in Geocentric)
        let equatorial: cartesian::Direction<Geocentric, Equatorial> = self.transform(jd);
        // 2. Remove aberration
        let deaberrated =  remove_aberration_from_direction(
            equatorial, jd
        );
        // 3. Recover target Frame
        let target_center: cartesian::Direction::<Geocentric, F> = deaberrated.transform(jd);
        // 4. Transform target Center
        cartesian::Direction::<C, F>::from_vec3(target_center.as_vec3())
    }
}

// ------------- If None of the centers are geocentric, we can just pass the spherical coordinates through -------------
impl<C1, C2, F> Transform<cartesian::Direction<C1, F>> for cartesian::Direction<C2, F>
where
    C1: ReferenceCenter + NonGeocentric,
    C2: ReferenceCenter + NonGeocentric,
    F: MutableFrame,
{
    #[inline]
    fn transform(&self, _jd: JulianDay) -> cartesian::Direction<C1, F> {
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

// ------------- If we transform TO a Geocentric ReferenceCenter, we need to apply aberration -------------------------
impl<C: ReferenceCenter, F: ReferenceFrame> Transform<spherical::Direction<Geocentric, F>> for spherical::Direction<C, F>
where
    C: NonGeocentric,
    F: MutableFrame,
    cartesian::Direction<Geocentric, F>: Transform<cartesian::Direction<Geocentric, Equatorial>>,
    cartesian::Direction<Geocentric, Equatorial>: Transform<cartesian::Direction<Geocentric, F>>,
{
    #[inline]
    fn transform(&self, jd: JulianDay) -> spherical::Direction<Geocentric, F> {
        let cart = self.to_cartesian();
        let cart_transformed = Transform::<cartesian::Direction<Geocentric, F>>::transform(&cart, jd);
        cart_transformed.to_spherical()
    }
}

// ------------- If we transform FROM a Geocentric ReferenceCenter, we need to remove aberration -----------------------
impl<C, F> Transform<spherical::Direction<C, F>> for spherical::Direction<Geocentric, F>
where
    C: ReferenceCenter + NonGeocentric,
    F: MutableFrame,
    cartesian::Direction<Geocentric, F>: Transform<cartesian::Direction<Geocentric, Equatorial>>,
    cartesian::Direction<Geocentric, Equatorial>: Transform<cartesian::Direction<Geocentric, F>>,
{
    #[inline]
    fn transform(&self, jd: JulianDay) -> spherical::Direction<C, F> {
        let cart = self.to_cartesian();
        let cart_transformed = Transform::<cartesian::Direction<C, F>>::transform(&cart, jd);
        cart_transformed.to_spherical()
    }
}


// ------------- If None of the centers are geocentric, we can just pass the spherical coordinates through -------------
impl<C1, C2, F> Transform<spherical::Direction<C1, F>> for spherical::Direction<C2, F>
where
    C1: ReferenceCenter + NonGeocentric,
    C2: ReferenceCenter + NonGeocentric,
    F: MutableFrame,
{
    #[inline]
    fn transform(&self, _jd: JulianDay) -> spherical::Direction<C1, F> {
        spherical::Direction::<C1, F>::new_spherical_coord(
            self.polar,
            self.azimuth,
            None
        )
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

}
