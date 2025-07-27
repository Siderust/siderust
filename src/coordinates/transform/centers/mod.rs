pub mod position;
pub mod direction;
pub mod from_geocentric;

use crate::astro::JulianDate;
use crate::coordinates::{
    frames, centers::*,
    cartesian, spherical,
};
use crate::coordinates::transform::Transform;

// ------------- Identity Transform -------------
impl Transform<cartesian::direction::Equatorial>
    for cartesian::direction::Equatorial
{
    #[inline]
    fn transform(&self, _jd: JulianDate) -> cartesian::direction::Equatorial {
        cartesian::Direction::from_vec3(self.as_vec3())
    }
}

// ------------- If None of the centers are geocentric, we can just pass the spherical coordinates through -------------
impl<C1, C2, F> Transform<cartesian::Direction<C1, F>> for cartesian::Direction<C2, F>
where
    C1: ReferenceCenter + NonGeocentric,
    C2: ReferenceCenter + NonGeocentric,
    F: frames::MutableFrame,
{
    #[inline]
    fn transform(&self, _jd: JulianDate) -> cartesian::Direction<C1, F> {
        cartesian::Direction::from_vec3(self.as_vec3())
    }
}

impl<C1, C2, F> Transform<spherical::Direction<C1, F>> for spherical::Direction<C2, F>
where
    C1: ReferenceCenter + NonGeocentric,
    C2: ReferenceCenter + NonGeocentric,
    F: frames::MutableFrame,
{
    #[inline]
    fn transform(&self, _jd: JulianDate) -> spherical::Direction<C1, F> {
        spherical::Direction::new_raw(
            self.polar,
            self.azimuth,
            self.distance
        )
    }
}


impl<C, F, U> cartesian::Vector<C, F, U>
where
    C: ReferenceCenter,
    F: frames::ReferenceFrame,
    U: crate::units::Unit,
{
    pub fn to_center<C2: ReferenceCenter>(&self, jd: JulianDate) -> cartesian::Vector<C2, F, U>
    where
        cartesian::Vector<C, F, U>: Transform<cartesian::Vector<C2, F, U>>,
    {
        self.transform(jd)
    }
}

impl<C, F, U> spherical::SphericalCoord<C, F, U>
where
    C: ReferenceCenter,
    F: frames::ReferenceFrame,
    U: crate::units::Unit,
{
    pub fn to_center<C2: ReferenceCenter>(&self, jd: JulianDate) -> spherical::SphericalCoord<C2, F, U>
    where
        spherical::SphericalCoord<C, F, U>: Transform<spherical::SphericalCoord<C2, F, U>>,
    {
        self.transform(jd)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coordinates::spherical::direction::*;
    use crate::coordinates::frames;
    use crate::units::Degrees;

    #[test]
    fn test_from_heliocentric_to_geocentric() {
        let dir = Direction::<Heliocentric, frames::Equatorial>::new(Degrees::new(120.0), Degrees::new(-30.0));
        let jd = crate::astro::JulianDate::J2000;
        let transformed: Equatorial = dir.transform(jd);
        assert!(!transformed.polar.value().is_nan(), "Polar should not be NaN");
        assert!(!transformed.azimuth.value().is_nan(), "Azimuth should not be NaN");
    }

    #[test]
    fn test_from_barycentric_to_geocentric() {
        let dir = Direction::<Barycentric, frames::Equatorial>::new(Degrees::new(120.0), Degrees::new(-30.0));
        let jd = crate::astro::JulianDate::J2000;
        let transformed: Equatorial = dir.transform(jd);
        assert!(!transformed.polar.value().is_nan(), "Polar should not be NaN");
        assert!(!transformed.azimuth.value().is_nan(), "Azimuth should not be NaN");
    }

    #[test]
    fn test_from_geocentric_to_heliocentric() {
        let dir = Equatorial::new(Degrees::new(10.0), Degrees::new(20.0));
        let jd = crate::astro::JulianDate::J2000;
        let transformed: Direction<Heliocentric, frames::Equatorial> = dir.transform(jd);
        // Should be close in direction, but not identical due to aberration, so just check type and not NaN
        assert!(!transformed.polar.value().is_nan(), "Polar should not be NaN");
        assert!(!transformed.azimuth.value().is_nan(), "Azimuth should not be NaN");
    }

    #[test]
    fn test_barycentric_to_heliocentric() {
        let dir = ICRS::new(Degrees::new(200.0), Degrees::new(45.0));
        let jd = crate::astro::JulianDate::J2000;
        let transformed: HCRS = dir.transform(jd);
        assert_eq!(dir.polar.value(), transformed.polar.value(), "Polar should not change");
        assert_eq!(dir.azimuth.value(), transformed.azimuth.value(), "Azimuth should not change");
    }

    #[test]
    fn test_heliocentric_to_barycentric() {
        let dir = HCRS::new(Degrees::new(100.0), Degrees::new(-10.0));
        let jd = crate::astro::JulianDate::J2000;
        let transformed: ICRS = dir.transform(jd);
        assert_eq!(dir.polar.value(), transformed.polar.value(), "Polar should not change");
        assert_eq!(dir.azimuth.value(), transformed.azimuth.value(), "Azimuth should not change");
    }

    #[test]
    fn test_to_center() {
        let icrs = ICRS::new(Degrees::new(100.0), Degrees::new(-10.0));

        let gcrs: GCRS = icrs.to_center::<Geocentric>(JulianDate::J2000);
        let expected: GCRS = icrs.transform(JulianDate::J2000);
        assert_eq!(gcrs.polar, expected.polar);
        assert_eq!(gcrs.azimuth, expected.azimuth);
        assert_eq!(gcrs.distance, expected.distance);
    }
}
