pub mod to_barycentric;
pub mod to_heliocentric;
pub mod to_geocentric;
pub mod from_geocentric;

use crate::units::JulianDay;
use crate::coordinates::{
    frames::*, centers::*,
    cartesian, spherical,
};
use crate::coordinates::transform::Transform;
use crate::units::Unit;

// ------------- Identity Transform -------------
impl<U: Unit> Transform<cartesian::Direction<Geocentric, Equatorial, U>>
    for cartesian::Direction<Geocentric, Equatorial, U>
{
    #[inline]
    fn transform(&self, _jd: JulianDay) -> cartesian::Direction<Geocentric, Equatorial, U> {
        cartesian::Direction::from_vec3(self.as_vec3())
    }
}

// ------------- If None of the centers are geocentric, we can just pass the spherical coordinates through -------------
impl<C1, C2, F, U> Transform<cartesian::Direction<C1, F, U>> for cartesian::Direction<C2, F, U>
where
    C1: ReferenceCenter + NonGeocentric,
    C2: ReferenceCenter + NonGeocentric,
    F: MutableFrame,
    U: Unit,
{
    #[inline]
    fn transform(&self, _jd: JulianDay) -> cartesian::Direction<C1, F, U> {
        cartesian::Direction::from_vec3(self.as_vec3())
    }
}

impl<C1, C2, F, U> Transform<spherical::Direction<C1, F, U>> for spherical::Direction<C2, F, U>
where
    C1: ReferenceCenter + NonGeocentric,
    C2: ReferenceCenter + NonGeocentric,
    F: MutableFrame,
    U: Unit,
{
    #[inline]
    fn transform(&self, _jd: JulianDay) -> spherical::Direction<C1, F, U> {
        spherical::Direction::<C1, F, U>::new_spherical_coord(
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
        assert_eq!(dir.polar.as_f64(), transformed.polar.as_f64(), "Polar should not change");
        assert_eq!(dir.azimuth.as_f64(), transformed.azimuth.as_f64(), "Azimuth should not change");
    }

    #[test]
    fn test_heliocentric_to_barycentric() {
        let dir = Direction::<Heliocentric, ICRS>::new(Degrees::new(100.0), Degrees::new(-10.0));
        let jd = crate::units::JulianDay::J2000;
        let transformed: Direction<Barycentric, ICRS> = dir.transform(jd);
        assert_eq!(dir.polar.as_f64(), transformed.polar.as_f64(), "Polar should not change");
        assert_eq!(dir.azimuth.as_f64(), transformed.azimuth.as_f64(), "Azimuth should not change");
    }
}
