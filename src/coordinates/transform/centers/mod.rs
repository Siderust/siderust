pub mod direction;
pub mod position;

use crate::astro::JulianDate;
use crate::coordinates::{cartesian::Vector, centers::*, frames};
use crate::units::Unit;

pub trait TransformCenter<Coord> {
    fn to_center(&self, jd: crate::astro::JulianDate) -> Coord;
}

impl<C, F, U> TransformCenter<Vector<C, F, U>> for Vector<C, F, U>
where
    C: ReferenceCenter,
    F: frames::ReferenceFrame,
    U: Unit,
{
    fn to_center(&self, _jd: JulianDate) -> Vector<C, F, U> {
        Vector::<C, F, U>::from_vec3(self.as_vec3())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coordinates::frames;
    use crate::coordinates::spherical::direction::*;
    use crate::coordinates::transform::Transform;
    use crate::units::Degrees;

    #[test]
    fn test_from_heliocentric_to_geocentric() {
        let dir = Direction::<Heliocentric, frames::Equatorial>::new(
            Degrees::new(120.0),
            Degrees::new(-30.0),
        );
        let jd = crate::astro::JulianDate::J2000;
        let transformed: Equatorial = dir.transform(jd);
        assert!(
            !transformed.polar.value().is_nan(),
            "Polar should not be NaN"
        );
        assert!(
            !transformed.azimuth.value().is_nan(),
            "Azimuth should not be NaN"
        );
    }

    #[test]
    fn test_from_barycentric_to_geocentric() {
        let dir = Direction::<Barycentric, frames::Equatorial>::new(
            Degrees::new(120.0),
            Degrees::new(-30.0),
        );
        let jd = crate::astro::JulianDate::J2000;
        let transformed: Equatorial = dir.transform(jd);
        assert!(
            !transformed.polar.value().is_nan(),
            "Polar should not be NaN"
        );
        assert!(
            !transformed.azimuth.value().is_nan(),
            "Azimuth should not be NaN"
        );
    }

    #[test]
    fn test_from_geocentric_to_heliocentric() {
        let dir = Equatorial::new(Degrees::new(10.0), Degrees::new(20.0));
        let jd = crate::astro::JulianDate::J2000;
        let transformed: Direction<Heliocentric, frames::Equatorial> = dir.transform(jd);
        // Should be close in direction, but not identical due to aberration, so just check type and not NaN
        assert!(
            !transformed.polar.value().is_nan(),
            "Polar should not be NaN"
        );
        assert!(
            !transformed.azimuth.value().is_nan(),
            "Azimuth should not be NaN"
        );
    }

    #[test]
    fn test_heliocentric_to_barycentric() {
        let dir = HCRS::new(Degrees::new(100.0), Degrees::new(-10.0));
        let jd = crate::astro::JulianDate::J2000;
        let transformed: ICRS = dir.transform(jd);
        assert_eq!(
            dir.polar.value(),
            transformed.polar.value(),
            "Polar should not change"
        );
        assert_eq!(
            dir.azimuth.value(),
            transformed.azimuth.value(),
            "Azimuth should not change"
        );
    }
}
