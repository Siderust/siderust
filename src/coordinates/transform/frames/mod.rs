pub mod to_ecliptic;
pub mod to_equatorial;
pub mod to_horizontal;
pub mod to_icrs;

use crate::astro::JulianDate;
use crate::coordinates::cartesian::Vector;
use crate::coordinates::centers::ReferenceCenter;
use crate::coordinates::frames::MutableFrame;
use crate::coordinates::spherical::SphericalCoord;
use crate::coordinates::{cartesian, spherical};
use qtty::Unit;

use crate::coordinates::transform::Transform;

pub trait TransformFrame<Coord> {
    fn to_frame(&self) -> Coord;
}

// Implement Identity frame transform
impl<C, F, U> TransformFrame<Vector<C, F, U>> for Vector<C, F, U>
where
    U: Unit,
    F: MutableFrame,
    C: ReferenceCenter,
{
    fn to_frame(&self) -> Vector<C, F, U> {
        Vector::from_vec3(self.center_params().clone(), self.as_vec3())
    }
}

impl<C, F1, F2, U> TransformFrame<spherical::SphericalCoord<C, F2, U>>
    for spherical::SphericalCoord<C, F1, U>
where
    cartesian::Vector<C, F1, U>: TransformFrame<cartesian::Vector<C, F2, U>>,
    C: ReferenceCenter,
    F1: MutableFrame,
    F2: MutableFrame,
    U: Unit,
{
    fn to_frame(&self) -> spherical::SphericalCoord<C, F2, U> {
        self.to_cartesian().to_frame().to_spherical()
    }
}

impl<C, F, U> SphericalCoord<C, F, U>
where
    C: ReferenceCenter,
    F: MutableFrame,
    U: qtty::Unit,
{
    pub fn to_frame<F2: MutableFrame>(&self) -> SphericalCoord<C, F2, U>
    where
        Vector<C, F, U>: Transform<Vector<C, F2, U>>,
        Vector<C, F, U>: for<'a> From<&'a SphericalCoord<C, F, U>>, // to_cartesian
        SphericalCoord<C, F2, U>: for<'a> From<&'a Vector<C, F2, U>>, // to_spherical
    {
        self.to_cartesian()
            .transform(JulianDate::J2000)
            .to_spherical()
    }
}

/*
#[cfg(test)]
mod tests {
    use crate::coordinates::centers;
    use crate::coordinates::frames;
    use crate::coordinates::spherical::direction::Direction;
    use qtty::DEG;

    #[test]
    fn test_to_center() {
        let ecl  =  Direction::<centers::Heliocentric, frames::Ecliptic>::new(90.0*DEG, 45.0*DEG);
        let equ: Direction<centers::Heliocentric, frames::Equatorial> = ecl.to_frame::<frames::Equatorial>();
        let expected: Direction<centers::Heliocentric, frames::Equatorial> = (&ecl).into();

        assert_eq!(equ.polar, expected.polar);
        assert_eq!(equ.azimuth, expected.azimuth);
        assert_eq!(equ.distance, expected.distance);
    }
}
*/
