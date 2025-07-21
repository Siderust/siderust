pub mod to_ecliptic;
pub mod to_icrs;
pub mod to_equatorial;

use crate::coordinates::cartesian;
use crate::coordinates::spherical::SphericalCoord;
use crate::coordinates::frames::MutableFrame;
use crate::coordinates::centers::ReferenceCenter;
use crate::coordinates::transform::Transform;
use crate::astro::JulianDate;

impl<C, F, U> cartesian::Vector<C, F, U>
where
    C: ReferenceCenter,
    F: MutableFrame,
    U: crate::units::Unit,
{
    pub fn to_frame<F2: MutableFrame>(&self) -> cartesian::Vector<C, F2, U>
    where
        cartesian::Vector<C, F, U>: Transform<cartesian::Vector<C, F2, U>>,
    {
        self.transform(JulianDate::J2000)
    }
}

impl<C, F, U> SphericalCoord<C, F, U>
where
    C: ReferenceCenter,
    F: MutableFrame,
    U: crate::units::Unit,
{
    pub fn to_frame<F2: MutableFrame>(&self) -> SphericalCoord<C, F2, U>
    where
        cartesian::Vector<C, F, U>: Transform<cartesian::Vector<C, F2, U>>,
        cartesian::Vector<C, F, U>: for<'a> From<&'a SphericalCoord<C, F, U>>, // to_cartesian
        SphericalCoord<C, F2, U>: for<'a> From<&'a cartesian::Vector<C, F2, U>>, // to_spherical
    {
        self.to_cartesian()
            .transform(JulianDate::J2000)
            .to_spherical()
    }
}


#[cfg(test)]
mod tests {
    use crate::coordinates::centers;
    use crate::coordinates::frames;
    use crate::coordinates::spherical::direction::Direction;
    use crate::units::DEG;

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