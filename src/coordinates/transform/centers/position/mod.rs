pub mod to_barycentric;
pub mod to_heliocentric;
pub mod to_geocentric;

use crate::astro::JulianDate;
use crate::units::LengthUnit;
use crate::coordinates::{
    cartesian,
    frames::ReferenceFrame,
    centers::ReferenceCenter,
    transform::Transform,
    transform::TransformCenter,
};


impl<C1, C2, F, U> Transform<cartesian::Position<C1, F, U>> for cartesian::Position<C2, F, U>
where
    C1: ReferenceCenter,
    C2: ReferenceCenter,
    F: ReferenceFrame,
    U: LengthUnit,
    cartesian::Position<C2, F, U>: TransformCenter<cartesian::Position<C1, F, U>>
{
    fn transform(&self, jd: JulianDate) -> cartesian::Position<C1, F, U> {
        self.to_center(jd)
    }
}
