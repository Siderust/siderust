pub mod to_barycentric;
pub mod to_heliocentric;
pub mod to_geocentric;

use crate::astro::JulianDate;
use crate::coordinates::{
    cartesian,
    frames::ReferenceFrame,
    centers::ReferenceCenter,
    transform::Transform,
    transform::TransformCenter,
};

impl<C1, C2, F> Transform<cartesian::Direction<C1, F>> for cartesian::Direction<C2, F>
where
    C1: ReferenceCenter,
    C2: ReferenceCenter,
    F: ReferenceFrame,
    cartesian::Direction<C2, F>: TransformCenter<cartesian::Direction<C1, F>>
{
    fn transform(&self, jd: JulianDate) -> cartesian::Direction<C1, F> {
        self.to_center(jd)
    }
}
