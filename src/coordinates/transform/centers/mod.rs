pub mod to_barycentric;
pub mod to_heliocentric;
pub mod to_geocentric;

use crate::units::JulianDay;
use crate::coordinates::{
    frames, centers, cartesian,
};
use crate::coordinates::transform::Transform;

impl<C1, C2, F> Transform<cartesian::Direction<C1, F>> for cartesian::Direction<C2, F>
where 
    C1: centers::ReferenceCenter,
    C2: centers::ReferenceCenter,
    F: frames::ReferenceFrame,
{
    fn transform(&self, _jd: JulianDay) -> cartesian::Direction<C1, F> {
        // Direction does not change with respect to the center, so we can just copy the vector
        cartesian::Direction::from_vec3(self.as_vec3())
    }
}
