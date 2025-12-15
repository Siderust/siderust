use crate::astro::aberration::remove_aberration_from_direction;
use crate::astro::JulianDate;
use crate::coordinates::spherical::direction::DirectionUnit;
use crate::coordinates::transform::centers::TransformCenter;
use crate::coordinates::{
    cartesian::direction::Equatorial, cartesian::Vector, centers::*, frames::MutableFrame,
    transform::TransformFrame,
};

// Barycentric To Heliocentric (Direction only - uses DirectionUnit)
impl<F: MutableFrame> TransformCenter<Vector<Heliocentric, F, DirectionUnit>>
    for Vector<Barycentric, F, DirectionUnit>
{
    fn to_center(&self, _jd: JulianDate) -> Vector<Heliocentric, F, DirectionUnit> {
        Vector::from_vec3_origin(self.as_vec3())
    }
}

// Geocentric To Heliocentric (Direction only - uses DirectionUnit)
impl<F: MutableFrame> TransformCenter<Vector<Heliocentric, F, DirectionUnit>>
    for Vector<Geocentric, F, DirectionUnit>
where
    Vector<Geocentric, F, DirectionUnit>: TransformFrame<Equatorial>, // ToEquatorial
    Equatorial: TransformFrame<Vector<Geocentric, F, DirectionUnit>>, // FromEquatorial
{
    #[inline]
    fn to_center(&self, jd: JulianDate) -> Vector<Heliocentric, F, DirectionUnit> {
        // 1. Transform to Equatorial (already in Geocentric)
        let equatorial: Equatorial = self.to_frame();
        // 2. Remove aberration
        let deaberrated = remove_aberration_from_direction(equatorial, jd);
        // 3. Recover target Frame
        let target_center: Vector<Geocentric, F, DirectionUnit> = deaberrated.to_frame();
        // 4. Transform target Center
        Vector::<Heliocentric, F, DirectionUnit>::from_vec3_origin(target_center.as_vec3())
    }
}
