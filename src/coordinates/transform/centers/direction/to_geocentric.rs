use crate::astro::aberration::apply_aberration_to_direction;
use crate::astro::JulianDate;
use crate::coordinates::spherical::direction::DirectionUnit;
use crate::coordinates::{
    cartesian::Vector,
    cartesian::direction::Equatorial,
    centers::*,
    frames::MutableFrame,
    transform::centers::TransformCenter,
    transform::TransformFrame,
};

// Heliocentric To Geocentric (Direction only - uses DirectionUnit)
impl<F: MutableFrame> TransformCenter<Vector<Geocentric, F, DirectionUnit>> for Vector<Heliocentric, F, DirectionUnit>
where
    Vector<Geocentric, F, DirectionUnit>: TransformFrame<Equatorial>,
    Equatorial: TransformFrame<Vector<Geocentric, F, DirectionUnit>>,
{
    fn to_center(&self, jd: JulianDate) -> Vector<Geocentric, F, DirectionUnit> {
        // 1. Convert to Geocentric coordinates
        let geocentric = Vector::<Geocentric, F, DirectionUnit>::from_vec3(self.as_vec3());
        // 2. Transform to Equatorial
        let equatorial: Equatorial<Geocentric> = geocentric.to_frame();
        // 3. Apply aberration
        let aberrated = apply_aberration_to_direction(
            Equatorial::<Geocentric>::from_vec3(equatorial.as_vec3()),
            jd,
        );
        // 4. Recover target Frame
        aberrated.to_frame()
    }
}

// Barycentric To Geocentric (Direction only - uses DirectionUnit)
impl<F: MutableFrame> TransformCenter<Vector<Geocentric, F, DirectionUnit>> for Vector<Barycentric, F, DirectionUnit>
where
    Vector<Geocentric, F, DirectionUnit>: TransformFrame<Equatorial>,
    Equatorial: TransformFrame<Vector<Geocentric, F, DirectionUnit>>,
{
    fn to_center(&self, jd: JulianDate) -> Vector<Geocentric, F, DirectionUnit> {
        // 1. Convert to Geocentric coordinates
        let geocentric = Vector::<Geocentric, F, DirectionUnit>::from_vec3(self.as_vec3());
        // 2. Transform to Equatorial
        let equatorial: Equatorial<Geocentric> = geocentric.to_frame();
        // 3. Apply aberration
        let aberrated = apply_aberration_to_direction(
            Equatorial::<Geocentric>::from_vec3(equatorial.as_vec3()),
            jd,
        );
        // 4. Recover target Frame
        aberrated.to_frame()
    }
}
