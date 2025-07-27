use crate::astro::JulianDate;
use crate::coordinates::{
    centers::*,
    frames::MutableFrame,
    transform::TransformFrame,
    cartesian::direction::{Direction, Equatorial},
};
use crate::coordinates::transform::Transform;
use crate::coordinates::transform::centers::TransformCenter;
use crate::astro::aberration::remove_aberration_from_direction;

// Heliocentric To Barycentric
impl<F: MutableFrame> TransformCenter<Direction<Barycentric, F>>
    for Direction<Heliocentric, F>
where
{
    fn to_center(&self, _jd: JulianDate) -> Direction<Barycentric, F> {
        Direction::from_vec3(self.as_vec3())
    }
}

// Geocentric To Barycentric
impl<F: MutableFrame> TransformCenter<Direction<Barycentric, F>>
    for Direction<Geocentric, F>
where
    Direction<Geocentric, F>: TransformFrame<Equatorial>, // ToEquatorial
    Equatorial: TransformFrame<Direction<Geocentric, F>>, // FromEquatorial
{
    #[inline]
    fn to_center(&self, jd: JulianDate) -> Direction<Barycentric, F> {
        // 1. Transform to Equatorial (already in Geocentric)
        let equatorial: Equatorial = self.to_frame();
        // 2. Remove aberration
        let deaberrated =  remove_aberration_from_direction(equatorial, jd);
        // 3. Recover target Frame
        let target_center: Direction::<Geocentric, F> = deaberrated.to_frame();
        // 4. Transform target Center
        Direction::<Barycentric, F>::from_vec3(target_center.as_vec3())
    }
}


//TODO: REMOVE ME
impl<F: MutableFrame> Transform<Direction<Barycentric, F>>
    for Direction<Heliocentric, F>
where
{
    fn transform(&self, _jd: JulianDate) -> Direction<Barycentric, F> {
        self.to_center(_jd)
    }
}

impl<F: MutableFrame> Transform<Direction<Barycentric, F>>
    for Direction<Geocentric, F>
where
    Direction<Geocentric, F>: TransformFrame<Equatorial>, // ToEquatorial
    Equatorial: TransformFrame<Direction<Geocentric, F>>, // FromEquatorial
{
    #[inline]
    fn transform(&self, jd: JulianDate) -> Direction<Barycentric, F> {
        self.to_center(jd)
    }
}
