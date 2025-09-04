use crate::astro::aberration::remove_aberration_from_direction;
use crate::astro::JulianDate;
use crate::coordinates::transform::centers::TransformCenter;
use crate::coordinates::{
    cartesian::direction::{Direction, Equatorial},
    centers::*,
    frames::MutableFrame,
    transform::TransformFrame,
};

// Barycentric To Heliocentric
impl<F: MutableFrame> TransformCenter<Direction<Heliocentric, F>> for Direction<Barycentric, F> {
    fn to_center(&self, _jd: JulianDate) -> Direction<Heliocentric, F> {
        Direction::from_vec3(self.as_vec3())
    }
}

// Geocentric To Heliocentric
impl<F: MutableFrame> TransformCenter<Direction<Heliocentric, F>> for Direction<Geocentric, F>
where
    Direction<Geocentric, F>: TransformFrame<Equatorial>, // ToEquatorial
    Equatorial: TransformFrame<Direction<Geocentric, F>>, // FromEquatorial
{
    #[inline]
    fn to_center(&self, jd: JulianDate) -> Direction<Heliocentric, F> {
        // 1. Transform to Equatorial (already in Geocentric)
        let equatorial: Equatorial = self.to_frame();
        // 2. Remove aberration
        let deaberrated = remove_aberration_from_direction(equatorial, jd);
        // 3. Recover target Frame
        let target_center: Direction<Geocentric, F> = deaberrated.to_frame();
        // 4. Transform target Center
        Direction::<Heliocentric, F>::from_vec3(target_center.as_vec3())
    }
}
