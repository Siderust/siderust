use crate::astro::JulianDate;
use crate::coordinates::{
    centers::*,
    frames::MutableFrame,
    cartesian::direction::{Direction, Equatorial}
};
use crate::coordinates::transform::Transform;
use crate::coordinates::transform::centers::TransformCenter;
use crate::astro::aberration::remove_aberration_from_direction;

// Barycentric To Heliocentric
impl<F: MutableFrame> TransformCenter<Direction<Heliocentric, F>>
    for Direction<Barycentric, F>
where
{
    fn to_center(&self, _jd: JulianDate) -> Direction<Heliocentric, F> {
        Direction::from_vec3(self.as_vec3())
    }
}

// Geocentric To Heliocentric
impl<F: MutableFrame> TransformCenter<Direction<Heliocentric, F>>
    for Direction<Geocentric, F>
where
    Direction<Geocentric, F>: Transform<Equatorial>, // ToEquatorial
    Equatorial: Transform<Direction<Geocentric, F>>, // FromEquatorial
{
    #[inline]
    fn to_center(&self, jd: JulianDate) -> Direction<Heliocentric, F> {
        // 1. Transform to Equatorial (already in Geocentric)
        let equatorial: Equatorial = self.transform(jd);
        // 2. Remove aberration
        let deaberrated =  remove_aberration_from_direction(
            equatorial, jd
        );
        // 3. Recover target Frame
        let target_center: Direction::<Geocentric, F> = deaberrated.transform(jd);
        // 4. Transform target Center
        Direction::<Heliocentric, F>::from_vec3(target_center.as_vec3())
    }
}

//TODO: REMOVE ME

// Barycentric To Heliocentric
impl<F: MutableFrame> Transform<Direction<Heliocentric, F>>
    for Direction<Barycentric, F>
where
{
    fn transform(&self, _jd: JulianDate) -> Direction<Heliocentric, F> {
        self.to_center(_jd)
    }
}

impl<F: MutableFrame> Transform<Direction<Heliocentric, F>>
    for Direction<Geocentric, F>
where
    Direction<Geocentric, F>: Transform<Equatorial>, // ToEquatorial
    Equatorial: Transform<Direction<Geocentric, F>>, // FromEquatorial
{
    #[inline]
    fn transform(&self, jd: JulianDate) -> Direction<Heliocentric, F> {
        self.to_center(jd)
    }
}
