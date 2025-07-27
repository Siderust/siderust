use crate::astro::aberration::apply_aberration_to_direction;
use crate::coordinates::{
    transform::Transform,
    transform::centers::TransformCenter,
    frames::MutableFrame,
    centers::*,
    cartesian::direction::{Direction, Equatorial}
};
use crate::astro::JulianDate;

// Heliocentric To Geocentric
impl<F: MutableFrame> TransformCenter<Direction<Geocentric, F>>
    for Direction<Heliocentric, F>
where
    Direction<Geocentric, F>: Transform<Equatorial>,
    Equatorial: Transform<Direction<Geocentric, F>>,
{
    fn to_center(&self, jd: JulianDate) -> Direction<Geocentric, F> {
        // 1. Convert to Geocentric coordinates
        let geocentric = Direction::<Geocentric, F>::from_vec3(self.as_vec3());
        // 2. Transform to Equatorial
        let equatorial: Equatorial<Geocentric> = geocentric.transform(jd);
        // 3. Apply aberration
        let aberrated = apply_aberration_to_direction(
            Equatorial::<Geocentric>::from_vec3(equatorial.as_vec3()),
            jd,
        );
        // 4. Recover target Frame
        Transform::<Direction<Geocentric, F>>::transform(&aberrated, jd)
    }
}

// Barycentric To Geocentric
impl<F: MutableFrame> TransformCenter<Direction<Geocentric, F>>
    for Direction<Barycentric, F>
where
    Direction<Geocentric, F>: Transform<Equatorial>,
    Equatorial: Transform<Direction<Geocentric, F>>,
{
    fn to_center(&self, jd: JulianDate) -> Direction<Geocentric, F> {
        // 1. Convert to Geocentric coordinates
        let geocentric = Direction::<Geocentric, F>::from_vec3(self.as_vec3());
        // 2. Transform to Equatorial
        let equatorial: Equatorial<Geocentric> = geocentric.transform(jd);
        // 3. Apply aberration
        let aberrated = apply_aberration_to_direction(
            Equatorial::<Geocentric>::from_vec3(equatorial.as_vec3()),
            jd,
        );
        // 4. Recover target Frame
        Transform::<Direction<Geocentric, F>>::transform(&aberrated, jd)
    }
}
