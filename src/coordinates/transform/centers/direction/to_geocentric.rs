use crate::astro::aberration::apply_aberration_to_direction;
use crate::coordinates::{
    transform::Transform,
    frames::MutableFrame,
    centers::*,
    spherical,
    cartesian::direction::{Direction, Equatorial}
};
use crate::astro::JulianDate;

// ------------- If we transform TO a Geocentric Direction, we only need to apply aberration ------------------
impl<C, F> Transform<Direction<Geocentric, F>> for Direction<C, F>
where
    C: ReferenceCenter + NonGeocentric,
    F: MutableFrame,
    Direction<Geocentric, F>: Transform<Equatorial>, // Required by Aberration
    Equatorial: Transform<Direction<Geocentric, F>>,
{
    #[inline]
    fn transform(&self, jd: JulianDate) -> Direction<Geocentric, F> {
        // 1. Convert to Geocentric Equatorial coordinates
        let geocentric = Direction::<Geocentric, F>::from_vec3(self.as_vec3());
        // 2. Transform to Geocentric Equatorial
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

impl<C, F> Transform<spherical::Direction<Geocentric, F>> for spherical::Direction<C, F>
where
    C: ReferenceCenter + NonGeocentric,
    F: MutableFrame,
    Direction<Geocentric, F>: Transform<Equatorial>,
    Equatorial: Transform<Direction<Geocentric, F>>,
{
    #[inline]
    fn transform(&self, jd: JulianDate) -> spherical::Direction<Geocentric, F> {
        let cart = self.to_cartesian();
        let cart_transformed = Transform::<Direction<Geocentric, F>>::transform(&cart, jd);
        cart_transformed.to_spherical()
    }
}
