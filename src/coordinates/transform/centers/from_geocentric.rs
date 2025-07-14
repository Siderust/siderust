use crate::astro::JulianDate;
use crate::coordinates::{
    centers::*,
    frames::MutableFrame,
    spherical, cartesian,
    cartesian::direction::Equatorial
};
use crate::coordinates::transform::Transform;
use crate::astro::aberration::remove_aberration_from_direction;

// ------------- If we transform FROM a Geocentric ReferenceCenter, we need to remove aberration ------------------
impl<C, F> Transform<cartesian::Direction<C, F>> for cartesian::Direction<Geocentric, F>
where
    C: ReferenceCenter + NonGeocentric,
    F: MutableFrame,
    cartesian::Direction<Geocentric, F>: Transform<Equatorial>, // ToEquatorial
    Equatorial: Transform<cartesian::Direction<Geocentric, F>>, // FromEquatorial
{
    #[inline]
    fn transform(&self, jd: JulianDate) -> cartesian::Direction<C, F> {
        // 1. Transform to Equatorial (already in Geocentric)
        let equatorial: Equatorial = self.transform(jd);
        // 2. Remove aberration
        let deaberrated =  remove_aberration_from_direction(
            equatorial, jd
        );
        // 3. Recover target Frame
        let target_center: cartesian::Direction::<Geocentric, F> = deaberrated.transform(jd);
        // 4. Transform target Center
        cartesian::Direction::<C, F>::from_vec3(target_center.as_vec3())
    }
}


// ------------- If we transform FROM a Geocentric ReferenceCenter, we need to remove aberration -----------------------
impl<C, F> Transform<spherical::Direction<C, F>> for spherical::Direction<Geocentric, F>
where
    C: ReferenceCenter + NonGeocentric,
    F: MutableFrame,
    cartesian::Direction<Geocentric, F>: Transform<Equatorial>, // ToEquatorial
    Equatorial: Transform<cartesian::Direction<Geocentric, F>>, // FromEquatorial
{
    #[inline]
    fn transform(&self, jd: JulianDate) -> spherical::Direction<C, F> {
        let cart = self.to_cartesian();
        let cart_transformed = Transform::<cartesian::Direction<C, F>>::transform(&cart, jd);
        cart_transformed.to_spherical()
    }
}
