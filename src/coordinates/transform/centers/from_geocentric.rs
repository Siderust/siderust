use crate::units::{JulianDay, Unit};
use crate::coordinates::{
    frames::*, centers::*,
    cartesian, spherical,
};
use crate::coordinates::transform::Transform;
use crate::astro::aberration::remove_aberration_from_direction;


// ------------- If we transform FROM a Geocentric ReferenceCenter, we need to remove aberration ------------------
impl<C, F, U> Transform<cartesian::Direction<C, F, U>> for cartesian::Direction<Geocentric, F, U>
where
    C: ReferenceCenter + NonGeocentric,
    F: MutableFrame,
    U: Unit,
    cartesian::Direction<Geocentric, F, U>: Transform<cartesian::Direction<Geocentric, Equatorial, U>>, // Required by Aberration
    cartesian::Direction<Geocentric, Equatorial, U>: Transform<cartesian::Direction<Geocentric, F, U>>,
{
    #[inline]
    fn transform(&self, jd: JulianDay) -> cartesian::Direction<C, F, U> {
        // 1. Transform to Equatorial (already in Geocentric)
        let equatorial: cartesian::Direction<Geocentric, Equatorial, U> = self.transform(jd);
        // 2. Remove aberration
        let deaberrated =  remove_aberration_from_direction(
            equatorial, jd
        );
        // 3. Recover target Frame
        let target_center: cartesian::Direction::<Geocentric, F, U> = deaberrated.transform(jd);
        // 4. Transform target Center
        cartesian::Direction::<C, F, U>::from_vec3(target_center.as_vec3())
    }
}


// ------------- If we transform FROM a Geocentric ReferenceCenter, we need to remove aberration -----------------------
impl<C, F, U> Transform<spherical::Direction<C, F, U>> for spherical::Direction<Geocentric, F, U>
where
    C: ReferenceCenter + NonGeocentric,
    F: MutableFrame,
    U: Unit,
    cartesian::Direction<Geocentric, F, U>: Transform<cartesian::Direction<Geocentric, Equatorial, U>>,
    cartesian::Direction<Geocentric, Equatorial, U>: Transform<cartesian::Direction<Geocentric, F, U>>,
{
    #[inline]
    fn transform(&self, jd: JulianDay) -> spherical::Direction<C, F, U> {
        let cart = self.to_cartesian();
        let cart_transformed = Transform::<cartesian::Direction<C, F, U>>::transform(&cart, jd);
        cart_transformed.to_spherical()
    }
}
