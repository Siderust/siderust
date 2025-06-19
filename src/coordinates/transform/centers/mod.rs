pub mod to_barycentric;
pub mod to_heliocentric;
pub mod to_geocentric;

use crate::units::JulianDay;
use crate::coordinates::{
    frames::*, centers::*,
    cartesian::Direction,
};
use crate::coordinates::frames;
use crate::coordinates::transform::Transform;
use crate::astro::aberration::{
    apply_aberration_to_direction,
    remove_aberration_from_direction,
};

#[inline]
#[must_use]
fn aberrate<F, Op>(
    dir: &Direction<Geocentric, F>,
    jd: JulianDay,
    op: Op,
) -> Direction<Geocentric, F>
where
    F: frames::ReferenceFrame,
    Op: FnOnce(Direction<Geocentric, Equatorial>, JulianDay)
        -> Direction<Geocentric, Equatorial>,
    Direction<Geocentric, Equatorial>: for<'a> From<&'a Direction<Geocentric, F>>,
    Direction<Geocentric, F>: for<'a> From<&'a Direction<Geocentric, Equatorial>>,
{
    let geo_eq = Direction::<Geocentric, Equatorial>::from(dir);
    let geo_eq = op(geo_eq, jd);
    Direction::<Geocentric, F>::from(&geo_eq)
}

impl<C1, C2, F> Transform<Direction<C1, F>> for Direction<C2, F>
where
    C1: ReferenceCenter,
    C2: ReferenceCenter,
    F: frames::ReferenceFrame,
    Direction<Geocentric, Equatorial>: for<'a> From<&'a Direction<Geocentric, F>>,
    Direction<Geocentric, F>: for<'a> From<&'a Direction<Geocentric, Equatorial>>,
    Direction<Geocentric, F>: for<'a> From<&'a Direction<C2, F>>,
    Direction<C1, F>: for<'a> From<&'a Direction<Geocentric, F>>,
{
    #[inline]
    fn transform(&self, jd: JulianDay) -> Direction<C1, F> {
        // Same center No aberration needed.
        if C1::IS_GEOCENTRIC == C2::IS_GEOCENTRIC {
            return Direction::from_vec3(self.as_vec3());
        }

        // Convert source vector to Geocentric+F
        let geo_f: Direction<Geocentric, F> = self.into();

        // Apply or remove aberration depending on the direction.
        let geo_f = if C2::IS_GEOCENTRIC {
            aberrate(&geo_f, jd, remove_aberration_from_direction)
        } else {
            aberrate(&geo_f, jd, apply_aberration_to_direction)
        };

        // Convert back to the target center.
        Direction::<C1, F>::from(&geo_f)
    }
}
