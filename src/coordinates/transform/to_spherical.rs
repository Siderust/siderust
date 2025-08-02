use crate::coordinates::{
    cartesian, spherical,
    centers::*, frames::*,
};
use crate::units::*;

fn cartesian_to_spherical_pos<C, F, U>(
    cart: &cartesian::Position<C, F, U>,
) -> spherical::Position<C, F, U>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
    U: LengthUnit,
{
    let r = cart.distance();
    let x = cart.x();
    let y = cart.y();
    let z = cart.z();

    if r == 0.0 {
        return spherical::Position::CENTER;
    }

    let polar = Degrees::new((z / r).asin().to_degrees());
    let azimuth = Degrees::new(y.value().atan2(x.value()).to_degrees());
    spherical::Position::new_raw(polar, azimuth, cart.distance())
}

fn cartesian_to_spherical_dir<C, F>(
    cart: &cartesian::Direction<C, F>,
) -> spherical::Direction<C, F>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
{
    let x = cart.x().value();
    let y = cart.y().value();
    let z = cart.z().value();

    debug_assert!(
        (cart.distance().value() - 1.0).abs() < 1e-12,
        "A Vector<…, DirectionKind> must have a magnitude ≈ 1.0"
    );

    let polar = Degrees::new(z.asin().to_degrees());
    let azimuth = Degrees::new(y.atan2(x).to_degrees());

    spherical::Direction::new_raw(polar, azimuth, Quantity::<Unitless>::new(1.0))
}

/// Implements conversion from Cartesian to Spherical coordinates
/// by borrowing a `&Vector` reference.
///
/// The conversion uses the following formulas:
/// - `r = sqrt(x² + y² + z²)`
/// - `polar = atan2(y, x)` (angle in the XY plane)
/// - `azimuth = asin(z / r)` (angle from the XY plane to the vector)
///
/// Returns zero values when the Cartesian vector has zero length.
///
/// # Type Parameters
/// - `Center`: The reference center (e.g., Geocentric).
/// - `Frame`: The reference frame (e.g., ICRS).
impl<C, F, U> From<&cartesian::Position<C, F, U>> for spherical::Position<C, F, U>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
    U: LengthUnit,
{
    fn from(cart: &cartesian::Position<C, F, U>) -> Self {
        cartesian_to_spherical_pos(cart)
    }
}

impl<C, F> From<&cartesian::Direction<C, F>> for spherical::Direction<C, F>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
{
    fn from(cart: &cartesian::Direction<C, F>) -> Self {
        cartesian_to_spherical_dir(cart)
    }
}

impl<C, F, U> crate::coordinates::spherical::Position<C, F, U>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
    U: LengthUnit,
{
    pub fn from_cartesian(cart: &cartesian::Position<C, F, U>) -> Self
    {
        cartesian_to_spherical_pos(cart)
    }
}

impl<C, F> crate::coordinates::spherical::Direction<C, F>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
{
    pub fn from_cartesian(cart: &cartesian::Direction<C, F>) -> Self
    {
        cartesian_to_spherical_dir(cart)
    }
}


impl<C, F, U> crate::coordinates::cartesian::Position<C, F, U>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
    U: LengthUnit,
{
    pub fn to_spherical(&self) -> spherical::Position<C, F, U>
    {
        cartesian_to_spherical_pos(self)
    }
}

impl<C, F> crate::coordinates::cartesian::Direction<C, F>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
{
    pub fn to_spherical(&self) -> spherical::Direction<C, F>
    {
        cartesian_to_spherical_dir(self)
    }
}

#[cfg(test)]
mod tests {
    use crate::units::{AstronomicalUnit, Degrees};
    use crate::coordinates::{cartesian, spherical};

    #[test]
    fn test_cartesian_to_spherical() {
        use crate::macros::assert_spherical_eq;
        let cart = cartesian::position::GCRS::<AstronomicalUnit>::new(1.0, 1.0, 1.0);
        let sph: spherical::position::GCRS::<AstronomicalUnit> = cart.to_spherical();
        let expected = spherical::position::GCRS::<AstronomicalUnit>::new(
            Degrees::new(45.0),
            Degrees::new(35.26438968275466),
            1.7320508075688772,
        );
        assert_spherical_eq!(&sph, &expected, 1e-6, "Spherical coordinates do not match expected values");
    }

    #[test]
    fn test_cartesian_spherical_round_trip() {
        use crate::macros::assert_cartesian_eq;
        let cart_original = cartesian::position::GCRS::<AstronomicalUnit>::new(2.0, 3.0, 4.0,);
        let sph = cart_original.to_spherical();
        let cart_converted = cartesian::position::GCRS::<AstronomicalUnit>::from_spherical(&sph);
        assert_cartesian_eq!(&cart_original, &cart_converted, 1e-6, "Cartesian coordinates do not match expected values");
    }

}
