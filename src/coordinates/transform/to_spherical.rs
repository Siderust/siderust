use crate::coordinates::{cartesian, centers::*, frames::*, spherical};
use qtty::*;

fn cartesian_to_spherical<C, F, U>(
    cart: &cartesian::Vector<C, F, U>,
) -> spherical::Position<C, F, U>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
    U: Unit,
{
    let r = cart.distance();
    let x = cart.x().value();
    let y = cart.y().value();
    let z = cart.z().value();
    let r_val = r.value();

    if r_val == 0.0 {
        return spherical::Position::new_raw_with_params(
            cart.center_params().clone(),
            Degrees::new(0.0),
            Degrees::new(0.0),
            Quantity::<U>::new(0.0),
        );
    }

    let polar = Degrees::new((z / r_val).asin().to_degrees());
    let azimuth = Degrees::new(y.atan2(x).to_degrees());
    spherical::Position::new_raw_with_params(cart.center_params().clone(), polar, azimuth, r)
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
impl<C, F, U> From<&cartesian::Vector<C, F, U>> for spherical::Position<C, F, U>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
    U: Unit,
{
    fn from(cart: &cartesian::Vector<C, F, U>) -> Self {
        cartesian_to_spherical(cart)
    }
}

impl<C, F, U> crate::coordinates::spherical::Position<C, F, U>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
    U: Unit,
{
    pub fn from_cartesian(cart: &cartesian::Vector<C, F, U>) -> Self {
        cartesian_to_spherical(cart)
    }
}

impl<C, F, U> crate::coordinates::cartesian::Vector<C, F, U>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
    U: Unit,
{
    pub fn to_spherical(&self) -> spherical::Position<C, F, U> {
        cartesian_to_spherical(self)
    }
}

#[cfg(test)]
mod tests {
    use crate::coordinates::{cartesian, spherical};
    use qtty::{AstronomicalUnit, Degrees};

    #[test]
    fn test_cartesian_to_spherical() {
        use crate::macros::assert_spherical_eq;
        let cart = cartesian::position::GCRS::<AstronomicalUnit>::new(1.0, 1.0, 1.0);
        let sph: spherical::position::GCRS<AstronomicalUnit> = cart.to_spherical();
        let expected = spherical::position::GCRS::<AstronomicalUnit>::new(
            Degrees::new(45.0),
            Degrees::new(35.26438968275466),
            1.7320508075688772,
        );
        assert_spherical_eq!(
            &sph,
            &expected,
            1e-6,
            "Spherical coordinates do not match expected values"
        );
    }

    #[test]
    fn test_cartesian_spherical_round_trip() {
        use crate::macros::assert_cartesian_eq;
        let cart_original = cartesian::position::GCRS::<AstronomicalUnit>::new(2.0, 3.0, 4.0);
        let sph = cart_original.to_spherical();
        let cart_converted = cartesian::position::GCRS::<AstronomicalUnit>::from_spherical(&sph);
        assert_cartesian_eq!(
            &cart_original,
            &cart_converted,
            1e-6,
            "Cartesian coordinates do not match expected values"
        );
    }
}
