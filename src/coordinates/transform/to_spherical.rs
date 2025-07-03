use crate::units::Degrees;
use crate::coordinates::{
    cartesian, spherical,
    centers::*, frames::*,
    kinds::Kind,
};
use crate::units::Unit;

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
    U: Unit,
{
    fn from(cart: &cartesian::Position<C, F, U>) -> Self {
        let r:f64 = cart.distance().into();
        let x:f64 = cart.x().into();
        let y:f64 = cart.y().into();
        let z:f64 = cart.z().into();

        if r == 0.0 {
            return Self::CENTER;
        }

        let polar = Degrees::new((z / r).asin().to_degrees());
        let azimuth = Degrees::new(y.atan2(x).to_degrees());
        Self::new_spherical_coord(polar, azimuth, Some(cart.distance()))
    }
}

impl<C, F> From<&cartesian::Direction<C, F>> for spherical::Direction<C, F>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
{
    fn from(cart: &cartesian::Direction<C, F>) -> Self {
        let x:f64 = cart.x().into();
        let y:f64 = cart.y().into();
        let z:f64 = cart.z().into();

        debug_assert!(
            (cart.as_vec3().magnitude() - 1.0).abs() < 1e-12,
            "A Vector<…, DirectionKind> must have a magnitude ≈ 1.0"
        );

        let polar   = Degrees::new(z.asin().to_degrees());
        let azimuth = Degrees::new(y.atan2(x).to_degrees());

        Self::new_spherical_coord(polar, azimuth, None)
    }
}

impl<C, F, U, K> cartesian::Vector<C, F, U, K>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
    U: Unit,
    K: Kind,
    spherical::SphericalCoord<C, F, U, K>: for<'a> From<&'a cartesian::Vector<C, F, U, K>>,
{
    pub fn to_spherical(&self) -> spherical::SphericalCoord<C, F, U, K> { self.into() }
}


impl<C, F, U, K> cartesian::Vector<C, F, U, K>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
    U: Unit,
    K: Kind,
    cartesian::Vector<C, F, U, K>: for<'a> From<&'a spherical::SphericalCoord<C, F, U, K>>,
{
    pub fn from_spherical(sph: &spherical::SphericalCoord<C, F, U, K>) -> Self { Self::from(&sph) }
}

#[cfg(test)]
mod tests {
    use crate::units::Degrees;
    use crate::macros::{assert_cartesian_eq, assert_spherical_eq};
    use crate::coordinates::{cartesian, spherical};

    #[test]
    fn test_cartesian_to_spherical() {
        let cart = cartesian::position::GCRS::new(1.0, 1.0, 1.0);
        let sph: spherical::position::GCRS = cart.to_spherical();
        let expected = spherical::position::GCRS::new(
            Degrees::new(45.0),
            Degrees::new(35.26438968275466),
            1.7320508075688772,
        );
        assert_spherical_eq!(&sph, &expected, 1e-6, "Spherical coordinates do not match expected values");
    }

    #[test]
    fn test_spherical_to_cartesian() {
        let sph = spherical::position::GCRS::new(
            Degrees::new(45.0),
            Degrees::new(35.26438968275466), 
            1.7320508075688772,
        );
        let cart = cartesian::position::GCRS::from_spherical(&sph);
        let expected = cartesian::position::GCRS::new(1.0, 1.0, 1.0);
        assert_cartesian_eq!(&cart, &expected, 1e-6, "Cartesian coordinates do not match expected values");
    }

    #[test]
    fn test_cartesian_spherical_round_trip() {
        let cart_original = cartesian::position::GCRS::new(2.0, 3.0, 4.0,);
        let sph = cart_original.to_spherical();
        let cart_converted = cartesian::Position::from_spherical(&sph);
        assert_cartesian_eq!(&cart_original, &cart_converted, 1e-6, "Cartesian coordinates do not match expected values");
    }

    #[test]
    fn test_spherical_cartesian_round_trip() {
        let sph_original = spherical::position::GCRS::new(
            Degrees::new(30.0),
            Degrees::new(60.0),
            5.0,
        );
        let cart = cartesian::Vector::from_spherical(&sph_original);
        let sph_converted = cart.to_spherical();
        assert_spherical_eq!(&sph_original, &sph_converted, 1e-6, "Spherical coordinates do not match expected values");
    }
}
