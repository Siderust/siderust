use crate::units::Degrees;
use crate::coordinates::{
    cartesian, spherical,
    centers::*, frames::*,
    kinds::Kind,
};


/// Implements conversion from Cartesian to Spherical coordinates
/// by borrowing a `&CartesianCoord` reference.
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
impl<C, F> From<&cartesian::Position<C, F>> for spherical::Position<C, F>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
{
    fn from(cart: &cartesian::Position<C, F>) -> Self {
        let r = cart.distance_from_origin();
        if r == 0.0 {
            return Self::CENTER;
        }

        let polar = Degrees::new((cart.z() / r).asin().to_degrees());
        let azimuth = Degrees::new(cart.y().atan2(cart.x()).to_degrees());
        Self::new_spherical_coord(polar, azimuth, Some(r))
    }
}

impl<C, F> From<&cartesian::Direction<C, F>> for spherical::Direction<C, F>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
{
    fn from(cart: &cartesian::Direction<C, F>) -> Self {
        debug_assert!(
            (cart.distance_from_origin() - 1.0).abs() < 1e-12,
            "A CartesianCoord<…, DirectionKind> must have a magnitude ≈ 1.0"
        );

        let polar   = Degrees::new(cart.z().asin().to_degrees());
        let azimuth = Degrees::new(cart.y().atan2(cart.x()).to_degrees());

        Self::new_spherical_coord(polar, azimuth, None)
    }
}

impl<C, F, K> cartesian::CartesianCoord<C, F, K>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
    K: Kind,
{
    /// Converts this Cartesian coordinate into its equivalent spherical representation.
    pub fn to_spherical(&self)   -> spherical::SphericalCoord<C, F, K> { self.into() }
    pub fn from_spherical(&self) -> spherical::SphericalCoord<C, F, K> { self.into() }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::units::Degrees;

    #[test]
    fn test_cartesian_to_spherical() {
        let cart = cartesian::Position::<Geocentric, ICRS>::new(1.0, 1.0, 1.0);
        let sph: spherical::Position<Geocentric, ICRS> = cart.to_spherical();

        assert!((sph.distance - 1.7320508075688772).abs() < 1e-6);
        assert!((sph.azimuth - Degrees::new(45.0)).abs() < Degrees::new(1e-6));
        assert!((sph.polar - Degrees::new(35.26438968275466)).abs() < Degrees::new(1e-6));
    }

    #[test]
    fn test_spherical_to_cartesian() {
        let sph = spherical::Position::<Geocentric, ICRS>::new(
            Degrees::new(45.0),  
            Degrees::new(35.26438968275466), 
            1.7320508075688772,
        );
        let cart = cartesian::Position::<Geocentric, ICRS>::from_spherical(&sph);

        assert!((cart.x() - 1.0).abs() < 1e-6);
        assert!((cart.y() - 1.0).abs() < 1e-6);
        assert!((cart.z() - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_cartesian_spherical_round_trip() {
        let cart_original = cartesian::Position::<Geocentric, ICRS>::new(2.0, 3.0, 4.0,);
        let sph = cart_original.to_spherical();
        let cart_converted = cartesian::Position::from_spherical(&sph);

        assert!((cart_original.x() - cart_converted.x()).abs() < 1e-6);
        assert!((cart_original.y() - cart_converted.y()).abs() < 1e-6);
        assert!((cart_original.z() - cart_converted.z()).abs() < 1e-6);
    }

    #[test]
    fn test_spherical_cartesian_round_trip() {
        let sph_original = spherical::Position::<Geocentric, ICRS>::new(
            Degrees::new(30.0),
            Degrees::new(60.0),
            5.0,
        );
        let cart = cartesian::CartesianCoord::from_spherical(&sph_original);
        let sph_converted = cart.to_spherical();

        assert!((sph_original.distance - sph_converted.distance).abs() < 1e-6);
        assert!((sph_original.azimuth - sph_converted.azimuth).abs() < Degrees::new(1e-6));
        assert!((sph_original.polar - sph_converted.polar).abs() < Degrees::new(1e-6));
    }
}
