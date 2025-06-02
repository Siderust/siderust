use crate::units::Degrees;
use crate::coordinates::{
    CartesianCoord, SphericalCoord,
    centers::*, frames::*
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
impl<Center, Frame> From<&CartesianCoord<Center, Frame>> for SphericalCoord<Center, Frame>
where
    Center: ReferenceCenter,
    Frame: ReferenceFrame,
{
    fn from(cart: &CartesianCoord<Center, Frame>) -> Self {
        let r = cart.distance_from_origin();
        if r == 0.0 {
            return SphericalCoord::<Center, Frame>::new_spherical_coord(
                Degrees::new(0.0),
                Degrees::new(0.0),
                0.0
            );
        }

        let polar = Degrees::new((cart.z() / r).asin().to_degrees());
        let azimuth = Degrees::new(cart.y().atan2(cart.x()).to_degrees());

        SphericalCoord::<Center, Frame>::new_spherical_coord(polar, azimuth, r)
    }
}

impl<Center, Frame> CartesianCoord<Center, Frame>
where
    Center: ReferenceCenter,
    Frame: ReferenceFrame,
{
    /// Converts this Cartesian coordinate into its equivalent spherical representation.
    pub fn to_spherical(&self) -> SphericalCoord<Center, Frame> {
        self.into()
    }

    /// Constructs a Cartesian coordinate from a spherical one.
    pub fn from_spherical(sph: &SphericalCoord<Center, Frame>) -> Self {
        sph.into()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::units::Degrees;

    #[test]
    fn test_cartesian_to_spherical() {
        let cart = CartesianCoord::<Geocentric, ICRS>::new(1.0, 1.0, 1.0);
        let sph: SphericalCoord<Geocentric, ICRS> = cart.to_spherical();

        assert!((sph.distance - 1.7320508075688772).abs() < 1e-6);
        assert!((sph.azimuth - Degrees::new(45.0)).abs() < Degrees::new(1e-6));
        assert!((sph.polar - Degrees::new(35.26438968275466)).abs() < Degrees::new(1e-6));
    }

    #[test]
    fn test_spherical_to_cartesian() {
        let sph = SphericalCoord::<Geocentric, ICRS>::new(
            Degrees::new(45.0),  
            Degrees::new(35.26438968275466), 
            1.7320508075688772,
        );
        let cart: CartesianCoord<Geocentric, ICRS> = CartesianCoord::from_spherical(&sph);

        assert!((cart.x() - 1.0).abs() < 1e-6);
        assert!((cart.y() - 1.0).abs() < 1e-6);
        assert!((cart.z() - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_cartesian_spherical_round_trip() {
        let cart_original = CartesianCoord::<Geocentric, ICRS>::new(2.0, 3.0, 4.0,);
        let sph = cart_original.to_spherical();
        let cart_converted = CartesianCoord::from_spherical(&sph);

        assert!((cart_original.x() - cart_converted.x()).abs() < 1e-6);
        assert!((cart_original.y() - cart_converted.y()).abs() < 1e-6);
        assert!((cart_original.z() - cart_converted.z()).abs() < 1e-6);
    }

    #[test]
    fn test_spherical_cartesian_round_trip() {
        let sph_original = SphericalCoord::<Geocentric, ICRS>::new(
            Degrees::new(30.0),
            Degrees::new(60.0),
            5.0,
        );
        let cart = CartesianCoord::from_spherical(&sph_original);
        let sph_converted = cart.to_spherical();

        assert!((sph_original.distance - sph_converted.distance).abs() < 1e-6);
        assert!((sph_original.azimuth - sph_converted.azimuth).abs() < Degrees::new(1e-6));
        assert!((sph_original.polar - sph_converted.polar).abs() < Degrees::new(1e-6));
    }
}
