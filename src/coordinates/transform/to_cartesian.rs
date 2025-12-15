use crate::coordinates::{cartesian, centers::ReferenceCenter, frames::ReferenceFrame, spherical};
use qtty::*;

fn spherical_to_cartesian<C, F, U>(
    sph: &spherical::Position<C, F, U>,
) -> cartesian::Position<C, F, U>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
    U: Unit,
{
    let ra_rad = sph.azimuth.to::<Radian>();
    let dec_rad = sph.polar.to::<Radian>();
    let r = sph.distance;
    cartesian::Position::new_with_params(
        sph.center_params().clone(),
        r * dec_rad.cos() * ra_rad.cos(),
        r * dec_rad.cos() * ra_rad.sin(),
        r * dec_rad.sin(),
    )
}

/// Implements conversion from a spherical coordinate to a cartesian coordinate.
///
/// # Formula
/// Given a spherical coordinate `(azimuth, polar, radius)`, the cartesian coordinates `(x, y, z)` are calculated as:
///
/// - `x = r * cos(polar) * cos(azimuth)`
/// - `y = r * cos(polar) * sin(azimuth)`
/// - `z = r * sin(polar)`
impl<C, F, U> From<&spherical::SphericalCoord<C, F, U>> for cartesian::Vector<C, F, U>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
    U: Unit,
{
    fn from(sph: &spherical::SphericalCoord<C, F, U>) -> Self {
        spherical_to_cartesian(sph)
    }
}

impl<C, F, U> crate::coordinates::spherical::SphericalCoord<C, F, U>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
    U: Unit,
{
    pub fn to_cartesian(&self) -> cartesian::Vector<C, F, U> {
        spherical_to_cartesian(self)
    }
}

impl<C, F, U> crate::coordinates::cartesian::Vector<C, F, U>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
    U: Unit,
{
    pub fn from_spherical(sph: &spherical::SphericalCoord<C, F, U>) -> Self {
        spherical_to_cartesian(sph)
    }
}

#[cfg(test)]
mod tests {
    use crate::coordinates::{cartesian, spherical};
    use qtty::{AstronomicalUnit, Degrees};

    #[test]
    fn test_spherical_to_cartesian() {
        use crate::macros::assert_cartesian_eq;
        let sph = spherical::position::GCRS::<AstronomicalUnit>::new(
            Degrees::new(45.0),
            Degrees::new(35.26438968275466),
            1.7320508075688772,
        );
        let cart = cartesian::position::GCRS::<AstronomicalUnit>::from_spherical(&sph);
        let expected = cartesian::position::GCRS::new(1.0, 1.0, 1.0);
        assert_cartesian_eq!(
            &cart,
            &expected,
            1e-6,
            "Cartesian coordinates do not match expected values"
        );
    }

    #[test]
    fn test_spherical_cartesian_round_trip() {
        use crate::macros::assert_spherical_eq;
        let sph_original = spherical::position::GCRS::<AstronomicalUnit>::new(
            Degrees::new(30.0),
            Degrees::new(60.0),
            5.0,
        );
        let cart = cartesian::position::GCRS::<AstronomicalUnit>::from_spherical(&sph_original);
        let sph_converted = cart.to_spherical();
        assert_spherical_eq!(
            &sph_original,
            &sph_converted,
            1e-6,
            "Spherical coordinates do not match expected values"
        );
    }
}
