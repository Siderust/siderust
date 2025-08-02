use crate::coordinates::{
    centers::ReferenceCenter, frames::ReferenceFrame, spherical, cartesian,
};
use crate::units::*;

pub fn spherical_to_cartesian_pos<C, F, U>(
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
    cartesian::Position::new(
        r * dec_rad.cos() * ra_rad.cos(),
        r * dec_rad.cos() * ra_rad.sin(),
        r * dec_rad.sin(),
    )
}

pub fn spherical_to_cartesian_dir<C, F>(
    sph: &spherical::Direction<C, F>,
) -> cartesian::Direction<C, F>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
{
    let ra_rad = sph.azimuth.to::<Radian>();
    let dec_rad = sph.polar.to::<Radian>();

    cartesian::Direction::new(
        dec_rad.cos() * ra_rad.cos(),
        dec_rad.cos() * ra_rad.sin(),
        dec_rad.sin(),
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
impl<C, F, U> From<&spherical::Position<C, F, U>> for cartesian::Position<C, F, U>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
    U: LengthUnit
{
    fn from(sph: &spherical::Position<C, F, U>) -> Self {
        spherical_to_cartesian_pos(sph)
    }
}

impl<C, F> From<&spherical::Direction<C, F>> for cartesian::Direction<C, F>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
{
    fn from(sph: &spherical::Direction<C, F>) -> Self {
        spherical_to_cartesian_dir(sph)
    }
}

impl<C, F, U> crate::coordinates::spherical::Position<C, F, U>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
    U: LengthUnit,
{
    pub fn to_cartesian(&self) -> cartesian::Position<C, F, U>
    {
        spherical_to_cartesian_pos(self)
    }
}

impl<C, F> crate::coordinates::spherical::Direction<C, F>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
{
    pub fn to_cartesian(&self) -> cartesian::Direction<C, F>
    {
        spherical_to_cartesian_dir(self)
    }
}

impl<C, F, U> crate::coordinates::cartesian::Position<C, F, U>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
    U: LengthUnit,
{
    pub fn from_spherical(sph: &spherical::Position<C, F, U>) -> Self
    {
        spherical_to_cartesian_pos(sph).into()
    }
}

impl<C, F> crate::coordinates::cartesian::Direction<C, F>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
{
    pub fn from_spherical(sph: &spherical::Direction<C, F>) -> Self
    {
        spherical_to_cartesian_dir(sph).into()
    }
}


#[cfg(test)]
mod tests {
    use crate::units::{AstronomicalUnit, Degrees};
    use crate::coordinates::{cartesian, spherical};

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
        assert_cartesian_eq!(&cart, &expected, 1e-6, "Cartesian coordinates do not match expected values");
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
        assert_spherical_eq!(&sph_original, &sph_converted, 1e-6, "Spherical coordinates do not match expected values");
    }
}
