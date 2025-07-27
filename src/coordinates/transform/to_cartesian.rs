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

impl<C, F, U> crate::coordinates::spherical::SphericalCoord<C, F, U>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
    U: Unit,
{
    pub fn to_cartesian(&self) -> cartesian::Position<C, F, U>
    {
        spherical_to_cartesian_pos(self)
    }
}


impl<C, F, U> crate::coordinates::cartesian::Vector<C, F, U>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
    U: Unit,
{
    pub fn from_spherical(sph: &spherical::Position<C, F, U>) -> Self
    {
        spherical_to_cartesian_pos(sph).into()
    }
}
