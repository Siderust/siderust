use crate::coordinates::{
    centers::ReferenceCenter, frames::ReferenceFrame, kinds::Kind, spherical, cartesian,
};
use crate::units::Unit;

/// Implements conversion from a spherical coordinate to a cartesian coordinate.
///
/// # Formula
/// Given a spherical coordinate `(azimuth, polar, radius)`, the cartesian coordinates `(x, y, z)` are calculated as:
///
/// - `x = r * cos(polar) * cos(azimuth)`
/// - `y = r * cos(polar) * sin(azimuth)`
/// - `z = r * sin(polar)`
impl<C, F, U> From<&spherical::Position<C, F>> for cartesian::Position<C, F, U>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
    U: Unit
{
    fn from(sph: &spherical::Position<C, F>) -> Self {
        let ra_rad = sph.azimuth.to_radians();
        let dec_rad = sph.polar.to_radians();
        let r = sph.distance.expect("SphericalCoord must have a distance");
        let x = r * dec_rad.cos() * ra_rad.cos();
        let y = r * dec_rad.cos() * ra_rad.sin();
        let z = r * dec_rad.sin();
        Self::new(x, y, z)
    }
}

impl<C, F, U> From<&spherical::Direction<C, F, U>> for cartesian::Direction<C, F, U>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
    U: Unit,
{
    fn from(sph: &spherical::Direction<C, F, U>) -> Self {
        let ra_rad = sph.azimuth.to_radians();
        let dec_rad = sph.polar.to_radians();
        let x = dec_rad.cos() * ra_rad.cos();
        let y = dec_rad.cos() * ra_rad.sin();
        let z = dec_rad.sin();
        Self::new(x, y, z)
    }
}

// TODO: can we simply say impl<F> ...?
impl<C, F, U, K> spherical::SphericalCoord<C, F, U, K>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
    U: Unit,
    K: Kind,
    cartesian::Vector<C, F, U, K>: for<'a> From<&'a spherical::SphericalCoord<C, F, U, K>>,
{
    pub fn to_cartesian(&self) -> cartesian::Vector<C, F, U, K> { self.into() }
}

impl<C, F, U, K> spherical::SphericalCoord<C, F, U, K>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
    U: Unit,
    K: Kind,
    spherical::SphericalCoord<C, F, U, K>: for<'a> From<&'a cartesian::Vector<C, F, U, K>>,
{
    pub fn from_cartesian(cart: &cartesian::Vector<C, F, U, K>) -> Self { Self::from(&cart) }
}
