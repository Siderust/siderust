use crate::coordinates::{
    centers::ReferenceCenter, frames::ReferenceFrame, spherical, cartesian,
};
use crate::units::*;

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
        let ra_rad = sph.azimuth.to::<Radian>();
        let dec_rad = sph.polar.to::<Radian>();
        let r = sph.distance.expect("SphericalCoord must have a distance");
        let x = r * dec_rad.cos() * ra_rad.cos();
        let y = r * dec_rad.cos() * ra_rad.sin();
        let z = r * dec_rad.sin();
        Self::new(x, y, z)
    }
}

impl<C, F> From<&spherical::Direction<C, F>> for cartesian::Direction<C, F>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
{
    fn from(sph: &spherical::Direction<C, F>) -> Self {
        let ra_rad = sph.azimuth.to::<Radian>();
        let dec_rad = sph.polar.to::<Radian>();

        Self::new(
            Quantity::<f64>::new(dec_rad.cos() * ra_rad.cos()),
            Quantity::<f64>::new(dec_rad.cos() * ra_rad.sin()),
            Quantity::<f64>::new(dec_rad.sin())
        )
    }
}

// TODO: can we simply say impl<F> ...?
impl<C, F, U> spherical::SphericalCoord<C, F, U>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
    U: Unit,
    cartesian::Vector<C, F, U>: for<'a> From<&'a spherical::SphericalCoord<C, F, U>>,
{
    pub fn to_cartesian(&self) -> cartesian::Vector<C, F, U> { self.into() }
}

impl<C, F, U> spherical::SphericalCoord<C, F, U>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
    U: Unit,
    spherical::SphericalCoord<C, F, U>: for<'a> From<&'a cartesian::Vector<C, F, U>>,
{
    pub fn from_cartesian(cart: &cartesian::Vector<C, F, U>) -> Self { Self::from(&cart) }
}
