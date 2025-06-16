use crate::coordinates::{
    centers::ReferenceCenter, frames::ReferenceFrame, kinds::Kind, spherical, cartesian,
};

/// Implements conversion from a spherical coordinate to a cartesian coordinate.
///
/// # Formula
/// Given a spherical coordinate `(azimuth, polar, radius)`, the cartesian coordinates `(x, y, z)` are calculated as:
///
/// - `x = r * cos(polar) * cos(azimuth)`
/// - `y = r * cos(polar) * sin(azimuth)`
/// - `z = r * sin(polar)`
impl<C, F> From<&spherical::Position<C, F>> for cartesian::Position<C, F>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
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

impl<C, F> From<&spherical::Direction<C, F>> for cartesian::Direction<C, F>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
{
    fn from(sph: &spherical::Direction<C, F>) -> Self {
        let ra_rad = sph.azimuth.to_radians();
        let dec_rad = sph.polar.to_radians();
        let x = dec_rad.cos() * ra_rad.cos();
        let y = dec_rad.cos() * ra_rad.sin();
        let z = dec_rad.sin();
        Self::new(x, y, z)
    }
}

impl<C, F, K> spherical::SphericalCoord<C, F, K>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
    K: Kind,
    cartesian::CartesianCoord<C, F, K>: for<'a> From<&'a spherical::SphericalCoord<C, F, K>>,
{
    pub fn to_cartesian(&self) -> cartesian::CartesianCoord<C, F, K> { self.into() }
}

/*
impl<C: ReferenceCenter, F: ReferenceFrame> spherical::Position<C, F>
{
    pub fn to_cartesian(&self) -> cartesian::Position<C, F> { self.into() }
}

impl<C: ReferenceCenter, F: ReferenceFrame> spherical::Direction<C, F>
{
    pub fn to_cartesian(&self) -> cartesian::Direction<C, F> { self.into() }
} */

impl<C, F, K> spherical::SphericalCoord<C, F, K>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
    K: Kind,
    spherical::SphericalCoord<C, F, K>: for<'a> From<&'a cartesian::CartesianCoord<C, F, K>>,
{
    pub fn from_cartesian(cart: &cartesian::CartesianCoord<C, F, K>) -> Self { Self::from(&cart) }
}
