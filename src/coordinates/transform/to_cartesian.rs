use crate::coordinates::{
    CartesianCoord, SphericalCoord, SphericalBuilder,
    frames::ReferenceFrame,
    centers::ReferenceCenter
};

/// Implements conversion from a spherical coordinate to a cartesian coordinate.
///
/// # Formula
/// Given a spherical coordinate `(azimuth, polar, radius)`, the cartesian coordinates `(x, y, z)` are calculated as:
///
/// - `x = r * cos(polar) * cos(azimuth)`
/// - `y = r * cos(polar) * sin(azimuth)`
/// - `z = r * sin(polar)`
impl<Center, Frame> From<&SphericalCoord<Center, Frame>> for CartesianCoord<Center, Frame>
where
    Center: ReferenceCenter,
    Frame: ReferenceFrame,
{
    fn from(sph: &SphericalCoord<Center, Frame>) -> Self {
        let ra_rad = sph.azimuth.to_radians();
        let dec_rad = sph.polar.to_radians();
        let r = sph.radial_distance;
        let x = r * dec_rad.cos() * ra_rad.cos();
        let y = r * dec_rad.cos() * ra_rad.sin();
        let z = r * dec_rad.sin();
        CartesianCoord::new(x, y, z)
    }
}

impl<Center, Frame> SphericalCoord<Center, Frame>
where
    Center: ReferenceCenter,
    Frame: ReferenceFrame,
    SphericalCoord<Center, Frame>: SphericalBuilder<Center, Frame>,
{
    pub fn to_cartesian(&self) -> CartesianCoord<Center, Frame> { self.into() }
    pub fn from_cartesian(cart: &CartesianCoord<Center, Frame>) -> Self { cart.into() }
}
