use crate::coordinates::{
    CartesianCoord,
    centers::ReferenceCenter,
    frames
};
use crate::coordinates::transform::Transform;

// Implement Transform trait for Ecliptic -> ICRS
impl<C: ReferenceCenter> Transform<CartesianCoord<C, frames::ICRS>> for CartesianCoord<C, frames::Ecliptic> {
    fn transform(&self, _jd: crate::units::JulianDay) -> CartesianCoord<C, frames::ICRS> {
        let eps = 23.439281_f64.to_radians(); // obliquity in radians
        let cos_e = eps.cos();
        let sin_e = eps.sin();

        let x_icrs = self.x();
        let y_icrs = cos_e * self.y() - sin_e * self.z();
        let z_icrs = sin_e * self.y() + cos_e * self.z();

        CartesianCoord::new(x_icrs, y_icrs, z_icrs)
    }
}

// Implement Transform trait for Equatorial -> ICRS (identity)
impl<C: ReferenceCenter> Transform<CartesianCoord<C, frames::ICRS>> for CartesianCoord<C, frames::Equatorial> {
    fn transform(&self, _jd: crate::units::JulianDay) -> CartesianCoord<C, frames::ICRS> {
        CartesianCoord::new(self.x(), self.y(), self.z())
    }
}
