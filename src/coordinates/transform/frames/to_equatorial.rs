use crate::coordinates::{
    CartesianCoord,
    centers::ReferenceCenter,
    frames
};
use crate::coordinates::transform::Transform;

// Implement Transform trait for Ecliptic -> Equatorial
impl<C: ReferenceCenter> Transform<CartesianCoord<C, frames::Equatorial>> for CartesianCoord<C, frames::Ecliptic> {
    fn transform(&self, _jd: crate::units::JulianDay) -> CartesianCoord<C, frames::Equatorial> {
        let eps = 23.439281_f64.to_radians(); // obliquity in radians
        let (sin_e, cos_e) = (eps.sin(), eps.cos());

        let x_eq = self.x();
        let y_eq = cos_e * self.y() - sin_e * self.z();
        let z_eq = sin_e * self.y() + cos_e * self.z();

        CartesianCoord::new(x_eq, y_eq, z_eq)
    }
}

// Implement Transform trait for ICRS -> Equatorial (identity)
impl<C: ReferenceCenter> Transform<CartesianCoord<C, frames::Equatorial>> for CartesianCoord<C, frames::ICRS> {
    fn transform(&self, _jd: crate::units::JulianDay) -> CartesianCoord<C, frames::Equatorial> {
        CartesianCoord::new(self.x(), self.y(), self.z())
    }
}
