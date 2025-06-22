use crate::coordinates::{
    cartesian::Vector,
    centers::ReferenceCenter,
    kinds::Kind,
    frames
};
use crate::coordinates::transform::Transform;
use crate::units::Unit;

// Implement Transform trait for Ecliptic -> ICRS
impl<C: ReferenceCenter, K: Kind, U: Unit> Transform<Vector<C, frames::ICRS, U, K>> for Vector<C, frames::Ecliptic, U, K> {
    fn transform(&self, _jd: crate::units::JulianDay) -> Vector<C, frames::ICRS, U, K> {
        let eps = 23.439281_f64.to_radians(); // obliquity in radians
        let cos_e = eps.cos();
        let sin_e = eps.sin();

        let x_icrs = self.x();
        let y_icrs = cos_e * self.y() - sin_e * self.z();
        let z_icrs = sin_e * self.y() + cos_e * self.z();

        Vector::new(x_icrs, y_icrs, z_icrs)
    }
}

// Implement Transform trait for Equatorial -> ICRS (identity)
impl<C: ReferenceCenter, K: Kind, U: Unit> Transform<Vector<C, frames::ICRS, U, K>> for Vector<C, frames::Equatorial, U, K> {
    fn transform(&self, _jd: crate::units::JulianDay) -> Vector<C, frames::ICRS, U, K> {
        Vector::new(self.x(), self.y(), self.z())
    }
}
