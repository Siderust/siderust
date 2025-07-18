use crate::coordinates::{
    cartesian::Vector,
    centers::ReferenceCenter,
    frames
};
use crate::coordinates::transform::Transform;
use crate::units::Unit;

// Implement Transform trait for Ecliptic -> ICRS
impl<C: ReferenceCenter, U: Unit> Transform<Vector<C, frames::ICRS, U>> for Vector<C, frames::Ecliptic, U> {
    fn transform(&self, _jd: crate::astro::JulianDate) -> Vector<C, frames::ICRS, U> {
        let eps = 23.439281_f64.to_radians(); // obliquity in radians
        let cos_e = eps.cos();
        let sin_e = eps.sin();

        let y = self.y();
        let z = self.z();
        Vector::new(
            self.x(),
            cos_e * y - sin_e * z,
            sin_e * y + cos_e * z
        )
    }
}

// Implement Transform trait for Equatorial -> ICRS (identity)
impl<C: ReferenceCenter, U: Unit> Transform<Vector<C, frames::ICRS, U>> for Vector<C, frames::Equatorial, U> {
    fn transform(&self, _jd: crate::astro::JulianDate) -> Vector<C, frames::ICRS, U> {
        Vector::new(self.x(), self.y(), self.z())
    }
}
