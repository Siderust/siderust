use crate::coordinates::{
    cartesian::Vector,
    centers::ReferenceCenter,
    kinds::Kind,
    frames
};
use crate::coordinates::transform::Transform;
use crate::units::LengthUnit;

// Implement Transform trait for Ecliptic -> ICRS
impl<C: ReferenceCenter, K: Kind, U: LengthUnit> Transform<Vector<C, frames::ICRS, U, K>> for Vector<C, frames::Ecliptic, U, K> {
    fn transform(&self, _jd: crate::units::JulianDay) -> Vector<C, frames::ICRS, U, K> {
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
impl<C: ReferenceCenter, K: Kind, U: LengthUnit> Transform<Vector<C, frames::ICRS, U, K>> for Vector<C, frames::Equatorial, U, K> {
    fn transform(&self, _jd: crate::units::JulianDay) -> Vector<C, frames::ICRS, U, K> {
        Vector::new(self.x(), self.y(), self.z())
    }
}
