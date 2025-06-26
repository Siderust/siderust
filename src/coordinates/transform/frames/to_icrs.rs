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
        let x:f64 = self.x().into();
        let y:f64 = self.y().into();
        let z:f64 = self.z().into();
        let eps = 23.439281_f64.to_radians(); // obliquity in radians
        let cos_e = eps.cos();
        let sin_e = eps.sin();

        let x_icrs = x;
        let y_icrs = cos_e * y - sin_e * z;
        let z_icrs = sin_e * y + cos_e * z;

        Vector::new(x_icrs.into(), y_icrs.into(), z_icrs.into())
    }
}

// Implement Transform trait for Equatorial -> ICRS (identity)
impl<C: ReferenceCenter, K: Kind, U: Unit> Transform<Vector<C, frames::ICRS, U, K>> for Vector<C, frames::Equatorial, U, K> {
    fn transform(&self, _jd: crate::units::JulianDay) -> Vector<C, frames::ICRS, U, K> {
        Vector::new(self.x(), self.y(), self.z())
    }
}
