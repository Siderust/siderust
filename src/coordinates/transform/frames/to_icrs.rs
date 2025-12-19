use super::TransformFrame;
use crate::coordinates::{cartesian::Vector, centers::ReferenceCenter, frames};
use qtty::Unit;

// Implement Transform trait for Ecliptic -> ICRS
impl<C: ReferenceCenter, U: Unit> TransformFrame<Vector<C, frames::ICRS, U>>
    for Vector<C, frames::Ecliptic, U>
{
    fn to_frame(&self) -> Vector<C, frames::ICRS, U> {
        let eps = 23.439281_f64.to_radians(); // obliquity in radians
        let cos_e = eps.cos();
        let sin_e = eps.sin();

        let y = self.y();
        let z = self.z();
        Vector::from_vec3(
            self.center_params().clone(),
            nalgebra::Vector3::new(self.x(), cos_e * y - sin_e * z, sin_e * y + cos_e * z),
        )
    }
}

// Implement Transform trait for Equatorial -> ICRS (identity)
impl<C: ReferenceCenter, U: Unit> TransformFrame<Vector<C, frames::ICRS, U>>
    for Vector<C, frames::Equatorial, U>
{
    fn to_frame(&self) -> Vector<C, frames::ICRS, U> {
        Vector::from_vec3(self.center_params().clone(), self.as_vec3())
    }
}
