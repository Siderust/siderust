use super::bias::frame_bias_j2000_to_icrs;
use super::TransformFrame;
use crate::coordinates::{cartesian::Position, centers::ReferenceCenter, frames};
use affn::Rotation3;
use qtty::LengthUnit;

// Implement Transform trait for Ecliptic -> ICRS
impl<C: ReferenceCenter, U: LengthUnit> TransformFrame<Position<C, frames::ICRS, U>>
    for Position<C, frames::Ecliptic, U>
{
    fn to_frame(&self) -> Position<C, frames::ICRS, U> {
        // J2000 mean obliquity ε₀ (IAU 2006): 84381.406″
        let eps = (84381.406_f64 / 3600.0).to_radians();
        let cos_e = eps.cos();
        let sin_e = eps.sin();

        let x0 = self.x().value();
        let y0 = self.y().value();
        let z0 = self.z().value();
        let mean_y = cos_e * y0 - sin_e * z0;
        let mean_z = sin_e * y0 + cos_e * z0;
        let rot: Rotation3 = frame_bias_j2000_to_icrs();
        let [x, y, z] = rot.apply_array([x0, mean_y, mean_z]);
        Position::from_vec3(
            self.center_params().clone(),
            nalgebra::Vector3::new(x.into(), y.into(), z.into()),
        )
    }
}

// Implement Transform trait for EquatorialMeanJ2000 -> ICRS (frame bias inverse)
impl<C: ReferenceCenter, U: LengthUnit> TransformFrame<Position<C, frames::ICRS, U>>
    for Position<C, frames::EquatorialMeanJ2000, U>
{
    fn to_frame(&self) -> Position<C, frames::ICRS, U> {
        let rot: Rotation3 = frame_bias_j2000_to_icrs();
        let [x, y, z] = rot.apply_array([self.x().value(), self.y().value(), self.z().value()]);
        Position::from_vec3(
            self.center_params().clone(),
            nalgebra::Vector3::new(x.into(), y.into(), z.into()),
        )
    }
}
