// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

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

        let x0 = self.x();
        let y0 = self.y();
        let z0 = self.z();
        let mean_y = cos_e * y0 - sin_e * z0;
        let mean_z = sin_e * y0 + cos_e * z0;
        let rot: Rotation3 = frame_bias_j2000_to_icrs();
        let [x, y, z] = rot * [x0, mean_y, mean_z];
        Position::from_vec3(
            self.center_params().clone(),
            nalgebra::Vector3::new(x, y, z),
        )
    }
}

// Implement Transform trait for EquatorialMeanJ2000 -> ICRS (frame bias inverse)
impl<C: ReferenceCenter, U: LengthUnit> TransformFrame<Position<C, frames::ICRS, U>>
    for Position<C, frames::EquatorialMeanJ2000, U>
{
    fn to_frame(&self) -> Position<C, frames::ICRS, U> {
        let rot: Rotation3 = frame_bias_j2000_to_icrs();
        let [x, y, z] = rot * [self.x(), self.y(), self.z()];
        Position::from_vec3(
            self.center_params().clone(),
            nalgebra::Vector3::new(x, y, z),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coordinates::centers::Barycentric;
    use crate::macros::assert_cartesian_eq;
    use qtty::AstronomicalUnit;

    #[test]
    fn ecliptic_roundtrip_through_icrs_preserves_vector() {
        let ecl = Position::<Barycentric, frames::Ecliptic, AstronomicalUnit>::new(0.0, 1.0, 1.0);
        let icrs: Position<_, frames::ICRS, AstronomicalUnit> = ecl.to_frame();
        let back: Position<_, frames::Ecliptic, AstronomicalUnit> = icrs.to_frame();

        assert_cartesian_eq!(ecl, back, 1e-12);
    }

    #[test]
    fn equatorial_mean_j2000_to_icrs_is_bias_only() {
        let eq = Position::<Barycentric, frames::EquatorialMeanJ2000, AstronomicalUnit>::new(
            1.0, -2.0, 0.5,
        );
        let icrs: Position<_, frames::ICRS, AstronomicalUnit> = eq.to_frame();
        let back: Position<_, frames::EquatorialMeanJ2000, AstronomicalUnit> = icrs.to_frame();

        assert_cartesian_eq!(eq, back, 1e-12);
    }
}
