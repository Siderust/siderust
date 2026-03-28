// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Frame transformations for Direction types.
//!
//! Directions are free unit vectors, so frame transformations are pure rotations
//! without any translation or center dependence.

use crate::coordinates::cartesian::Direction;
use crate::coordinates::frames::{self, MutableFrame};
use crate::coordinates::transform::frames::bias;
use crate::coordinates::transform::TransformFrame;
use affn::Rotation3;
use nalgebra::Vector3;

/// Identity frame transform for directions.
impl<F: MutableFrame> TransformFrame<Direction<F>> for Direction<F> {
    fn to_frame(&self) -> Direction<F> {
        Direction::<F>::from_vec3(self.as_vec3())
    }
}

/// Frame transform from EclipticMeanJ2000 to EquatorialMeanJ2000 for directions.
impl TransformFrame<Direction<frames::EquatorialMeanJ2000>>
    for Direction<frames::EclipticMeanJ2000>
{
    fn to_frame(&self) -> Direction<frames::EquatorialMeanJ2000> {
        let rot = bias::obliquity_ecl_to_eq();
        let [x, y, z] = rot.apply_array([self.x(), self.y(), self.z()]);
        Direction::<frames::EquatorialMeanJ2000>::from_vec3(Vector3::new(x, y, z))
    }
}

/// Frame transform from EquatorialMeanJ2000 to EclipticMeanJ2000 for directions.
impl TransformFrame<Direction<frames::EclipticMeanJ2000>>
    for Direction<frames::EquatorialMeanJ2000>
{
    fn to_frame(&self) -> Direction<frames::EclipticMeanJ2000> {
        let rot = bias::obliquity_eq_to_ecl();
        let [x, y, z] = rot.apply_array([self.x(), self.y(), self.z()]);
        Direction::<frames::EclipticMeanJ2000>::from_vec3(Vector3::new(x, y, z))
    }
}

/// Frame transform from ICRS to EquatorialMeanJ2000 for directions (frame bias).
impl TransformFrame<Direction<frames::EquatorialMeanJ2000>> for Direction<frames::ICRS> {
    fn to_frame(&self) -> Direction<frames::EquatorialMeanJ2000> {
        let rot: Rotation3 = bias::frame_bias_icrs_to_j2000();
        let [x, y, z] = rot.apply_array([self.x(), self.y(), self.z()]);
        Direction::<frames::EquatorialMeanJ2000>::from_vec3(Vector3::new(x, y, z))
    }
}

/// Frame transform from EquatorialMeanJ2000 to ICRS for directions (frame bias).
impl TransformFrame<Direction<frames::ICRS>> for Direction<frames::EquatorialMeanJ2000> {
    fn to_frame(&self) -> Direction<frames::ICRS> {
        let rot: Rotation3 = bias::frame_bias_j2000_to_icrs();
        let [x, y, z] = rot.apply_array([self.x(), self.y(), self.z()]);
        Direction::<frames::ICRS>::from_vec3(Vector3::new(x, y, z))
    }
}

/// Frame transform from EclipticMeanJ2000 to ICRS for directions.
/// Composed via EquatorialMeanJ2000 intermediate frame.
impl TransformFrame<Direction<frames::ICRS>> for Direction<frames::EclipticMeanJ2000> {
    fn to_frame(&self) -> Direction<frames::ICRS> {
        // EclipticMeanJ2000 -> EquatorialMeanJ2000 -> ICRS
        let eq: Direction<frames::EquatorialMeanJ2000> = self.to_frame();
        eq.to_frame()
    }
}

/// Frame transform from ICRS to EclipticMeanJ2000 for directions.
/// Composed via EquatorialMeanJ2000 intermediate frame.
impl TransformFrame<Direction<frames::EclipticMeanJ2000>> for Direction<frames::ICRS> {
    fn to_frame(&self) -> Direction<frames::EclipticMeanJ2000> {
        // ICRS -> EquatorialMeanJ2000 -> EclipticMeanJ2000
        let eq: Direction<frames::EquatorialMeanJ2000> = self.to_frame();
        eq.to_frame()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;

    fn assert_vec3_close(a: &Vector3<f64>, b: &Vector3<f64>, eps: f64) {
        assert!((a.x - b.x).abs() < eps);
        assert!((a.y - b.y).abs() < eps);
        assert!((a.z - b.z).abs() < eps);
    }

    #[test]
    fn identity_direction_transform_is_noop() {
        let dir = Direction::<frames::EclipticMeanJ2000>::from_vec3(Vector3::new(0.0, 0.6, 0.8));
        let same: Direction<frames::EclipticMeanJ2000> = dir.to_frame();
        assert_vec3_close(&dir.as_vec3(), &same.as_vec3(), 1e-15);
    }

    #[test]
    fn ecliptic_to_icrs_roundtrip_matches_original() {
        let dir_ecl =
            Direction::<frames::EclipticMeanJ2000>::from_vec3(Vector3::new(0.1, 0.2, 0.3));
        let icrs: Direction<frames::ICRS> = dir_ecl.to_frame();
        let back: Direction<frames::EclipticMeanJ2000> = icrs.to_frame();

        assert_vec3_close(&dir_ecl.as_vec3(), &back.as_vec3(), 1e-12);
    }

    #[test]
    fn icrs_to_equatorial_bias_and_back_is_stable() {
        let icrs = Direction::<frames::ICRS>::from_vec3(Vector3::new(-0.3, 0.4, -0.5));
        let eq: Direction<frames::EquatorialMeanJ2000> = icrs.to_frame();
        let back: Direction<frames::ICRS> = eq.to_frame();

        assert_vec3_close(&icrs.as_vec3(), &back.as_vec3(), 1e-12);
    }
}
