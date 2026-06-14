// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Frame rotations for [`Direction`](crate::coordinates::cartesian::Direction) types.
//!
//! ## Scientific scope
//!
//! Directions are unit vectors representing a pointing on the celestial sphere.
//! They are **translation-invariant**: only the orientation of the coordinate
//! axes matters, not the origin. Consequently, frame transformations for
//! directions are pure rotations — no center-shift correction is needed or
//! defined.
//!
//! ## Technical scope
//!
//! This module provides `TransformFrame` impls for `cartesian::Direction<F>`
//! between the J2000.0 fixed frames (ICRS, EquatorialMeanJ2000,
//! EclipticMeanJ2000) using the constant rotation matrices from [`bias`].
//!
//! The higher-level runtime transforms (precession, nutation, GAST) for
//! `Direction` types are provided by the provider system in
//! [`crate::coordinates::transform::ext`].
//!
//! ## References
//!
//! - Capitaine, N. & Wallace, P. T. (2006). *Astronomical Journal*, 132, 2922.
//! - SOFA routines `iauObl06`, `iauBp06`.

use crate::coordinates::cartesian::Direction;
use crate::coordinates::frames::{self, MutableFrame};
use crate::coordinates::transform::frames::bias;
use crate::coordinates::transform::TransformFrame;
use affn::Rotation3;

/// Identity frame transform for directions.
impl<F: MutableFrame> TransformFrame<Direction<F>> for Direction<F> {
    fn to_frame(&self) -> Direction<F> {
        Direction::<F>::from_array(self.as_array())
    }
}

/// Frame transform from EclipticMeanJ2000 to EquatorialMeanJ2000 for directions.
impl TransformFrame<Direction<frames::EquatorialMeanJ2000>>
    for Direction<frames::EclipticMeanJ2000>
{
    fn to_frame(&self) -> Direction<frames::EquatorialMeanJ2000> {
        let rot = bias::obliquity_ecl_to_eq();
        let [x, y, z] = rot.apply_array([self.x(), self.y(), self.z()]);
        Direction::<frames::EquatorialMeanJ2000>::from_array([x, y, z])
    }
}

/// Frame transform from EquatorialMeanJ2000 to EclipticMeanJ2000 for directions.
impl TransformFrame<Direction<frames::EclipticMeanJ2000>>
    for Direction<frames::EquatorialMeanJ2000>
{
    fn to_frame(&self) -> Direction<frames::EclipticMeanJ2000> {
        let rot = bias::obliquity_eq_to_ecl();
        let [x, y, z] = rot.apply_array([self.x(), self.y(), self.z()]);
        Direction::<frames::EclipticMeanJ2000>::from_array([x, y, z])
    }
}

/// Frame transform from ICRS to EquatorialMeanJ2000 for directions (frame bias).
impl TransformFrame<Direction<frames::EquatorialMeanJ2000>> for Direction<frames::ICRS> {
    fn to_frame(&self) -> Direction<frames::EquatorialMeanJ2000> {
        let rot: Rotation3 = bias::frame_bias_icrs_to_j2000();
        let [x, y, z] = rot.apply_array([self.x(), self.y(), self.z()]);
        Direction::<frames::EquatorialMeanJ2000>::from_array([x, y, z])
    }
}

/// Frame transform from EquatorialMeanJ2000 to ICRS for directions (frame bias).
impl TransformFrame<Direction<frames::ICRS>> for Direction<frames::EquatorialMeanJ2000> {
    fn to_frame(&self) -> Direction<frames::ICRS> {
        let rot: Rotation3 = bias::frame_bias_j2000_to_icrs();
        let [x, y, z] = rot.apply_array([self.x(), self.y(), self.z()]);
        Direction::<frames::ICRS>::from_array([x, y, z])
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

    fn assert_arr3_close(a: [f64; 3], b: [f64; 3], eps: f64) {
        assert!((a[0] - b[0]).abs() < eps);
        assert!((a[1] - b[1]).abs() < eps);
        assert!((a[2] - b[2]).abs() < eps);
    }

    #[test]
    fn identity_direction_transform_is_noop() {
        let dir = Direction::<frames::EclipticMeanJ2000>::from_array([0.0, 0.6, 0.8]);
        let same: Direction<frames::EclipticMeanJ2000> = dir.to_frame();
        assert_arr3_close(dir.as_array(), same.as_array(), 1e-15);
    }

    #[test]
    fn ecliptic_to_icrs_roundtrip_matches_original() {
        let dir_ecl = Direction::<frames::EclipticMeanJ2000>::from_array([0.1, 0.2, 0.3]);
        let icrs: Direction<frames::ICRS> = dir_ecl.to_frame();
        let back: Direction<frames::EclipticMeanJ2000> = icrs.to_frame();

        assert_arr3_close(dir_ecl.as_array(), back.as_array(), 1e-12);
    }

    #[test]
    fn icrs_to_equatorial_bias_and_back_is_stable() {
        let icrs = Direction::<frames::ICRS>::from_array([-0.3, 0.4, -0.5]);
        let eq: Direction<frames::EquatorialMeanJ2000> = icrs.to_frame();
        let back: Direction<frames::ICRS> = eq.to_frame();

        assert_arr3_close(icrs.as_array(), back.as_array(), 1e-12);
    }
}
