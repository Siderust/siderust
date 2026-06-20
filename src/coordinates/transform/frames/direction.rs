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
//! EclipticMeanJ2000, Galactic) using the constant rotation matrices from
//! [`bias`] plus the IAU ICRS↔Galactic rotation.
//!
//! The higher-level runtime transforms (precession, nutation, GAST) for
//! `Direction` types are provided by the provider system in
//! [`crate::coordinates::transform::ext`].
//!
//! ## References
//!
//! - Capitaine, N. & Wallace, P. T. (2006). *Astronomical Journal*, 132, 2922.
//! - SOFA routines `iauObl06`, `iauBp06`.
//! - Blaauw, A., Gum, C. S., Pawsey, J. L., & Westerhout, G. (1960).
//!   "The new I.A.U. system of galactic coordinates". *MNRAS* 121, 123.

use crate::coordinates::cartesian::Direction;
use crate::coordinates::frames::{self, MutableFrame};
use crate::coordinates::transform::frames::bias;
use crate::coordinates::transform::TransformFrame;
use affn::Rotation3;

/// ICRS → Galactic direction rotation.
///
/// This is the standard orthonormal IAU Galactic rotation used by modern
/// catalogue tooling. It maps ICRS unit vectors to Galactic `(l, b)` axes using
/// the Galactic north pole at approximately `(α, δ) = (192.85948°, 27.12825°)`.
const ICRS_TO_GALACTIC: Rotation3 = Rotation3::from_matrix_unchecked([
    [
        -0.054_875_560_416_215_4,
        -0.873_437_090_234_885_1,
        -0.483_835_015_548_713_2,
    ],
    [
        0.494_109_427_875_583_7,
        -0.444_829_629_960_011_2,
        0.746_982_244_580_286_6,
    ],
    [
        -0.867_666_149_019_004_7,
        -0.198_076_373_431_201_5,
        0.455_983_776_175_066_9,
    ],
]);

#[inline]
fn icrs_to_galactic() -> Rotation3 {
    ICRS_TO_GALACTIC
}

#[inline]
fn galactic_to_icrs() -> Rotation3 {
    ICRS_TO_GALACTIC.inverse()
}

#[inline]
fn equatorial_j2000_to_galactic() -> Rotation3 {
    icrs_to_galactic() * bias::frame_bias_j2000_to_icrs()
}

#[inline]
fn galactic_to_equatorial_j2000() -> Rotation3 {
    equatorial_j2000_to_galactic().inverse()
}

#[inline]
fn ecliptic_j2000_to_galactic() -> Rotation3 {
    icrs_to_galactic() * bias::ecliptic_j2000_to_icrs()
}

#[inline]
fn galactic_to_ecliptic_j2000() -> Rotation3 {
    ecliptic_j2000_to_galactic().inverse()
}

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

/// Frame transform from ICRS to Galactic for directions.
impl TransformFrame<Direction<frames::Galactic>> for Direction<frames::ICRS> {
    fn to_frame(&self) -> Direction<frames::Galactic> {
        let [x, y, z] = icrs_to_galactic().apply_array([self.x(), self.y(), self.z()]);
        Direction::<frames::Galactic>::from_array([x, y, z])
    }
}

/// Frame transform from Galactic to ICRS for directions.
impl TransformFrame<Direction<frames::ICRS>> for Direction<frames::Galactic> {
    fn to_frame(&self) -> Direction<frames::ICRS> {
        let [x, y, z] = galactic_to_icrs().apply_array([self.x(), self.y(), self.z()]);
        Direction::<frames::ICRS>::from_array([x, y, z])
    }
}

/// Frame transform from EquatorialMeanJ2000 to Galactic for directions.
impl TransformFrame<Direction<frames::Galactic>> for Direction<frames::EquatorialMeanJ2000> {
    fn to_frame(&self) -> Direction<frames::Galactic> {
        let [x, y, z] = equatorial_j2000_to_galactic().apply_array([self.x(), self.y(), self.z()]);
        Direction::<frames::Galactic>::from_array([x, y, z])
    }
}

/// Frame transform from Galactic to EquatorialMeanJ2000 for directions.
impl TransformFrame<Direction<frames::EquatorialMeanJ2000>> for Direction<frames::Galactic> {
    fn to_frame(&self) -> Direction<frames::EquatorialMeanJ2000> {
        let [x, y, z] = galactic_to_equatorial_j2000().apply_array([self.x(), self.y(), self.z()]);
        Direction::<frames::EquatorialMeanJ2000>::from_array([x, y, z])
    }
}

/// Frame transform from EclipticMeanJ2000 to Galactic for directions.
impl TransformFrame<Direction<frames::Galactic>> for Direction<frames::EclipticMeanJ2000> {
    fn to_frame(&self) -> Direction<frames::Galactic> {
        let [x, y, z] = ecliptic_j2000_to_galactic().apply_array([self.x(), self.y(), self.z()]);
        Direction::<frames::Galactic>::from_array([x, y, z])
    }
}

/// Frame transform from Galactic to EclipticMeanJ2000 for directions.
impl TransformFrame<Direction<frames::EclipticMeanJ2000>> for Direction<frames::Galactic> {
    fn to_frame(&self) -> Direction<frames::EclipticMeanJ2000> {
        let [x, y, z] = galactic_to_ecliptic_j2000().apply_array([self.x(), self.y(), self.z()]);
        Direction::<frames::EclipticMeanJ2000>::from_array([x, y, z])
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
    use crate::coordinates::spherical;
    use crate::qtty::Degrees;

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

    #[test]
    fn icrs_to_galactic_roundtrip_matches_original() {
        let icrs = Direction::<frames::ICRS>::from_array([0.2, -0.4, 0.7]);
        let galactic: Direction<frames::Galactic> = icrs.to_frame();
        let back: Direction<frames::ICRS> = galactic.to_frame();

        assert_arr3_close(icrs.as_array(), back.as_array(), 1e-12);
    }

    #[test]
    fn galactic_north_pole_maps_to_positive_galactic_z() {
        // IAU Galactic north-pole direction in ICRS coordinates.
        let ra = 192.859_48_f64.to_radians();
        let dec = 27.128_25_f64.to_radians();
        let icrs = Direction::<frames::ICRS>::from_array([
            dec.cos() * ra.cos(),
            dec.cos() * ra.sin(),
            dec.sin(),
        ]);

        let galactic: Direction<frames::Galactic> = icrs.to_frame();
        assert!((galactic.x()).abs() < 1e-9);
        assert!((galactic.y()).abs() < 1e-9);
        assert!((galactic.z() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn spherical_direction_blanket_impl_supports_galactic_frame() {
        let icrs = spherical::Direction::<frames::ICRS>::new_unchecked(
            Degrees::new(27.128_25),
            Degrees::new(192.859_48),
        );
        let galactic: spherical::Direction<frames::Galactic> = icrs.to_frame();

        assert!((galactic.polar.value() - 90.0).abs() < 1e-9);
    }
}
