// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! IAU/SOFA ICRS↔Galactic fixed frame rotation.
//!
//! ## Scientific scope
//!
//! The IAU Galactic coordinate system was originally defined in 1958 with
//! respect to FK4 B1950.0. In modern reductions, SOFA follows the Hipparcos
//! realization and uses a direct ICRS↔Galactic rotation that already accounts
//! for the historical FK4/FK5/ICRS interpretation details, including the
//! frame-bias relationship to the J2000.0 mean-place system.
//!
//! ## Technical scope
//!
//! This module mirrors SOFA's `iauIcrs2g` / `iauG2icrs` design: the canonical
//! definition is encoded as a fixed ICRS→Galactic rotation matrix derived from
//! three exact-by-convention angles:
//!
//! - `P = 192.85948°`: ICRS right ascension of the Galactic north pole.
//! - `Q =  27.12825°`: ICRS declination of the Galactic north pole.
//! - `R =  32.93192°`: Galactic longitude of the ascending node of the
//!   Galactic equator on the ICRS equator.
//!
//! SOFA recomputes the matrix from those canonical angles and publishes the
//! resulting elements to 30 decimal places. We store the same matrix as a
//! `Rotation3` constant because it is the actual rotation applied by the SOFA
//! routines, while keeping the canonical angles documented here as the source
//! of the convention.
//!
//! ## References
//!
//! - SOFA `iauIcrs2g` and `iauG2icrs`, release 2023-10-11.
//! - Perryman, M. A. C. & ESA (1997). *The Hipparcos and Tycho catalogues*,
//!   ESA SP-1200.
//! - Blaauw, A., Gum, C. S., Pawsey, J. L., & Westerhout, G. (1960).
//!   "The new I.A.U. system of galactic coordinates". *MNRAS* 121, 123.

use super::bias;
use affn::Rotation3;

/// ICRS right ascension of the Galactic north pole, degrees.
#[allow(dead_code)]
const GALACTIC_NORTH_POLE_RA_DEG: f64 = 192.859_48;

/// ICRS declination of the Galactic north pole, degrees.
#[allow(dead_code)]
const GALACTIC_NORTH_POLE_DEC_DEG: f64 = 27.128_25;

/// Galactic longitude of the ascending node on the ICRS equator, degrees.
#[allow(dead_code)]
const GALACTIC_ASCENDING_NODE_LON_DEG: f64 = 32.931_92;

/// ICRS → Galactic rotation matrix.
///
/// SOFA obtains this by computing:
///
/// ```text
/// R_3(-R) · R_1(π/2 - Q) · R_3(π/2 + P)
/// ```
///
/// from the canonical Hipparcos/IAU angles listed above.
const ICRS_TO_GALACTIC: Rotation3 = Rotation3::from_matrix_unchecked([
    [
        -0.054_875_560_416_215_368_492_398_900_454,
        -0.873_437_090_234_885_048_760_383_168_409,
        -0.483_835_015_548_713_226_831_774_175_116,
    ],
    [
        0.494_109_427_875_583_673_525_222_371_358,
        -0.444_829_629_960_011_178_146_614_061_616,
        0.746_982_244_497_218_890_527_388_004_556,
    ],
    [
        -0.867_666_149_019_004_701_181_616_534_570,
        -0.198_076_373_431_201_528_180_486_091_412,
        0.455_983_776_175_066_922_272_100_478_348,
    ],
]);

/// ICRS → Galactic.
#[inline]
pub(crate) fn icrs_to_galactic() -> Rotation3 {
    ICRS_TO_GALACTIC
}

/// Galactic → ICRS.
#[inline]
pub(crate) fn galactic_to_icrs() -> Rotation3 {
    ICRS_TO_GALACTIC.inverse()
}

/// EquatorialMeanJ2000 → Galactic.
#[inline]
pub(crate) fn equatorial_j2000_to_galactic() -> Rotation3 {
    icrs_to_galactic() * bias::frame_bias_j2000_to_icrs()
}

/// Galactic → EquatorialMeanJ2000.
#[inline]
pub(crate) fn galactic_to_equatorial_j2000() -> Rotation3 {
    equatorial_j2000_to_galactic().inverse()
}

/// EclipticMeanJ2000 → Galactic.
#[inline]
pub(crate) fn ecliptic_j2000_to_galactic() -> Rotation3 {
    icrs_to_galactic() * bias::ecliptic_j2000_to_icrs()
}

/// Galactic → EclipticMeanJ2000.
#[inline]
pub(crate) fn galactic_to_ecliptic_j2000() -> Rotation3 {
    ecliptic_j2000_to_galactic().inverse()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn icrs_galactic_roundtrip_is_identity() {
        let v = [0.3, -0.4, 0.5];
        let fwd = icrs_to_galactic().apply_array(v);
        let back = galactic_to_icrs().apply_array(fwd);
        for k in 0..3 {
            assert!((back[k] - v[k]).abs() < 1e-14);
        }
    }

    #[test]
    fn galactic_north_pole_maps_to_positive_galactic_z() {
        let ra = GALACTIC_NORTH_POLE_RA_DEG.to_radians();
        let dec = GALACTIC_NORTH_POLE_DEC_DEG.to_radians();
        let v = [dec.cos() * ra.cos(), dec.cos() * ra.sin(), dec.sin()];
        let out = icrs_to_galactic().apply_array(v);

        assert!(out[0].abs() < 1e-9);
        assert!(out[1].abs() < 1e-9);
        assert!((out[2] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn matrix_matches_sofa_icrs2g_reference_values() {
        let m = icrs_to_galactic();
        let expected = [
            [
                -0.054_875_560_416_215_368_492_398_900_454,
                -0.873_437_090_234_885_048_760_383_168_409,
                -0.483_835_015_548_713_226_831_774_175_116,
            ],
            [
                0.494_109_427_875_583_673_525_222_371_358,
                -0.444_829_629_960_011_178_146_614_061_616,
                0.746_982_244_497_218_890_527_388_004_556,
            ],
            [
                -0.867_666_149_019_004_701_181_616_534_570,
                -0.198_076_373_431_201_528_180_486_091_412,
                0.455_983_776_175_066_922_272_100_478_348,
            ],
        ];

        for i in 0..3 {
            for j in 0..3 {
                assert!((m.as_matrix()[i][j] - expected[i][j]).abs() < 1e-16);
            }
        }
    }
}
