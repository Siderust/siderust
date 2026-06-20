// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Single source of truth for **fixed J2000 rotations**.
//!
//! ## Scientific scope
//!
//! The J2000.0 epoch defines several fixed rotations that appear throughout
//! astrometric software:
//!
//! 1. **Frame bias** (`rb`): a ~80 µas rotation from the ICRS to the J2000.0
//!    mean equator/equinox. It arises because the ICRS pole and equinox do not
//!    coincide exactly with the classical J2000.0 definitions.
//! 2. **Mean obliquity** (ε₀ = 84 381.406″): the tilt between Earth's mean
//!    equatorial and ecliptic planes at J2000.0, used to rotate between
//!    equatorial and ecliptic coordinates.
//! 3. **IAU Galactic orientation**: the fixed ICRS↔Galactic rotation used by
//!    modern catalogues and sky-map tooling.
//!
//! ## Technical scope
//!
//! Every fixed frame rotation involving ICRS, EquatorialMeanJ2000,
//! EclipticMeanJ2000, and Galactic is derived from the constants in this
//! module:
//!
//! 1. [`FRAME_BIAS_ICRS_TO_J2000`] – the IAU 2006 frame bias matrix `rb`
//!    (equivalent to the output of SOFA `iauBp06` at J2000.0).
//! 2. The J2000 mean obliquity ε₀ = 84 381.406″ (IAU 2006, `iauObl06`).
//! 3. The orthonormal ICRS→Galactic rotation matrix.
//!
//! Both the legacy `TransformFrame` impls and the `FrameRotationProvider`
//! pipeline use this module so that the fixed rotations cannot diverge.
//!
//! ## References
//!
//! - Capitaine, N. & Wallace, P. T. (2006). *Astronomical Journal*, 132, 2922.
//! - SOFA routines `iauBp06`, `iauObl06`.
//! - Blaauw, A., Gum, C. S., Pawsey, J. L., & Westerhout, G. (1960).
//!   "The new I.A.U. system of galactic coordinates". *MNRAS* 121, 123.

use crate::astro::precession;
use affn::Rotation3;
use std::sync::OnceLock;

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

/// Frame bias rotation matrix from ICRS to mean equator/equinox of J2000.0.
///
/// This is the `rb` matrix produced by SOFA `iauBp06(2451545.0, 0.0, …)`,
/// which constructs the bias via the IAU 2006 Fukushima-Williams
/// parametrisation at epoch J2000.0:
///
///   rb = Rx(−εA) · Rz(−ψb) · Rx(φb) · Rz(γb)
///
/// where γb = −0.052928″, φb = 84381.412819″, ψb = −0.041775″, and
/// εA = obl06(J2000.0) = 84381.406″.
///
/// The off-diagonal signs follow the SOFA convention: `rb[0][1] < 0`.
fn frame_bias_matrix() -> Rotation3 {
    static FRAME_BIAS: OnceLock<Rotation3> = OnceLock::new();
    *FRAME_BIAS.get_or_init(|| precession::precession_matrix_iau2006(crate::J2000))
}

// ── Bias helpers ──────────────────────────────────────────────────────

/// Frame bias: ICRS → EquatorialMeanJ2000.
#[inline]
pub(crate) fn frame_bias_icrs_to_j2000() -> Rotation3 {
    frame_bias_matrix()
}

/// Inverse frame bias: EquatorialMeanJ2000 → ICRS.
#[inline]
pub(crate) fn frame_bias_j2000_to_icrs() -> Rotation3 {
    frame_bias_matrix().inverse()
}

// ── Obliquity helpers ─────────────────────────────────────────────────

/// J2000 mean obliquity ε₀ as `crate::qtty::Radians`.
///
/// 84381.406″ × π / 648000 (IAU 2006, matches SOFA `iauObl06` at J2000.0).
#[inline]
pub(crate) fn j2000_obliquity() -> crate::qtty::Radians {
    crate::qtty::Radians::new(
        precession::J2000_MEAN_OBLIQUITY_ARCSEC * std::f64::consts::PI / 648_000.0,
    )
}

/// Rx(−ε₀): EquatorialMeanJ2000 → EclipticMeanJ2000.
#[inline]
pub(crate) fn obliquity_eq_to_ecl() -> Rotation3 {
    Rotation3::rx(-j2000_obliquity())
}

/// Rx(+ε₀): EclipticMeanJ2000 → EquatorialMeanJ2000.
#[inline]
pub(crate) fn obliquity_ecl_to_eq() -> Rotation3 {
    Rotation3::rx(j2000_obliquity())
}

// ── Galactic helpers ──────────────────────────────────────────────────

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
    icrs_to_galactic() * frame_bias_j2000_to_icrs()
}

/// Galactic → EquatorialMeanJ2000.
#[inline]
pub(crate) fn galactic_to_equatorial_j2000() -> Rotation3 {
    equatorial_j2000_to_galactic().inverse()
}

/// EclipticMeanJ2000 → Galactic.
#[inline]
pub(crate) fn ecliptic_j2000_to_galactic() -> Rotation3 {
    icrs_to_galactic() * ecliptic_j2000_to_icrs()
}

/// Galactic → EclipticMeanJ2000.
#[inline]
pub(crate) fn galactic_to_ecliptic_j2000() -> Rotation3 {
    ecliptic_j2000_to_galactic().inverse()
}

// ── Composed helpers ──────────────────────────────────────────────────

/// ICRS → EclipticMeanJ2000: Rx(−ε₀) · bias.
#[inline]
pub(crate) fn icrs_to_ecliptic_j2000() -> Rotation3 {
    obliquity_eq_to_ecl() * frame_bias_icrs_to_j2000()
}

/// EclipticMeanJ2000 → ICRS: inverse of the above.
#[inline]
pub(crate) fn ecliptic_j2000_to_icrs() -> Rotation3 {
    icrs_to_ecliptic_j2000().inverse()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bias_inverse_roundtrip_is_identity() {
        let forward = frame_bias_icrs_to_j2000();
        let inverse = frame_bias_j2000_to_icrs();
        let v = [0.1, -0.2, 0.3];

        let rotated = forward.apply_array(v);
        let back = inverse.apply_array(rotated);

        let eps = 1e-10;
        assert!((back[0] - v[0]).abs() < eps);
        assert!((back[1] - v[1]).abs() < eps);
        assert!((back[2] - v[2]).abs() < eps);
    }

    #[test]
    fn bias_matrix_is_close_to_identity() {
        let rot = frame_bias_icrs_to_j2000();
        let mat = rot.as_matrix();
        let eps = 1e-7;
        for (i, row) in mat.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((val - expected).abs() < eps);
            }
        }
    }

    /// The bias matrix must match SOFA `bp06` rb to double-precision level.
    ///
    /// Reference values: `iauBp06(2451545.0, 0.0)` → `rb`.
    #[test]
    fn bias_matches_sofa_bp06_rb() {
        let rb = frame_bias_icrs_to_j2000();
        let m = rb.as_matrix();

        // SOFA bp06 rb at J2000.0 (IAU 2006 Fukushima-Williams parametrisation)
        let sofa_rb: [[f64; 3]; 3] = [
            [
                0.999_999_999_999_994_1,
                -7.078_368_960_971_556e-8,
                8.056_213_977_613_186e-8,
            ],
            [
                7.078_368_694_637_676e-8,
                0.999_999_999_999_996_9,
                3.305_943_735_432_137_5e-8,
            ],
            [
                -8.056_214_211_620_057e-8,
                -3.305_943_169_218_395e-8,
                0.999_999_999_999_996_2,
            ],
        ];

        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (m[i][j] - sofa_rb[i][j]).abs() < 1e-15,
                    "rb[{i}][{j}]: siderust {:.17e} vs SOFA {:.17e}",
                    m[i][j],
                    sofa_rb[i][j],
                );
            }
        }
    }

    /// Off-diagonal sign convention: rb[0][1] must be negative (SOFA).
    #[test]
    fn bias_off_diagonal_signs_match_sofa() {
        let rb = frame_bias_icrs_to_j2000();
        let m = rb.as_matrix();

        // SOFA convention: rb[0][1] < 0, rb[1][0] > 0
        assert!(
            m[0][1] < 0.0,
            "rb[0][1] should be negative, got {}",
            m[0][1]
        );
        assert!(
            m[1][0] > 0.0,
            "rb[1][0] should be positive, got {}",
            m[1][0]
        );
        assert!(
            m[0][2] > 0.0,
            "rb[0][2] should be positive, got {}",
            m[0][2]
        );
        assert!(
            m[2][0] < 0.0,
            "rb[2][0] should be negative, got {}",
            m[2][0]
        );
    }

    /// Applying bias to basis vectors must produce sub-mas agreement with SOFA.
    #[test]
    fn bias_on_basis_vectors_agrees_with_sofa() {
        let rb = frame_bias_icrs_to_j2000();

        // SOFA reference outputs for unit basis vectors
        let cases: [([f64; 3], [f64; 3]); 3] = [
            (
                [1.0, 0.0, 0.0],
                [
                    0.999_999_999_999_994_1,
                    7.078_368_694_637_676e-8,
                    -8.056_214_211_620_057e-8,
                ],
            ),
            (
                [0.0, 1.0, 0.0],
                [
                    -7.078_368_960_971_556e-8,
                    0.999_999_999_999_996_9,
                    -3.305_943_169_218_395e-8,
                ],
            ),
            (
                [0.0, 0.0, 1.0],
                [
                    8.056_213_977_613_186e-8,
                    3.305_943_735_432_137_5e-8,
                    0.999_999_999_999_996_2,
                ],
            ),
        ];

        for (vin, expected) in &cases {
            let out = rb.apply_array(*vin);
            for k in 0..3 {
                assert!(
                    (out[k] - expected[k]).abs() < 1e-15,
                    "v={vin:?} component {k}: {:.17e} vs {:.17e}",
                    out[k],
                    expected[k],
                );
            }
        }
    }

    #[test]
    fn obliquity_roundtrip() {
        let v = [0.2, -0.5, 0.7];
        let fwd = obliquity_eq_to_ecl().apply_array(v);
        let back = obliquity_ecl_to_eq().apply_array(fwd);
        for k in 0..3 {
            assert!((back[k] - v[k]).abs() < 1e-14);
        }
    }

    #[test]
    fn icrs_ecliptic_composed_roundtrip() {
        let v = [0.3, -0.4, 0.5];
        let fwd = icrs_to_ecliptic_j2000().apply_array(v);
        let back = ecliptic_j2000_to_icrs().apply_array(fwd);
        for k in 0..3 {
            assert!((back[k] - v[k]).abs() < 1e-14);
        }
    }

    #[test]
    fn icrs_galactic_composed_roundtrip() {
        let v = [0.3, -0.4, 0.5];
        let fwd = icrs_to_galactic().apply_array(v);
        let back = galactic_to_icrs().apply_array(fwd);
        for k in 0..3 {
            assert!((back[k] - v[k]).abs() < 1e-14);
        }
    }

    #[test]
    fn galactic_north_pole_matches_matrix_z_axis() {
        let ra = 192.859_48_f64.to_radians();
        let dec = 27.128_25_f64.to_radians();
        let v = [dec.cos() * ra.cos(), dec.cos() * ra.sin(), dec.sin()];
        let out = icrs_to_galactic().apply_array(v);

        assert!(out[0].abs() < 1e-9);
        assert!(out[1].abs() < 1e-9);
        assert!((out[2] - 1.0).abs() < 1e-12);
    }
}
