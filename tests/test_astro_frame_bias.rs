//! Tests for the IAU 2006 frame-bias rotation between GCRS/ICRS and EME2000.
//!
//! These tests were originally in `affn/tests/test_astro_frame_relationships.rs`
//! and were moved here when the frame-bias logic was relocated to
//! `siderust::astro::frame_bias` per the architecture rule that `affn` should
//! remain a pure-geometry crate.
//!
//! Reference: IERS Conventions (2010), §5.4.4 / Table 5.2b.

use affn::frames::{EME2000, GCRS, ICRS};
use affn::ops::Rotation3;
use siderust::astro::frame_bias::{Eme2000FrameBias, GcrsFrameBias, IcrsFrameBias};

/// Convert a `Rotation3` to its rotation angle (radians) via the trace.
fn rotation_angle_rad(r: &Rotation3) -> f64 {
    let m = r.as_matrix();
    let trace = m[0][0] + m[1][1] + m[2][2];
    let cos = ((trace - 1.0) * 0.5).clamp(-1.0, 1.0);
    cos.acos()
}

/// Maximum elementwise absolute difference between two 3×3 rotation matrices.
fn max_abs_diff(a: &Rotation3, b: &Rotation3) -> f64 {
    let am = a.as_matrix();
    let bm = b.as_matrix();
    let mut max = 0.0_f64;
    for i in 0..3 {
        for j in 0..3 {
            let d = (am[i][j] - bm[i][j]).abs();
            if d > max {
                max = d;
            }
        }
    }
    max
}

// ─────────────────────────────────────────────────────────────────────────────
// GCRS ↔ EME2000  (frame bias B; IERS Conventions 2010 §5.4.4 / Table 5.2b)
// ─────────────────────────────────────────────────────────────────────────────

#[test]
fn gcrs_to_eme2000_frame_bias_magnitude_is_about_23_mas() {
    let b = GCRS::frame_bias_to_eme2000();

    // ~23 mas ≈ 1.115e-7 rad. Bound generously (1e-8 .. 5e-7 rad) to lock
    // down the order of magnitude without pinning a specific PFW06/IERS variant.
    let theta = rotation_angle_rad(&b);
    assert!(
        (1.0e-8..5.0e-7).contains(&theta),
        "frame-bias angle {theta:e} rad is outside the documented ≈23 mas range",
    );

    // Tighter sanity: ~23 mas == 23e-3 * pi/648000 rad ≈ 1.115e-7.
    let nominal = 23.0e-3 * std::f64::consts::PI / 648_000.0;
    assert!(
        (theta - nominal).abs() < 0.5 * nominal,
        "frame-bias angle {theta:e} rad differs from the nominal 23 mas \
         ({nominal:e} rad) by more than 50%",
    );
}

#[test]
fn gcrs_to_eme2000_is_not_identity() {
    let b = GCRS::frame_bias_to_eme2000();
    assert_ne!(b, Rotation3::IDENTITY);
    assert!(max_abs_diff(&b, &Rotation3::IDENTITY) > 1.0e-9);
}

#[test]
fn gcrs_eme2000_roundtrip_is_identity_to_1e_minus_15() {
    let fwd = GCRS::frame_bias_to_eme2000();
    let inv = EME2000::frame_bias_to_gcrs();
    let round = fwd * inv;

    let diff = max_abs_diff(&round, &Rotation3::IDENTITY);
    assert!(
        diff < 1.0e-15,
        "GCRS→EME2000→GCRS deviates from identity by {diff:e} (> 1e-15)",
    );

    let round_other = inv * fwd;
    let diff_other = max_abs_diff(&round_other, &Rotation3::IDENTITY);
    assert!(
        diff_other < 1.0e-15,
        "EME2000→GCRS→EME2000 deviates from identity by {diff_other:e}",
    );
}

#[test]
fn icrs_and_gcrs_share_the_same_frame_bias_to_eme2000() {
    // GCRS and ICRS are direction-identical, so their bias matrices must match.
    let from_gcrs = GCRS::frame_bias_to_eme2000();
    let from_icrs = ICRS::frame_bias_to_eme2000();
    assert_eq!(from_gcrs, from_icrs);
}

#[test]
fn bias_off_diagonal_signs_match_sofa() {
    // IERS / SOFA convention for the IAU 2006 bias rb at J2000.0:
    //   rb[0][1] < 0, rb[1][0] > 0,
    //   rb[0][2] > 0, rb[2][0] < 0,
    //   rb[1][2] > 0, rb[2][1] < 0.
    let b = GCRS::frame_bias_to_eme2000();
    let m = b.as_matrix();
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
    assert!(
        m[1][2] > 0.0,
        "rb[1][2] should be positive, got {}",
        m[1][2]
    );
    assert!(
        m[2][1] < 0.0,
        "rb[2][1] should be negative, got {}",
        m[2][1]
    );
}

#[test]
fn eme2000_frame_bias_to_gcrs_and_to_icrs_are_equal() {
    // For direction purposes GCRS and ICRS are bit-identical, so the
    // inverse bias to each must be the same matrix.
    let to_gcrs = EME2000::frame_bias_to_gcrs();
    let to_icrs = EME2000::frame_bias_to_icrs();
    assert_eq!(to_gcrs, to_icrs);
}
