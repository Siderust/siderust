// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Angular Separation Audit Tests
//!
//! Verifies the `angular_separation` function across all frame variants:
//!
//! - **Formula**: all frame variants share the same underlying
//!   `affn::spherical::angular_separation_impl` which uses the Vincenty
//!   formula `atan2(|n×m|, n·m)` — NOT naive `acos(n·m)`.
//! - **Return type**: typed `Degrees` (`qtty::angular::Degrees`), not raw `f64`.
//! - **Identity**: `a.angular_separation(&a) == Degrees::new(0.0)` exactly.
//! - **Antipode**: `a.angular_separation(&antipode(a)) == 180°`.
//! - **Symmetry**: `sep(a,b) ≈ sep(b,a)` within a small tolerance
//!   (the computation paths differ in FP order when swapping arguments).
//! - **Cross-frame**: angular separation is a geometric invariant; equatorial
//!   and ecliptic frames must agree to 1 × 10⁻¹² rad.
//! - **Precision near 0**: two points 10⁻⁹ rad apart must yield a result > 0
//!   and within 1 % of the expected value.  Naive `acos` would lose all
//!   significant digits here.
//! - **Precision near π**: two near-antipodal points must yield a result
//!   within 10⁻¹² rad of the expected value.

use std::f64::consts::PI;

use siderust::coordinates::spherical::direction;
use siderust::coordinates::transform::TransformFrame;
use siderust::qtty::*;

// ---------------------------------------------------------------------------
// Internal helper: Lehmer LCG for deterministic pseudo-random test angles.
// ---------------------------------------------------------------------------
fn lcg_next(state: &mut u64) -> f64 {
    *state = state
        .wrapping_mul(6_364_136_223_846_793_005)
        .wrapping_add(1_442_695_040_888_963_407);
    // Upper 31 bits → [0, 1)
    ((*state >> 33) as f64) / (0x8000_0000_u64 as f64)
}

// ---------------------------------------------------------------------------
// 1. Identity – identical directions must yield exactly 0°
// ---------------------------------------------------------------------------

#[test]
fn angular_separation_identity_is_exact_zero() {
    let cases: &[(f64, f64)] = &[
        (0.0, 0.0),
        (37.9546, 89.2641),   // near Polaris
        (279.2347, 38.7836),  // near Vega
        (180.0, -45.0),
        (359.9999, -89.9999),
    ];

    for &(ra, dec) in cases {
        let d = direction::EquatorialMeanJ2000::new(Degrees::new(ra), Degrees::new(dec));
        let sep = d.angular_separation(&d).value();
        assert_eq!(
            sep, 0.0,
            "identity failed for (RA={ra}, Dec={dec}): got {sep}"
        );
    }
}

// ---------------------------------------------------------------------------
// 2. Antipode – must return exactly 180° (π rad) within floating-point limits
// ---------------------------------------------------------------------------

#[test]
fn angular_separation_poles_are_antipodal() {
    let north = direction::EquatorialMeanJ2000::new(0.0 * DEG, 90.0 * DEG);
    let south = direction::EquatorialMeanJ2000::new(0.0 * DEG, -90.0 * DEG);
    let sep = north.angular_separation(&south).value();
    assert!(
        (sep - 180.0).abs() < 1e-12,
        "pole antipode: expected 180°, got {sep}"
    );
}

#[test]
fn angular_separation_general_antipode() {
    // For a direction at (ra, dec), its antipode is (ra + 180°, -dec).
    let cases: &[(f64, f64)] = &[
        (45.0, 30.0),
        (270.0, -60.0),
        (120.0, 0.0),
        (0.0, 45.0),
        (10.0, 89.0),  // near pole
    ];

    for &(ra, dec) in cases {
        let a = direction::EquatorialMeanJ2000::new(Degrees::new(ra), Degrees::new(dec));
        let antipode = direction::EquatorialMeanJ2000::new(
            Degrees::new((ra + 180.0) % 360.0),
            Degrees::new(-dec),
        );
        let sep = a.angular_separation(&antipode).value();
        assert!(
            (sep - 180.0).abs() < 1e-11,
            "antipode ({ra},{dec}): expected 180°, got {sep}"
        );
    }
}

// ---------------------------------------------------------------------------
// 3. Symmetry – sep(a,b) ≈ sep(b,a) over 100 seeded random pairs.
//    The Vincenty formula is mathematically symmetric; in floating-point the
//    two evaluation paths differ slightly.  We accept up to 4 ulps.
// ---------------------------------------------------------------------------

#[test]
fn angular_separation_symmetry_100_random_pairs() {
    let mut state: u64 = 0xdead_beef_cafe_f00d;

    for i in 0u64..100 {
        let ra1 = Degrees::new(lcg_next(&mut state) * 360.0);
        let dec1 = Degrees::new(lcg_next(&mut state) * 180.0 - 90.0);
        let ra2 = Degrees::new(lcg_next(&mut state) * 360.0);
        let dec2 = Degrees::new(lcg_next(&mut state) * 180.0 - 90.0);

        let a = direction::EquatorialMeanJ2000::new(ra1, dec1);
        let b = direction::EquatorialMeanJ2000::new(ra2, dec2);

        let sep_ab = a.angular_separation(&b).value();
        let sep_ba = b.angular_separation(&a).value();

        // Angular separation is always non-negative; compare raw bit patterns.
        let ba = sep_ab.to_bits();
        let bb = sep_ba.to_bits();
        let ulp_diff = if ba >= bb { ba - bb } else { bb - ba };

        assert!(
            ulp_diff <= 4,
            "pair {i}: symmetry violated: sep(a,b)={sep_ab:.15e} sep(b,a)={sep_ba:.15e} ulps={ulp_diff}"
        );
    }
}

// ---------------------------------------------------------------------------
// 4. Cross-frame agreement – Vega vs Polaris: equatorial == ecliptic
//    Angular separation is a geometric invariant under any rigid rotation.
// ---------------------------------------------------------------------------

#[test]
fn angular_separation_cross_frame_equatorial_vs_ecliptic() {
    // Vega  J2000: RA=279.2347°, Dec=+38.7836°
    // Polaris J2000: RA=37.9546°, Dec=+89.2641°
    let vega_eq =
        direction::EquatorialMeanJ2000::new(Degrees::new(279.2347), Degrees::new(38.7836));
    let polaris_eq =
        direction::EquatorialMeanJ2000::new(Degrees::new(37.9546), Degrees::new(89.2641));

    let sep_eq_deg = vega_eq.angular_separation(&polaris_eq).value();

    // Convert both directions to ecliptic mean J2000.
    let vega_ecl: direction::EclipticMeanJ2000 = TransformFrame::to_frame(&vega_eq);
    let polaris_ecl: direction::EclipticMeanJ2000 = TransformFrame::to_frame(&polaris_eq);

    let sep_ecl_deg = vega_ecl.angular_separation(&polaris_ecl).value();

    let diff_rad = (sep_eq_deg - sep_ecl_deg).abs().to_radians();
    assert!(
        diff_rad < 1e-12,
        "cross-frame disagreement: equatorial={sep_eq_deg:.15}° ecliptic={sep_ecl_deg:.15}° diff={diff_rad:.2e} rad"
    );
}

// ---------------------------------------------------------------------------
// 5. Precision near 0 – two points 1×10⁻⁹ rad apart.
//    Naive acos(dot) loses all significant figures here; Vincenty does not.
// ---------------------------------------------------------------------------

#[test]
fn angular_separation_precision_near_zero() {
    let epsilon_rad = 1e-9_f64;
    let epsilon_deg = epsilon_rad.to_degrees(); // ≈ 5.73×10⁻⁸ °

    // Separation entirely in RA along the equator.
    let a = direction::EquatorialMeanJ2000::new(0.0 * DEG, 0.0 * DEG);
    let b = direction::EquatorialMeanJ2000::new(Degrees::new(epsilon_deg), 0.0 * DEG);

    let sep_deg = a.angular_separation(&b).value();

    assert!(
        sep_deg > 0.0,
        "precision near 0: Vincenty returned exactly 0 (expected ~{epsilon_deg:.2e}°)"
    );
    let rel_err = (sep_deg - epsilon_deg).abs() / epsilon_deg;
    assert!(
        rel_err < 0.01,
        "precision near 0: relative error {rel_err:.2e} (expected < 1 %)"
    );
}

// ---------------------------------------------------------------------------
// 6. Precision near π – two near-antipodal points separated by (π − 10⁻⁹) rad.
//    Vincenty handles this without precision loss.
// ---------------------------------------------------------------------------

#[test]
fn angular_separation_precision_near_pi() {
    let delta_rad = 1e-9_f64; // deviation from exact antipode
    let delta_deg = delta_rad.to_degrees();

    // North pole, and a point delta_rad above the south pole (dec = -(90 - δ)°).
    let north = direction::EquatorialMeanJ2000::new(0.0 * DEG, 90.0 * DEG);
    let near_south =
        direction::EquatorialMeanJ2000::new(0.0 * DEG, Degrees::new(-(90.0 - delta_deg)));

    let sep_rad = north.angular_separation(&near_south).value().to_radians();
    let expected_rad = PI - delta_rad;
    let error_rad = (sep_rad - expected_rad).abs();

    assert!(
        error_rad < 1e-12,
        "precision near pi: error={error_rad:.2e} rad (expected {expected_rad:.15}, got {sep_rad:.15})"
    );
}

// ---------------------------------------------------------------------------
// 7. Galactic and Horizontal frames share the same implementation.
//    Quick sanity checks to confirm the formula applies there too.
// ---------------------------------------------------------------------------

#[test]
fn angular_separation_galactic_identity_and_quarter_turn() {
    // Identity
    let gc = direction::Galactic::new(120.0 * DEG, 45.0 * DEG);
    assert_eq!(gc.angular_separation(&gc).value(), 0.0, "galactic identity");

    // Two galactic directions 90° apart in longitude at b=0
    let a = direction::Galactic::new(0.0 * DEG, 0.0 * DEG);
    let b = direction::Galactic::new(90.0 * DEG, 0.0 * DEG);
    let sep = a.angular_separation(&b).value();
    assert!(
        (sep - 90.0).abs() < 1e-12,
        "galactic 90° separation: expected 90°, got {sep}"
    );
}

#[test]
fn angular_separation_horizontal_identity_and_quarter_turn() {
    // Identity
    let h = direction::Horizontal::new(180.0 * DEG, 30.0 * DEG);
    assert_eq!(h.angular_separation(&h).value(), 0.0, "horizontal identity");

    // Zenith (alt=90) and horizon (alt=0) are 90° apart
    let zenith = direction::Horizontal::new(0.0 * DEG, 90.0 * DEG);
    let horizon = direction::Horizontal::new(0.0 * DEG, 0.0 * DEG);
    let sep = zenith.angular_separation(&horizon).value();
    assert!(
        (sep - 90.0).abs() < 1e-12,
        "horizontal 90° separation: expected 90°, got {sep}"
    );
}
