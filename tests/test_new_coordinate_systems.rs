// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Tests for the extended coordinate systems: FK4 B1950, TEME, Galactic,
//! planetary body-fixed frames, and planetocentric centers.

use qtty::*;
use siderust::coordinates::{
    cartesian,
    centers::*,
    frames::*,
    spherical,
    transform::{
        providers::{center_shift, frame_rotation},
        AstroContext,
    },
};
use siderust::time::JulianDate;

const EPSILON: f64 = 1e-10;
const LOOSE_EPS: f64 = 1e-6;

// =============================================================================
// Frame Rotation Tests: FK4 B1950
// =============================================================================

#[test]
fn fk4_b1950_to_icrs_roundtrip() {
    let ctx = AstroContext::default();
    let jd = JulianDate::J2000;

    let r = frame_rotation::<FK4B1950, ICRS>(jd, &ctx);
    let rinv = frame_rotation::<ICRS, FK4B1950>(jd, &ctx);

    let v = [0.3, -0.5, 0.7];
    let roundtrip = rinv.apply_array(r.apply_array(v));

    assert!((roundtrip[0] - v[0]).abs() < EPSILON);
    assert!((roundtrip[1] - v[1]).abs() < EPSILON);
    assert!((roundtrip[2] - v[2]).abs() < EPSILON);
}

#[test]
fn fk4_b1950_is_not_identity() {
    let ctx = AstroContext::default();
    let jd = JulianDate::J2000;

    let r = frame_rotation::<FK4B1950, ICRS>(jd, &ctx);
    let v = [1.0, 0.0, 0.0];
    let out = r.apply_array(v);

    // FK4→ICRS involves a significant frame rotation (degrees-level)
    let delta = (out[0] - v[0]).abs() + (out[1] - v[1]).abs() + (out[2] - v[2]).abs();
    assert!(delta > 1e-6, "FK4→ICRS should not be identity");
}

#[test]
fn fk4_b1950_preserves_vector_length() {
    let ctx = AstroContext::default();
    let jd = JulianDate::J2000;

    let r = frame_rotation::<FK4B1950, ICRS>(jd, &ctx);
    let v = [1.0, 2.0, 3.0];
    let w = r.apply_array(v);

    let nv = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    let nw = (w[0] * w[0] + w[1] * w[1] + w[2] * w[2]).sqrt();
    assert!((nv - nw).abs() < 1e-12);
}

#[test]
fn fk4_via_equatorial_mean_j2000_roundtrip() {
    let ctx = AstroContext::default();
    let jd = JulianDate::J2000;

    let r1 = frame_rotation::<FK4B1950, ICRS>(jd, &ctx);
    let r2 = frame_rotation::<ICRS, EquatorialMeanJ2000>(jd, &ctx);
    let r3 = frame_rotation::<FK4B1950, EquatorialMeanJ2000>(jd, &ctx);

    let v = [0.5, 0.3, 0.1];
    let via_icrs = r2.apply_array(r1.apply_array(v));
    let direct = r3.apply_array(v);

    assert!((via_icrs[0] - direct[0]).abs() < EPSILON);
    assert!((via_icrs[1] - direct[1]).abs() < EPSILON);
    assert!((via_icrs[2] - direct[2]).abs() < EPSILON);
}

// =============================================================================
// Frame Rotation Tests: TEME
// =============================================================================

#[test]
fn teme_to_tod_roundtrip() {
    let ctx = AstroContext::default();
    let jd = JulianDate::new(2_460_000.5);

    let r = frame_rotation::<TEME, EquatorialTrueOfDate>(jd, &ctx);
    let rinv = frame_rotation::<EquatorialTrueOfDate, TEME>(jd, &ctx);

    let v = [0.4, -0.6, 0.2];
    let roundtrip = rinv.apply_array(r.apply_array(v));

    assert!((roundtrip[0] - v[0]).abs() < EPSILON);
    assert!((roundtrip[1] - v[1]).abs() < EPSILON);
    assert!((roundtrip[2] - v[2]).abs() < EPSILON);
}

#[test]
fn teme_to_icrs_roundtrip() {
    let ctx = AstroContext::default();
    let jd = JulianDate::new(2_460_000.5);

    let r = frame_rotation::<TEME, ICRS>(jd, &ctx);
    let rinv = frame_rotation::<ICRS, TEME>(jd, &ctx);

    let v = [1.0, 0.0, 0.0];
    let roundtrip = rinv.apply_array(r.apply_array(v));

    assert!((roundtrip[0] - v[0]).abs() < EPSILON);
    assert!((roundtrip[1] - v[1]).abs() < EPSILON);
    assert!((roundtrip[2] - v[2]).abs() < EPSILON);
}

#[test]
fn teme_to_tod_is_z_rotation() {
    // TEME → TOD should be a pure z-axis rotation (equation of equinoxes).
    let ctx = AstroContext::default();
    let jd = JulianDate::new(2_460_000.5);

    let r = frame_rotation::<TEME, EquatorialTrueOfDate>(jd, &ctx);
    let v = [0.0, 0.0, 1.0]; // z-axis should be invariant under z-rotation
    let out = r.apply_array(v);

    assert!((out[0]).abs() < EPSILON, "z-axis rotated in x: {}", out[0]);
    assert!((out[1]).abs() < EPSILON, "z-axis rotated in y: {}", out[1]);
    assert!((out[2] - 1.0).abs() < EPSILON, "z-axis changed: {}", out[2]);
}

#[test]
fn teme_preserves_vector_length() {
    let ctx = AstroContext::default();
    let jd = JulianDate::new(2_460_000.5);

    let r = frame_rotation::<TEME, ICRS>(jd, &ctx);
    let v = [1.0, 2.0, 3.0];
    let w = r.apply_array(v);

    let nv = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    let nw = (w[0] * w[0] + w[1] * w[1] + w[2] * w[2]).sqrt();
    assert!((nv - nw).abs() < 1e-12);
}

// =============================================================================
// Frame Rotation Tests: Galactic
// =============================================================================

#[test]
fn galactic_to_icrs_roundtrip() {
    let ctx = AstroContext::default();
    let jd = JulianDate::J2000;

    let r = frame_rotation::<Galactic, ICRS>(jd, &ctx);
    let rinv = frame_rotation::<ICRS, Galactic>(jd, &ctx);

    let v = [0.1, 0.9, -0.3];
    let roundtrip = rinv.apply_array(r.apply_array(v));

    assert!((roundtrip[0] - v[0]).abs() < EPSILON);
    assert!((roundtrip[1] - v[1]).abs() < EPSILON);
    assert!((roundtrip[2] - v[2]).abs() < EPSILON);
}

#[test]
fn galactic_is_not_identity() {
    let ctx = AstroContext::default();
    let jd = JulianDate::J2000;

    let r = frame_rotation::<Galactic, ICRS>(jd, &ctx);
    let v = [1.0, 0.0, 0.0];
    let out = r.apply_array(v);

    // Galactic→ICRS involves a significant frame rotation (~60°)
    let delta = (out[0] - v[0]).abs() + (out[1] - v[1]).abs() + (out[2] - v[2]).abs();
    assert!(delta > 0.1, "Galactic→ICRS should be a major rotation");
}

#[test]
fn galactic_preserves_vector_length() {
    let ctx = AstroContext::default();
    let jd = JulianDate::J2000;

    let r = frame_rotation::<Galactic, ICRS>(jd, &ctx);
    let v = [1.0, 2.0, 3.0];
    let w = r.apply_array(v);

    let nv = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    let nw = (w[0] * w[0] + w[1] * w[1] + w[2] * w[2]).sqrt();
    // Hipparcos AG matrix has ~6e-11 non-orthogonality at f64 precision,
    // amplified by ||v|| ≈ 3.74 → ~2e-10 norm deviation.
    assert!((nv - nw).abs() < 1e-9);
}

#[test]
fn galactic_via_equatorial_j2000_roundtrip() {
    let ctx = AstroContext::default();
    let jd = JulianDate::J2000;

    // Galactic→ICRS→EquatorialMeanJ2000 composed vs direct
    let r1 = frame_rotation::<Galactic, ICRS>(jd, &ctx);
    let r2 = frame_rotation::<ICRS, EquatorialMeanJ2000>(jd, &ctx);
    let r3 = frame_rotation::<Galactic, EquatorialMeanJ2000>(jd, &ctx);

    let v = [0.5, 0.3, 0.1];
    let via_icrs = r2.apply_array(r1.apply_array(v));
    let direct = r3.apply_array(v);

    assert!((via_icrs[0] - direct[0]).abs() < EPSILON);
    assert!((via_icrs[1] - direct[1]).abs() < EPSILON);
    assert!((via_icrs[2] - direct[2]).abs() < EPSILON);
}

// =============================================================================
// Frame Rotation Tests: GCRS ↔ ICRS (near-identity)
// =============================================================================

#[test]
fn gcrs_to_icrs_is_near_identity() {
    let ctx = AstroContext::default();
    let jd = JulianDate::J2000;

    let r = frame_rotation::<GCRS, ICRS>(jd, &ctx);
    let v = [1.0, 0.0, 0.0];
    let out = r.apply_array(v);

    // GCRS ≈ ICRS (< 1 mas)
    assert!((out[0] - v[0]).abs() < LOOSE_EPS);
    assert!((out[1] - v[1]).abs() < LOOSE_EPS);
    assert!((out[2] - v[2]).abs() < LOOSE_EPS);
}

// =============================================================================
// Frame Rotation Tests: Planetary Body-Fixed
// =============================================================================

#[test]
fn mars_fixed_to_icrs_roundtrip() {
    let ctx = AstroContext::default();
    let jd = JulianDate::J2000;

    let r = frame_rotation::<planetary::MarsFixed, ICRS>(jd, &ctx);
    let rinv = frame_rotation::<ICRS, planetary::MarsFixed>(jd, &ctx);

    let v = [1.0, 0.0, 0.0];
    let roundtrip = rinv.apply_array(r.apply_array(v));

    assert!((roundtrip[0] - v[0]).abs() < EPSILON);
    assert!((roundtrip[1] - v[1]).abs() < EPSILON);
    assert!((roundtrip[2] - v[2]).abs() < EPSILON);
}

#[test]
fn mars_fixed_preserves_length() {
    let ctx = AstroContext::default();
    let jd = JulianDate::new(2_460_000.5);

    let r = frame_rotation::<planetary::MarsFixed, ICRS>(jd, &ctx);
    let v = [1.0, 2.0, 3.0];
    let w = r.apply_array(v);

    let nv = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    let nw = (w[0] * w[0] + w[1] * w[1] + w[2] * w[2]).sqrt();
    assert!((nv - nw).abs() < 1e-12);
}

#[test]
fn all_body_fixed_roundtrips() {
    let ctx = AstroContext::default();
    let jd = JulianDate::new(2_460_000.5);
    let v = [0.3, -0.7, 0.5];

    // Test all body-fixed frames roundtrip through ICRS
    macro_rules! test_roundtrip {
        ($frame:ty) => {{
            let r = frame_rotation::<$frame, ICRS>(jd, &ctx);
            let rinv = frame_rotation::<ICRS, $frame>(jd, &ctx);
            let rt = rinv.apply_array(r.apply_array(v));
            assert!(
                (rt[0] - v[0]).abs() < EPSILON
                    && (rt[1] - v[1]).abs() < EPSILON
                    && (rt[2] - v[2]).abs() < EPSILON,
                "{} roundtrip failed: {:?} != {:?}",
                stringify!($frame),
                rt,
                v
            );
        }};
    }

    test_roundtrip!(planetary::MercuryFixed);
    test_roundtrip!(planetary::VenusFixed);
    test_roundtrip!(planetary::MarsFixed);
    test_roundtrip!(planetary::MoonPrincipalAxes);
    test_roundtrip!(planetary::JupiterSystemIII);
    test_roundtrip!(planetary::SaturnFixed);
    test_roundtrip!(planetary::UranusFixed);
    test_roundtrip!(planetary::NeptuneFixed);
    test_roundtrip!(planetary::PlutoFixed);
}

#[test]
fn body_fixed_rotation_varies_with_time() {
    // Body-fixed frames rotate with the body, so the rotation should
    // differ at different epochs.
    let ctx = AstroContext::default();
    let jd1 = JulianDate::J2000;
    let jd2 = JulianDate::new(2_460_000.5); // ~2023

    let r1 = frame_rotation::<planetary::MarsFixed, ICRS>(jd1, &ctx);
    let r2 = frame_rotation::<planetary::MarsFixed, ICRS>(jd2, &ctx);

    let v = [1.0, 0.0, 0.0];
    let w1 = r1.apply_array(v);
    let w2 = r2.apply_array(v);

    let diff = (w1[0] - w2[0]).abs() + (w1[1] - w2[1]).abs() + (w1[2] - w2[2]).abs();
    assert!(diff > 0.01, "body-fixed rotation should vary with time");
}

// =============================================================================
// Center Shift Tests: Planetocentric
// =============================================================================

#[test]
fn planetocentric_to_bary_antisymmetry() {
    let ctx = AstroContext::default();
    let jd = JulianDate::J2000;

    macro_rules! test_antisymmetry {
        ($center:ty) => {{
            let fwd = center_shift::<$center, Barycentric, EclipticMeanJ2000>(jd, &ctx);
            let bwd = center_shift::<Barycentric, $center, EclipticMeanJ2000>(jd, &ctx);
            assert!(
                (fwd[0] + bwd[0]).abs() < EPSILON
                    && (fwd[1] + bwd[1]).abs() < EPSILON
                    && (fwd[2] + bwd[2]).abs() < EPSILON,
                "{} antisymmetry failed",
                stringify!($center)
            );
        }};
    }

    test_antisymmetry!(Mercurycentric);
    test_antisymmetry!(Venuscentric);
    test_antisymmetry!(Marscentric);
    test_antisymmetry!(Jovicentric);
    test_antisymmetry!(Saturnocentric);
    test_antisymmetry!(Uranocentric);
    test_antisymmetry!(Neptunocentric);
    test_antisymmetry!(Plutocentric);
    test_antisymmetry!(Selenocentric);
}

#[test]
fn selenocentric_to_bary_antisymmetry() {
    let ctx = AstroContext::default();
    let jd = JulianDate::J2000;

    let fwd = center_shift::<Selenocentric, Barycentric, EclipticMeanJ2000>(jd, &ctx);
    let bwd = center_shift::<Barycentric, Selenocentric, EclipticMeanJ2000>(jd, &ctx);

    assert!((fwd[0] + bwd[0]).abs() < EPSILON);
    assert!((fwd[1] + bwd[1]).abs() < EPSILON);
    assert!((fwd[2] + bwd[2]).abs() < EPSILON);
}

#[test]
fn planetocentric_center_shifts_are_nonzero() {
    let ctx = AstroContext::default();
    let jd = JulianDate::J2000;

    macro_rules! test_nonzero {
        ($center:ty) => {{
            let s = center_shift::<$center, Barycentric, EclipticMeanJ2000>(jd, &ctx);
            let norm = (s[0] * s[0] + s[1] * s[1] + s[2] * s[2]).sqrt();
            assert!(
                norm > 0.01,
                "{} center shift should be nonzero, got {}",
                stringify!($center),
                norm
            );
        }};
    }

    test_nonzero!(Mercurycentric);
    test_nonzero!(Venuscentric);
    test_nonzero!(Marscentric);
    test_nonzero!(Jovicentric);
    test_nonzero!(Saturnocentric);
    test_nonzero!(Uranocentric);
    test_nonzero!(Neptunocentric);
    test_nonzero!(Plutocentric);
    test_nonzero!(Selenocentric);
}

#[test]
fn marscentric_helio_geo_composition() {
    let ctx = AstroContext::default();
    let jd = JulianDate::J2000;

    // Marscentric → Geocentric should compose via Barycentric
    let mars_bary = center_shift::<Marscentric, Barycentric, EclipticMeanJ2000>(jd, &ctx);
    let bary_geo = center_shift::<Barycentric, Geocentric, EclipticMeanJ2000>(jd, &ctx);
    let mars_geo = center_shift::<Marscentric, Geocentric, EclipticMeanJ2000>(jd, &ctx);

    let composed = [
        mars_bary[0] + bary_geo[0],
        mars_bary[1] + bary_geo[1],
        mars_bary[2] + bary_geo[2],
    ];

    assert!((mars_geo[0] - composed[0]).abs() < EPSILON);
    assert!((mars_geo[1] - composed[1]).abs() < EPSILON);
    assert!((mars_geo[2] - composed[2]).abs() < EPSILON);
}

#[test]
fn selenocentric_geocentric_composition() {
    let ctx = AstroContext::default();
    let jd = JulianDate::J2000;

    // Selenocentric → Barycentric should equal
    // Selenocentric → Geocentric + Geocentric → Barycentric
    let sel_bary = center_shift::<Selenocentric, Barycentric, EclipticMeanJ2000>(jd, &ctx);
    let sel_geo = center_shift::<Selenocentric, Geocentric, EclipticMeanJ2000>(jd, &ctx);
    let geo_bary = center_shift::<Geocentric, Barycentric, EclipticMeanJ2000>(jd, &ctx);

    let composed = [
        sel_geo[0] + geo_bary[0],
        sel_geo[1] + geo_bary[1],
        sel_geo[2] + geo_bary[2],
    ];

    assert!(
        (sel_bary[0] - composed[0]).abs() < EPSILON,
        "x: {} vs {}",
        sel_bary[0],
        composed[0]
    );
    assert!(
        (sel_bary[1] - composed[1]).abs() < EPSILON,
        "y: {} vs {}",
        sel_bary[1],
        composed[1]
    );
    assert!(
        (sel_bary[2] - composed[2]).abs() < EPSILON,
        "z: {} vs {}",
        sel_bary[2],
        composed[2]
    );
}

#[test]
fn selenocentric_geocentric_distance_is_reasonable() {
    let ctx = AstroContext::default();
    let jd = JulianDate::J2000;

    let s = center_shift::<Geocentric, Selenocentric, EclipticMeanJ2000>(jd, &ctx);
    let dist_au = (s[0] * s[0] + s[1] * s[1] + s[2] * s[2]).sqrt();

    // Moon is ~0.00257 AU from Earth
    assert!(
        dist_au > 0.001 && dist_au < 0.005,
        "Moon distance should be ~0.0026 AU, got {} AU",
        dist_au
    );
}

#[test]
fn jupiter_distance_is_reasonable() {
    let ctx = AstroContext::default();
    let jd = JulianDate::J2000;

    let s = center_shift::<Heliocentric, Jovicentric, EclipticMeanJ2000>(jd, &ctx);
    let dist_au = (s[0] * s[0] + s[1] * s[1] + s[2] * s[2]).sqrt();

    // Jupiter is ~5.2 AU from the Sun
    assert!(
        dist_au > 4.0 && dist_au < 6.5,
        "Jupiter distance should be ~5.2 AU, got {} AU",
        dist_au
    );
}

// =============================================================================
// Type Alias Construction Tests
// =============================================================================

#[test]
fn can_construct_fk4_direction() {
    // FK4B1950 uses new_raw(polar=dec, azimuth=ra)
    let dir = spherical::direction::FK4B1950::new_raw(Degrees::new(45.0), Degrees::new(180.0));
    assert!((dir.polar.value() - 45.0).abs() < LOOSE_EPS);
    assert!((dir.azimuth.value() - 180.0).abs() < LOOSE_EPS);
}

#[test]
fn can_construct_teme_direction() {
    let dir = spherical::direction::TEME::new_raw(Degrees::new(-30.0), Degrees::new(90.0));
    assert!((dir.polar.value() - (-30.0)).abs() < LOOSE_EPS);
    assert!((dir.azimuth.value() - 90.0).abs() < LOOSE_EPS);
}

#[test]
fn can_construct_galactic_position() {
    let pos = spherical::position::Galactic::<AstronomicalUnit>::new(
        Degrees::new(120.0),
        Degrees::new(30.0),
        10.0,
    );
    assert!((pos.distance.value() - 10.0).abs() < LOOSE_EPS);
}

#[test]
fn can_construct_marsfixed_position() {
    let pos = cartesian::position::MarsFixed::<AstronomicalUnit>::new(1.0 * AU, 0.5 * AU, 0.2 * AU);
    assert!((pos.x().value() - 1.0).abs() < LOOSE_EPS);
}

#[test]
fn can_construct_moon_principal_axes_position() {
    let pos = cartesian::position::MoonPrincipalAxes::<AstronomicalUnit>::new(
        0.001 * AU,
        0.002 * AU,
        0.0005 * AU,
    );
    assert!((pos.x().value() - 0.001).abs() < LOOSE_EPS);
}

// =============================================================================
// IAU Rotation Parameters Tests
// =============================================================================

#[test]
fn mars_rotation_parameters_at_j2000() {
    use siderust::coordinates::frames::planetary::MARS_ROTATION;

    // At T=0 (J2000.0), IAU 2015: α₀ = 317.269, δ₀ = 54.432, W(d=0) = 176.049
    let alpha0 = MARS_ROTATION.alpha0(JulianDate::J2000);
    let delta0 = MARS_ROTATION.delta0(JulianDate::J2000);
    let w = MARS_ROTATION.w(JulianDate::J2000);

    assert!(
        (alpha0.value() - 317.269).abs() < 0.001,
        "Mars α₀ at J2000: {}",
        alpha0.value()
    );
    assert!(
        (delta0.value() - 54.432).abs() < 0.001,
        "Mars δ₀ at J2000: {}",
        delta0.value()
    );
    assert!(
        (w.value() - 176.049).abs() < 0.001,
        "Mars W at J2000: {}",
        w.value()
    );
}

#[test]
fn jupiter_rotation_parameters_at_j2000() {
    use siderust::coordinates::frames::planetary::JUPITER_ROTATION;

    let alpha0 = JUPITER_ROTATION.alpha0(JulianDate::J2000);
    let delta0 = JUPITER_ROTATION.delta0(JulianDate::J2000);
    let w = JUPITER_ROTATION.w(JulianDate::J2000);

    // α₀ = 268.056595, δ₀ = 64.495303
    assert!(
        (alpha0.value() - 268.056595).abs() < 0.001,
        "Jupiter α₀: {}",
        alpha0.value()
    );
    assert!(
        (delta0.value() - 64.495303).abs() < 0.001,
        "Jupiter δ₀: {}",
        delta0.value()
    );
    assert!(w.value().is_finite(), "Jupiter W should be finite");
}

// =============================================================================
// TransformFrame trait tests with new frames
// =============================================================================

#[test]
fn direction_transform_fk4_to_icrs() {
    // Create a direction in FK4, rotate to ICRS, and verify it's still unit
    let ctx = AstroContext::default();
    let jd = JulianDate::J2000;

    let r = frame_rotation::<FK4B1950, ICRS>(jd, &ctx);
    let v = [1.0, 0.0, 0.0];
    let w = r.apply_array(v);

    let norm = (w[0] * w[0] + w[1] * w[1] + w[2] * w[2]).sqrt();
    assert!((norm - 1.0).abs() < 1e-12, "should remain unit vector");
}

#[test]
fn direction_transform_galactic_to_icrs() {
    let ctx = AstroContext::default();
    let jd = JulianDate::J2000;

    // Galactic z-axis = North Galactic Pole direction in ICRS
    let r = frame_rotation::<Galactic, ICRS>(jd, &ctx);
    let ngp = r.apply_array([0.0, 0.0, 1.0]);

    let norm = (ngp[0] * ngp[0] + ngp[1] * ngp[1] + ngp[2] * ngp[2]).sqrt();
    // Hipparcos AG matrix column 2 norm deviates by ~6e-11 at f64 precision.
    assert!((norm - 1.0).abs() < 1e-10, "should remain unit vector");
    // The z component should be positive (NGP has positive declination ~27°)
    assert!(ngp[2] > 0.0, "NGP should have positive z in ICRS");
}

// =============================================================================
// Cross-frame composition tests
// =============================================================================

#[test]
fn galactic_to_ecliptic_roundtrip() {
    let ctx = AstroContext::default();
    let jd = JulianDate::J2000;

    let r = frame_rotation::<Galactic, EclipticMeanJ2000>(jd, &ctx);
    let rinv = frame_rotation::<EclipticMeanJ2000, Galactic>(jd, &ctx);

    let v = [0.5, -0.3, 0.8];
    let rt = rinv.apply_array(r.apply_array(v));

    // Hipparcos AG matrix has ~6e-11 non-orthogonality: R^T·R deviates
    // from I by ~4e-11 in off-diagonal elements. Roundtrip error ~2e-10.
    let gal_eps = 1e-9;
    assert!((rt[0] - v[0]).abs() < gal_eps);
    assert!((rt[1] - v[1]).abs() < gal_eps);
    assert!((rt[2] - v[2]).abs() < gal_eps);
}

#[test]
fn gcrs_to_ecliptic_matches_icrs_to_ecliptic() {
    let ctx = AstroContext::default();
    let jd = JulianDate::J2000;

    let via_gcrs = frame_rotation::<GCRS, EclipticMeanJ2000>(jd, &ctx);
    let via_icrs = frame_rotation::<ICRS, EclipticMeanJ2000>(jd, &ctx);

    // GCRS ≈ ICRS, so these should be nearly identical
    let v = [0.5, 0.3, 0.1];
    let a = via_gcrs.apply_array(v);
    let b = via_icrs.apply_array(v);

    assert!((a[0] - b[0]).abs() < LOOSE_EPS);
    assert!((a[1] - b[1]).abs() < LOOSE_EPS);
    assert!((a[2] - b[2]).abs() < LOOSE_EPS);
}
