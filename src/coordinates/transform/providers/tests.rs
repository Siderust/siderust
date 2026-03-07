// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

use super::*;

const EPSILON: f64 = 1e-10;

#[test]
fn test_identity_frame_rotation() {
    let rot = frame_rotation::<ICRS, ICRS>(JulianDate::J2000, &AstroContext::default());
    let v = [1.0, 2.0, 3.0];
    let result = rot.apply_array(v);
    assert!((result[0] - v[0]).abs() < EPSILON);
    assert!((result[1] - v[1]).abs() < EPSILON);
    assert!((result[2] - v[2]).abs() < EPSILON);
}

#[test]
fn test_icrs_to_ecliptic_rotation() {
    let rot =
        frame_rotation::<ICRS, EclipticMeanJ2000>(JulianDate::J2000, &AstroContext::default());
    let v = [1.0, 2.0, 3.0];
    let w = rot.apply_array(v);

    assert!(w[0].is_finite() && w[1].is_finite() && w[2].is_finite());

    let nv = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    let nw = (w[0] * w[0] + w[1] * w[1] + w[2] * w[2]).sqrt();
    assert!((nv - nw).abs() < 1e-12);
}

#[test]
fn test_ecliptic_icrs_roundtrip() {
    let ctx = AstroContext::default();
    let jd = JulianDate::J2000;

    let r1 = frame_rotation::<ICRS, EclipticMeanJ2000>(jd, &ctx);
    let r2 = frame_rotation::<EclipticMeanJ2000, ICRS>(jd, &ctx);

    let v = [1.0, 2.0, 3.0];
    let roundtrip = r2.apply_array(r1.apply_array(v));

    assert!((roundtrip[0] - v[0]).abs() < EPSILON);
    assert!((roundtrip[1] - v[1]).abs() < EPSILON);
    assert!((roundtrip[2] - v[2]).abs() < EPSILON);
}

#[test]
fn test_identity_center_shift() {
    let shift = center_shift::<Barycentric, Barycentric, EclipticMeanJ2000>(
        JulianDate::J2000,
        &AstroContext::default(),
    );
    assert!((shift[0]).abs() < EPSILON);
    assert!((shift[1]).abs() < EPSILON);
    assert!((shift[2]).abs() < EPSILON);
}

#[test]
fn test_helio_bary_geo_composition() {
    let ctx = AstroContext::default();
    let jd = JulianDate::J2000;

    let helio_geo = center_shift::<Heliocentric, Geocentric, EclipticMeanJ2000>(jd, &ctx);
    let helio_bary = center_shift::<Heliocentric, Barycentric, EclipticMeanJ2000>(jd, &ctx);
    let bary_geo = center_shift::<Barycentric, Geocentric, EclipticMeanJ2000>(jd, &ctx);

    let composed = [
        helio_bary[0] + bary_geo[0],
        helio_bary[1] + bary_geo[1],
        helio_bary[2] + bary_geo[2],
    ];

    assert!((helio_geo[0] - composed[0]).abs() < EPSILON);
    assert!((helio_geo[1] - composed[1]).abs() < EPSILON);
    assert!((helio_geo[2] - composed[2]).abs() < EPSILON);
}

#[test]
fn test_center_shift_antisymmetry() {
    let ctx = AstroContext::default();
    let jd = JulianDate::J2000;

    let forward = center_shift::<Heliocentric, Geocentric, EclipticMeanJ2000>(jd, &ctx);
    let backward = center_shift::<Geocentric, Heliocentric, EclipticMeanJ2000>(jd, &ctx);

    assert!((forward[0] + backward[0]).abs() < EPSILON);
    assert!((forward[1] + backward[1]).abs() < EPSILON);
    assert!((forward[2] + backward[2]).abs() < EPSILON);
}

#[test]
fn test_frame_bias_is_non_identity() {
    let ctx = AstroContext::default();
    let jd = JulianDate::J2000;

    let rot = frame_rotation::<ICRS, EquatorialMeanJ2000>(jd, &ctx);
    let v = [0.0, 1.0, 0.0];
    let out = rot.apply_array(v);

    let delta = (out[0] - v[0]).abs() + (out[1] - v[1]).abs() + (out[2] - v[2]).abs();
    assert!(delta > 1e-12, "frame bias should not be identity");
}

#[test]
fn test_icrs_ecliptic_roundtrip_is_identity() {
    let ctx = AstroContext::default();
    let jd = JulianDate::J2000;

    let r = frame_rotation::<ICRS, EclipticMeanJ2000>(jd, &ctx);
    let rinv = frame_rotation::<EclipticMeanJ2000, ICRS>(jd, &ctx);

    let v = [1.0, 0.0, 0.0];
    let round = rinv.apply_array(r.apply_array(v));
    let err = (round[0] - v[0]).abs() + (round[1] - v[1]).abs() + (round[2] - v[2]).abs();
    assert!(
        err < 1e-12,
        "ICRS↔EclipticMeanJ2000 roundtrip should be identity"
    );
}

#[test]
fn test_precession_identity_at_j2000() {
    let ctx = AstroContext::default();
    let jd = JulianDate::J2000;

    let rot = frame_rotation::<EquatorialMeanJ2000, EquatorialMeanOfDate>(jd, &ctx);
    let v = [1.0, 0.0, 0.0];
    let out = rot.apply_array(v);

    assert!((out[0] - v[0]).abs() < 1e-6);
    assert!((out[1] - v[1]).abs() < 1e-6);
    assert!((out[2] - v[2]).abs() < 1e-6);
}

#[test]
fn test_nutation_rotation_roundtrip() {
    let ctx = AstroContext::default();
    let jd = JulianDate::new(2_460_000.5);

    let rot = frame_rotation::<EquatorialMeanOfDate, EquatorialTrueOfDate>(jd, &ctx);
    let inv = frame_rotation::<EquatorialTrueOfDate, EquatorialMeanOfDate>(jd, &ctx);

    let v = [0.3, 0.4, 0.5];
    let roundtrip = inv.apply_array(rot.apply_array(v));

    assert!((roundtrip[0] - v[0]).abs() < 1e-12);
    assert!((roundtrip[1] - v[1]).abs() < 1e-12);
    assert!((roundtrip[2] - v[2]).abs() < 1e-12);
}

#[test]
fn test_icrf_icrs_identity_rotation() {
    let ctx = AstroContext::default();
    let jd = JulianDate::J2000;
    let rot = frame_rotation::<ICRF, ICRS>(jd, &ctx);
    let v = [0.1, -0.2, 0.3];
    let out = rot.apply_array(v);
    assert!((out[0] - v[0]).abs() < 1e-15);
    assert!((out[1] - v[1]).abs() < 1e-15);
    assert!((out[2] - v[2]).abs() < 1e-15);
}

#[test]
fn test_icrf_ecliptic_roundtrip() {
    let ctx = AstroContext::default();
    let jd = JulianDate::J2000;
    let r = frame_rotation::<ICRF, EclipticMeanJ2000>(jd, &ctx);
    let rinv = frame_rotation::<EclipticMeanJ2000, ICRF>(jd, &ctx);

    let v = [1.0, 2.0, 3.0];
    let round = rinv.apply_array(r.apply_array(v));
    assert!((round[0] - v[0]).abs() < 1e-12);
    assert!((round[1] - v[1]).abs() < 1e-12);
    assert!((round[2] - v[2]).abs() < 1e-12);
}

#[test]
fn test_icrf_to_ecliptic_matches_icrs_to_ecliptic() {
    let ctx = AstroContext::default();
    let jd = JulianDate::J2000;
    let via_icrf = frame_rotation::<ICRF, EclipticMeanJ2000>(jd, &ctx);
    let via_icrs = frame_rotation::<ICRS, EclipticMeanJ2000>(jd, &ctx);

    let v = [0.3, -0.1, 0.8];
    let a = via_icrf.apply_array(v);
    let b = via_icrs.apply_array(v);
    assert!((a[0] - b[0]).abs() < 1e-15);
    assert!((a[1] - b[1]).abs() < 1e-15);
    assert!((a[2] - b[2]).abs() < 1e-15);
}
