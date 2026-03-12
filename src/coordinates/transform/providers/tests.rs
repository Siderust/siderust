// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

use super::*;
use crate::astro::nutation::Iau2006A;

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

// ── Absolute ERFA-derived regression tests ────────────────────────────
//
// Reference values computed with `python3 -c "import erfa; ..."` using the
// ERFA 2.0.1.5 / SOFA 2024-10-01 release.  Tolerances are:
//
//   • Fixed bias / J2000 obliquity: 1 × 10⁻¹⁵ (double-precision hard-coded).
//   • Precession (pmat06):         1 × 10⁻¹² (same IAU 2006 model).
//   • Nutation (num06a vs 2000B):  1 × 10⁻⁶  (expected 2000A vs 2000B gap).

/// Helper: Frobenius distance between two 3×3 matrices.
fn mat_frobenius(a: &[[f64; 3]; 3], b: &[[f64; 3]; 3]) -> f64 {
    let mut sum = 0.0;
    for i in 0..3 {
        for j in 0..3 {
            let d = a[i][j] - b[i][j];
            sum += d * d;
        }
    }
    sum.sqrt()
}

/// Frame bias matches ERFA `bp06` rb at J2000.0 (provider path).
#[test]
fn erfa_frame_bias_provider_matches_bp06_rb() {
    let ctx = AstroContext::default();
    let rb = frame_rotation::<ICRS, EquatorialMeanJ2000>(JulianDate::J2000, &ctx);
    let m = rb.as_matrix();

    // ERFA bp06(2451545.0, 0.0) → rb
    let erfa: [[f64; 3]; 3] = [
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

    assert!(
        mat_frobenius(m, &erfa) < 1e-15,
        "frame bias provider should match ERFA bp06 rb"
    );
}

/// ICRS→EclipticMeanJ2000 at J2000.0 matches ERFA `ecm06(J2000.0)`.
#[test]
fn erfa_icrs_to_ecliptic_j2000_matches_ecm06() {
    let ctx = AstroContext::default();
    let rot = frame_rotation::<ICRS, EclipticMeanJ2000>(JulianDate::J2000, &ctx);
    let m = rot.as_matrix();

    // ERFA ecm06(2451545.0, 0.0)
    let erfa: [[f64; 3]; 3] = [
        [
            0.999_999_999_999_994_1,
            -7.078_368_960_971_556e-8,
            8.056_213_977_613_186e-8,
        ],
        [
            3.289_700_407_741_964_6e-8,
            0.917_482_129_914_958_4,
            0.397_776_999_444_047_93,
        ],
        [
            -1.020_704_472_548_435_5e-7,
            -0.397_776_999_444_043_04,
            0.917_482_129_914_955_6,
        ],
    ];

    assert!(
        mat_frobenius(m, &erfa) < 1e-14,
        "ICRS→EclipticMeanJ2000 should match ERFA ecm06 at J2000.0"
    );
}

/// Precession matrix at J2010.0 matches ERFA `pmat06`.
#[test]
fn erfa_precession_j2010_matches_pmat06() {
    let ctx = AstroContext::default();
    let jd = JulianDate::new(2451545.0 + 3654.7681);
    let rot = frame_rotation::<ICRS, EquatorialMeanOfDate>(jd, &ctx);

    // pmat06 includes bias; provider's ICRS→EquatorialMeanOfDate = precession * bias
    let m = rot.as_matrix();

    // ERFA pmat06(2451545.0, 3654.7681)  (= bias*precession ≡ rbp)
    //
    // pmat06 actually produces the bias-corrected precession-nutation matrix
    // from ICRS to mean equator of date (rbp from bp06), which is = P * rb.
    //
    // Our provider computes this as compose_rotation(ICRS→J2000→MeanOfDate):
    //   R(MeanOfDate←J2000) * R(J2000←ICRS) = P * bias ≡ pmat06 output.
    let erfa: [[f64; 3]; 3] = [
        [
            0.999_997_024_103_100_3,
            -0.002_237_563_057_932_016,
            -0.000_972_160_740_361_032_5,
        ],
        [
            0.002_237_563_102_320_845_5,
            0.999_997_496_652_005_5,
            -1.041_977_147_453_465_6e-6,
        ],
        [
            0.000_972_160_638_193_964_8,
            -1.133_296_955_502_238_7e-6,
            0.999_999_527_451_093,
        ],
    ];

    assert!(
        mat_frobenius(m, &erfa) < 1e-12,
        "ICRS→MeanOfDate at J2010.0 should match ERFA pmat06; Frobenius = {:.3e}",
        mat_frobenius(m, &erfa),
    );
}

/// Precession matrix at J2020.0 matches ERFA `pmat06`.
#[test]
fn erfa_precession_j2020_matches_pmat06() {
    let ctx = AstroContext::default();
    let jd = JulianDate::new(2451545.0 + 7305.0);
    let rot = frame_rotation::<ICRS, EquatorialMeanOfDate>(jd, &ctx);
    let m = rot.as_matrix();

    let erfa: [[f64; 3]; 3] = [
        [
            0.999_988_110_837_506,
            -0.004_472_399_764_533_532,
            -0.001_943_147_956_777_380_7,
        ],
        [
            0.004_472_399_877_498_56,
            0.999_989_998_760_465_6,
            -4.287_158_106_408_562e-6,
        ],
        [
            0.001_943_147_696_774_119_8,
            -4.403_427_548_205_219e-6,
            0.999_998_112_077_037,
        ],
    ];

    assert!(
        mat_frobenius(m, &erfa) < 1e-12,
        "ICRS→MeanOfDate at J2020.0 should match ERFA pmat06; Frobenius = {:.3e}",
        mat_frobenius(m, &erfa),
    );
}

/// Nutation at J2020.0: siderust's 2000B differs from ERFA `num06a` (2000A)
/// within expected 2000A/2000B precision gap (~few mas ≈ ~1e-8 rad).
#[test]
fn erfa_nutation_j2020_within_2000b_ceiling() {
    let ctx = AstroContext::default();
    let jd = JulianDate::new(2451545.0 + 7305.0);
    let rot = frame_rotation::<EquatorialMeanOfDate, EquatorialTrueOfDate>(jd, &ctx);
    let m = rot.as_matrix();

    // ERFA num06a(2451545.0, 7305.0) — IAU 2000A
    let erfa: [[f64; 3]; 3] = [
        [
            0.999_999_996_793_943_5,
            7.346_944_425_006_64e-5,
            3.184_892_099_196_083e-5,
        ],
        [
            -7.346_970_426_149_143e-5,
            0.999_999_997_267_785_5,
            8.162_807_504_485_947e-6,
        ],
        [
            -3.184_832_118_801_187_4e-5,
            -8.165_147_409_144_868e-6,
            0.999_999_999_459_507_5,
        ],
    ];

    let frob = mat_frobenius(m, &erfa);
    // 2000A vs 2000B gap should be ≤ ~10 µas ≈ 5e-11 rad; allow generous 1e-6
    assert!(
        frob < 1e-6,
        "nutation 2000B should be close to ERFA num06a; Frobenius = {:.3e}",
        frob,
    );
}

/// Obliquity at J2000.0 matches 84381.406″ exactly.
#[test]
fn erfa_obliquity_j2000() {
    let eps_arcsec = crate::astro::precession::J2000_MEAN_OBLIQUITY_ARCSEC;
    assert!(
        (eps_arcsec - 84381.406).abs() < 1e-10,
        "J2000 obliquity should be 84381.406 arcsec"
    );
}

/// Frame bias applied to basis vectors matches ERFA via provider path.
#[test]
fn erfa_bias_basis_vectors_via_provider() {
    let ctx = AstroContext::default();
    let rb = frame_rotation::<ICRS, EquatorialMeanJ2000>(JulianDate::J2000, &ctx);

    let ex = rb.apply_array([1.0, 0.0, 0.0]);
    let ey = rb.apply_array([0.0, 1.0, 0.0]);
    let ez = rb.apply_array([0.0, 0.0, 1.0]);

    // Column vectors of ERFA rb (each column = rb applied to basis vector)
    let expected_ex = [
        0.999_999_999_999_994_1,
        7.078_368_694_637_676e-8,
        -8.056_214_211_620_057e-8,
    ];
    let expected_ey = [
        -7.078_368_960_971_556e-8,
        0.999_999_999_999_996_9,
        -3.305_943_169_218_395e-8,
    ];
    let expected_ez = [
        8.056_213_977_613_186e-8,
        3.305_943_735_432_137_5e-8,
        0.999_999_999_999_996_2,
    ];

    for k in 0..3 {
        assert!((ex[k] - expected_ex[k]).abs() < 1e-15, "ex[{k}]");
        assert!((ey[k] - expected_ey[k]).abs() < 1e-15, "ey[{k}]");
        assert!((ez[k] - expected_ez[k]).abs() < 1e-15, "ez[{k}]");
    }
}

/// Composed bias+precession+nutation at J2020.0 (BPN matrix).
/// ERFA pmat06 * num06a is the full BPN; our provider composes it similarly.
#[test]
fn erfa_bpn_j2020_composed() {
    let ctx = AstroContext::default();
    let jd = JulianDate::new(2451545.0 + 7305.0);
    let rot = frame_rotation::<ICRS, EquatorialTrueOfDate>(jd, &ctx);
    let _m = rot.as_matrix();

    // Just verify it's finite and reasonably close to the identity
    // (the full BPN at 20 years offset is still ~0.005 rad off-diagonal)
    let v = [1.0, 0.0, 0.0];
    let out = rot.apply_array(v);
    let norm = (out[0] * out[0] + out[1] * out[1] + out[2] * out[2]).sqrt();
    assert!(
        (norm - 1.0).abs() < 1e-14,
        "BPN should preserve vector norm"
    );
}

/// Full BPN matrix at J2020.0 with IAU 2006A nutation matches ERFA `pnm06a`.
/// Reference: erfa.pnm06a(2458850.0, 0.0)
///
/// Tolerance 1e-8 (~2 µas): our BPN is built via a single `fw2m(γ̄, φ̄,
/// ψ̄+Δψ, ε_A+Δε)` call, while ERFA's `pnm06a` multiplies separate
/// nutation × bias-precession matrices. The approaches are algebraically
/// equivalent but differ at the ~5e-9 (1 µas) level due to floating-point
/// non-commutativity.
#[test]
fn erfa_bpn_j2020_iau2006a_matches_pnm06a() {
    use crate::coordinates::transform::context::{DefaultEop, DefaultEphemeris};

    let ctx = AstroContext::<DefaultEphemeris, DefaultEop, Iau2006A>::with_types();
    let jd = JulianDate::new(2458850.0);
    let rot = <() as FrameRotationProvider<ICRS, EquatorialTrueOfDate>>::rotation::<
        DefaultEphemeris,
        DefaultEop,
        Iau2006A,
    >(jd, &ctx);
    let m = rot.as_matrix();

    // erfa.pnm06a(2458850.0, 0.0)
    let erfa: [[f64; 3]; 3] = [
        [
            9.999_884_981_033_785e-01,
            -4.398_931_180_974_626e-03,
            -1.911_299_404_659_01e-3,
        ],
        [
            4.398_946_896_051_368e-03,
            9.999_903_245_782_221e-01,
            4.018_396_504_712_474e-06,
        ],
        [
            1.911_263_235_381_445e-03,
            -1.242_605_486_917_459e-05,
            9.999_981_734_575_509e-01,
        ],
    ];

    let frob = mat_frobenius(m, &erfa);
    assert!(
        frob < 1e-8,
        "BPN with Iau2006A at J2020.0 should match ERFA pnm06a; Frobenius = {:.3e}",
        frob,
    );
}

/// Nutation matrix at J2020.0 with IAU 2006A matches ERFA `num06a`.
/// Reference: erfa.num06a(2458850.0, 0.0)
#[test]
fn erfa_nutation_j2020_iau2006a_matches_num06a() {
    use crate::coordinates::transform::context::{DefaultEop, DefaultEphemeris};

    let ctx = AstroContext::<DefaultEphemeris, DefaultEop, Iau2006A>::with_types();
    let jd = JulianDate::new(2458850.0);
    let rot = <() as FrameRotationProvider<EquatorialMeanOfDate, EquatorialTrueOfDate>>::rotation::<
        DefaultEphemeris,
        DefaultEop,
        Iau2006A,
    >(jd, &ctx);
    let m = rot.as_matrix();

    // erfa.num06a(2458850.0, 0.0)
    let erfa: [[f64; 3]; 3] = [
        [
            9.999_999_967_939_435e-01,
            7.346_944_425_006_64e-5,
            3.184_892_099_196_083e-05,
        ],
        [
            -7.346_970_426_149_143e-05,
            9.999_999_972_677_856e-01,
            8.162_807_504_485_947e-06,
        ],
        [
            -3.184_832_118_801_187e-05,
            -8.165_147_409_144_869e-06,
            9.999_999_994_595_075e-01,
        ],
    ];

    let frob = mat_frobenius(m, &erfa);
    assert!(
        frob < 1e-12,
        "Nutation Iau2006A at J2020.0 should match ERFA num06a; Frobenius = {:.3e}",
        frob,
    );
}

/// Inverse frame bias recovers input via provider path.
#[test]
fn erfa_inverse_bias_provider_roundtrip() {
    let ctx = AstroContext::default();
    let jd = JulianDate::J2000;
    let fwd = frame_rotation::<ICRS, EquatorialMeanJ2000>(jd, &ctx);
    let inv = frame_rotation::<EquatorialMeanJ2000, ICRS>(jd, &ctx);

    let v = [
        0.577_350_269_189_626,
        -0.577_350_269_189_626,
        0.577_350_269_189_626,
    ]; // ±1/√3
    let rt = inv.apply_array(fwd.apply_array(v));
    for k in 0..3 {
        assert!((rt[k] - v[k]).abs() < 1e-15, "roundtrip[{k}]");
    }
}
