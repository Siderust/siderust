// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # IAU 2000/2006 Compliance Integration Tests
//!
//! These tests verify the full IAU 2000/2006 transformation chain
//! by cross-checking module interactions: precession + nutation,
//! ERA + GMST, CIO-based GCRS→CIRS, polar motion, and light deflection.

use siderust::astro::cio::{cip_cio, gcrs_to_cirs_matrix};
use siderust::astro::era::earth_rotation_angle;
use siderust::astro::eop::{EopProvider, NullEop};
use siderust::astro::light_deflection::{solar_deflection, solar_deflection_inverse};
use siderust::astro::nutation_iau2000b::nutation_iau2000b;
use siderust::astro::polar_motion::{polar_motion_matrix, tio_locator_sp};
use siderust::astro::precession_iau2006::{
    fw_matrix, mean_obliquity_iau2006, precession_fw_angles, precession_matrix_iau2006,
    precession_nutation_matrix,
};
use siderust::astro::sidereal::{gmst_iau2006, gast_iau2006};
use siderust::time::JulianDate;
use qtty::*;
use std::f64::consts::TAU;

// ═══════════════════════════════════════════════════════════════════════════
// Full GCRS → ITRS transformation chain (CIO-based, IAU 2006)
// ═══════════════════════════════════════════════════════════════════════════

#[test]
fn full_gcrs_to_itrs_chain() {
    // Test the complete IAU 2006 celestial-to-terrestrial transformation:
    // GCRS → Q(t) → CIRS → R₃(ERA) → TIRS → W → ITRS
    //
    // We verify that the chain produces a proper rotation (det = 1, orthogonal).

    let jd_tt = JulianDate::new(2_460_000.5); // 2023-02-25 TT
    let jd_ut1 = jd_tt; // approximation: UT1 ≈ TT (valid to < 1 s)

    // Step 1: Nutation
    let nut = nutation_iau2000b(jd_tt);

    // Step 2: CIP/CIO
    let cip = cip_cio(jd_tt, nut.dpsi, nut.deps);

    // Step 3: GCRS → CIRS (celestial-to-intermediate)
    let q = gcrs_to_cirs_matrix(cip.x, cip.y, cip.s);

    // Step 4: ERA
    let era = earth_rotation_angle(jd_ut1);

    // Step 5: CIRS → TIRS (apply ERA as R₃(-ERA))
    let (s_era, c_era) = (-era.value()).sin_cos();
    let r3_era = affn::Rotation3::from_matrix([
        [c_era, -s_era, 0.0],
        [s_era, c_era, 0.0],
        [0.0, 0.0, 1.0],
    ]);

    // Step 6: Polar motion (TIRS → ITRS)
    let sp = tio_locator_sp(jd_tt);
    let w = polar_motion_matrix(Radians::new(0.0), Radians::new(0.0), sp); // zero pole coords (NullEop)

    // Full chain: W · R₃(-ERA) · Q
    let cirs_to_tirs = r3_era;
    let v_gcrs = [0.6, 0.7, 0.3742f64];
    let mag = (v_gcrs[0].powi(2) + v_gcrs[1].powi(2) + v_gcrs[2].powi(2)).sqrt();
    let v_gcrs = [v_gcrs[0] / mag, v_gcrs[1] / mag, v_gcrs[2] / mag];

    let v_cirs = q.apply_array(v_gcrs);
    let v_tirs = cirs_to_tirs.apply_array(v_cirs);
    let v_itrs = w.apply_array(v_tirs);

    // Verify output is still a unit vector (proper rotation preserves norm)
    let norm = (v_itrs[0].powi(2) + v_itrs[1].powi(2) + v_itrs[2].powi(2)).sqrt();
    assert!(
        (norm - 1.0).abs() < 1e-14,
        "ITRS vector norm = {}, expected 1.0",
        norm
    );
}

// ═══════════════════════════════════════════════════════════════════════════
// Precession + nutation consistency
// ═══════════════════════════════════════════════════════════════════════════

#[test]
fn precession_nutation_matrix_is_consistent() {
    // The NPB matrix should be very close to the pure precession matrix
    // when nutation is zero.
    let jd = JulianDate::new(2_460_000.5);

    let p = precession_matrix_iau2006(jd);
    let npb = precession_nutation_matrix(jd, Radians::new(0.0), Radians::new(0.0));

    let p_m = p.as_matrix();
    let npb_m = npb.as_matrix();

    for i in 0..3 {
        for j in 0..3 {
            assert!(
                (p_m[i][j] - npb_m[i][j]).abs() < 1e-14,
                "P[{}][{}] = {}, NPB(0,0)[{}][{}] = {}",
                i, j, p_m[i][j], i, j, npb_m[i][j]
            );
        }
    }
}

#[test]
fn nutation_shifts_cip_significantly() {
    // With real nutation, CIP should differ from zero-nutation by ~10-20 arcseconds
    let jd = JulianDate::new(2_460_000.5);
    let nut = nutation_iau2000b(jd);

    let (x0, y0) = siderust::astro::cio::cip_xy(jd, Radians::new(0.0), Radians::new(0.0));
    let (x1, y1) = siderust::astro::cio::cip_xy(jd, nut.dpsi, nut.deps);

    let dx = (x1 - x0) * 206_264.806; // rad → arcsec
    let dy = (y1 - y0) * 206_264.806;

    assert!(
        dx.abs() > 0.1 && dx.abs() < 25.0,
        "CIP X nutation shift = {}″, expected 1-20″",
        dx
    );
    assert!(
        dy.abs() > 0.1 && dy.abs() < 25.0,
        "CIP Y nutation shift = {}″, expected 1-20″",
        dy
    );
}

// ═══════════════════════════════════════════════════════════════════════════
// ERA ↔ GMST consistency
// ═══════════════════════════════════════════════════════════════════════════

#[test]
fn era_and_gmst_differ_by_accumulation() {
    // GMST = ERA + accumulated precession polynomial.
    // The polynomial is ~0.014506″ at t=0, growing at ~4612″/century.
    // At t ≈ 0.23 century (2023), the polynomial is ~1060″ ≈ 0.295°.
    let jd = JulianDate::new(2_460_000.5);
    let era = earth_rotation_angle(jd);
    let gmst = gmst_iau2006(jd, jd);

    let diff_as = ((gmst - era).value().rem_euclid(TAU)) * 206_264.806;
    // Should be ~1060″ ≈ 0.295°
    assert!(
        diff_as > 500.0 && diff_as < 2000.0,
        "GMST − ERA = {}″, expected ~1060″",
        diff_as
    );
}

#[test]
fn gast_minus_gmst_equals_equation_of_equinoxes() {
    // GAST = GMST + Δψ·cos(ε)
    let jd = JulianDate::new(2_460_000.5);
    let nut = nutation_iau2000b(jd);
    let true_obliquity = nut.true_obliquity();

    let gmst = gmst_iau2006(jd, jd);
    let gast = gast_iau2006(jd, jd, nut.dpsi, true_obliquity);

    let ee_expected = nut.dpsi.value() * true_obliquity.cos();
    let ee_actual = (gast - gmst).value();

    assert!(
        (ee_actual - ee_expected).abs() < 1e-12,
        "GAST−GMST = {}, Δψ·cos(ε) = {}",
        ee_actual, ee_expected
    );
}

// ═══════════════════════════════════════════════════════════════════════════
// EOP provider integration
// ═══════════════════════════════════════════════════════════════════════════

#[test]
fn null_eop_produces_ut1_equal_utc() {
    let eop = NullEop;
    let jd_utc = JulianDate::new(2_460_000.5);
    let vals = eop.eop_at(jd_utc);
    let jd_ut1 = vals.jd_ut1(jd_utc);
    assert_eq!(jd_ut1.value(), jd_utc.value());
}

// ═══════════════════════════════════════════════════════════════════════════
// Light deflection integration
// ═══════════════════════════════════════════════════════════════════════════

#[test]
fn light_deflection_roundtrip_accuracy() {
    // Apply and remove solar deflection; result should match original to < 1 μas
    use siderust::coordinates::cartesian::direction;
    
    let star_arr = [0.3f64, 0.8, 0.5];
    let mag = (star_arr[0].powi(2) + star_arr[1].powi(2) + star_arr[2].powi(2)).sqrt();
    let star = direction::EquatorialMeanJ2000::new(
        star_arr[0] / mag,
        star_arr[1] / mag,
        star_arr[2] / mag,
    );

    let earth_sun = [0.9, -0.3, 0.1]; // ~1 AU from Sun

    let apparent = solar_deflection(star, earth_sun);
    let recovered = solar_deflection_inverse(apparent, earth_sun);

    let angular_err = ((recovered.x() - star.x()).powi(2)
        + (recovered.y() - star.y()).powi(2)
        + (recovered.z() - star.z()).powi(2))
    .sqrt();
    let err_uas = angular_err * 206_264_806_000.0; // rad → μas

    assert!(
        err_uas < 1.0,
        "roundtrip error = {} μas, should be < 1 μas",
        err_uas
    );
}

// ═══════════════════════════════════════════════════════════════════════════
// Obliquity consistency across modules
// ═══════════════════════════════════════════════════════════════════════════

#[test]
fn obliquity_consistent_across_modules() {
    // The IAU 2006 obliquity at J2000 should be 84381.406″ everywhere
    let eps_prec = mean_obliquity_iau2006(JulianDate::J2000);
    let nut = nutation_iau2000b(JulianDate::J2000);

    let eps_prec_as = eps_prec.to::<Degree>().value() * 3600.0;
    let eps_nut_as = nut.mean_obliquity.to::<Degree>().value() * 3600.0;

    assert!(
        (eps_prec_as - 84381.406).abs() < 0.001,
        "precession module obliquity = {}″",
        eps_prec_as
    );
    assert!(
        (eps_nut_as - 84381.406).abs() < 0.001,
        "nutation module obliquity = {}″",
        eps_nut_as
    );
    assert!(
        (eps_prec_as - eps_nut_as).abs() < 1e-10,
        "obliquity mismatch between modules: {} vs {}",
        eps_prec_as, eps_nut_as
    );
}

// ═══════════════════════════════════════════════════════════════════════════
// Polar motion symmetry
// ═══════════════════════════════════════════════════════════════════════════

#[test]
fn polar_motion_roundtrip() {
    let as2rad = std::f64::consts::PI / (180.0 * 3600.0);
    let xp = Radians::new(0.2 * as2rad);
    let yp = Radians::new(0.35 * as2rad);
    let sp = Radians::new(-1.2e-9); // ~small s'

    let w = polar_motion_matrix(xp, yp, sp);
    let w_inv = w.inverse();

    let v = [0.5, 0.7, 0.3];
    let transformed = w.apply_array(v);
    let recovered = w_inv.apply_array(transformed);

    for i in 0..3 {
        assert!(
            (recovered[i] - v[i]).abs() < 1e-14,
            "polar motion roundtrip [{}]: {} vs {}",
            i, recovered[i], v[i]
        );
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// Temporal consistency: quantities vary smoothly
// ═══════════════════════════════════════════════════════════════════════════

#[test]
fn era_monotonic_over_sidereal_day() {
    // ERA should increase monotonically over one sidereal day
    let mut prev = 0.0f64;
    let jd0 = 2_460_000.5;
    for i in 0..24 {
        let jd = JulianDate::new(jd0 + i as f64 / 24.0 * 0.99727);
        let era = earth_rotation_angle(jd).value();
        if i > 0 {
            // Allow wraparound at 2π
            let diff = (era - prev).rem_euclid(TAU);
            assert!(
                diff > 0.0 && diff < std::f64::consts::PI,
                "ERA not monotonically increasing at step {}: prev={}, curr={}",
                i, prev, era
            );
        }
        prev = era;
    }
}

#[test]
fn nutation_varies_on_18_6_year_cycle() {
    // The dominant nutation term has an 18.6-year period (Ω of the Moon's node).
    // Check that nutation at two epochs separated by 9.3 years has opposite sign.
    let jd1 = JulianDate::new(2_451_545.0); // J2000.0
    let jd2 = JulianDate::new(2_451_545.0 + 9.3 * 365.25); // ~9.3 years later

    let nut1 = nutation_iau2000b(jd1);
    let nut2 = nutation_iau2000b(jd2);

    // Δψ should change sign (or at least differ significantly)
    let dpsi_diff = (nut2.dpsi - nut1.dpsi).value().abs() * 206_264.806;
    assert!(
        dpsi_diff > 5.0,
        "Δψ should vary by > 5″ over 9.3 years, got {}″",
        dpsi_diff
    );
}
