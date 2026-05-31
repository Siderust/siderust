// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Regression tests for the dynamics FFI surface.
//!
//! These tests exercise the full FFI lifecycle using `extern "C"` function
//! calls, validating that:
//!
//! - Handle creation and cleanup round-trips correctly.
//! - `siderust_propagator_propagate` closes a circular LEO orbit to within
//!   1e-3 km over one full orbital period.

use siderust_ffi::{
    siderust_dynamics_context_free, siderust_dynamics_context_new, siderust_orbit_state_epoch_jd,
    siderust_orbit_state_free, siderust_orbit_state_new, siderust_orbit_state_position,
    siderust_orbit_state_velocity, siderust_propagator_free, siderust_propagator_propagate,
    siderust_propagator_two_body_earth_new, siderust_propagator_two_body_new,
    SiderustDynamicsStatus, SiderustStatus,
};

// Earth GM (EGM2008) in km³/s²
const MU: f64 = 398_600.441_8;
// Orbit radius km
const R: f64 = 7_000.0;

fn orbital_period(r: f64, mu: f64) -> f64 {
    2.0 * std::f64::consts::PI * (r.powi(3) / mu).sqrt()
}

// ── Lifecycle / null-pointer smoke tests ──────────────────────────────────

#[test]
fn dynamics_context_create_free() {
    let mut ctx: *mut siderust_ffi::SiderustDynamicsContext = std::ptr::null_mut();
    let status = siderust_dynamics_context_new(&mut ctx);
    assert_eq!(status, SiderustStatus::Ok, "context creation failed");
    assert!(!ctx.is_null());
    unsafe { siderust_dynamics_context_free(ctx) };
}

#[test]
fn dynamics_context_new_null_out_returns_null_pointer() {
    let status = siderust_dynamics_context_new(std::ptr::null_mut());
    assert_eq!(status, SiderustStatus::NullPointer);
}

#[test]
fn orbit_state_create_free() {
    let jd = 2_451_545.0_f64;
    let mut s: *mut siderust_ffi::SiderustOrbitState = std::ptr::null_mut();
    let status = siderust_orbit_state_new(jd, R, 0.0, 0.0, 0.0, 7.5, 0.0, &mut s);
    assert_eq!(status, SiderustStatus::Ok);
    assert!(!s.is_null());
    unsafe { siderust_orbit_state_free(s) };
}

#[test]
fn orbit_state_accessors_roundtrip() {
    let jd_in = 2_451_545.0_f64;
    let (x0, y0, z0) = (7_000.0, 100.0, -50.0);
    let (vx0, vy0, vz0) = (0.1, 7.5, 0.3);

    let mut s: *mut siderust_ffi::SiderustOrbitState = std::ptr::null_mut();
    siderust_orbit_state_new(jd_in, x0, y0, z0, vx0, vy0, vz0, &mut s);
    assert!(!s.is_null());

    let mut x = 0.0_f64;
    let mut y = 0.0_f64;
    let mut z = 0.0_f64;
    let ps = siderust_orbit_state_position(s, &mut x, &mut y, &mut z);
    assert_eq!(ps, SiderustStatus::Ok);
    assert!((x - x0).abs() < 1e-12);
    assert!((y - y0).abs() < 1e-12);
    assert!((z - z0).abs() < 1e-12);

    let mut vx = 0.0_f64;
    let mut vy = 0.0_f64;
    let mut vz = 0.0_f64;
    let vs = siderust_orbit_state_velocity(s, &mut vx, &mut vy, &mut vz);
    assert_eq!(vs, SiderustStatus::Ok);
    assert!((vx - vx0).abs() < 1e-12);
    assert!((vy - vy0).abs() < 1e-12);
    assert!((vz - vz0).abs() < 1e-12);

    let mut jd_out = 0.0_f64;
    let es = siderust_orbit_state_epoch_jd(s, &mut jd_out);
    assert_eq!(es, SiderustStatus::Ok);
    assert!((jd_out - jd_in).abs() < 1e-9, "JD roundtrip: {jd_out}");

    unsafe { siderust_orbit_state_free(s) };
}

#[test]
fn propagator_earth_create_free() {
    let mut p: *mut siderust_ffi::SiderustPropagator = std::ptr::null_mut();
    let status = siderust_propagator_two_body_earth_new(&mut p);
    assert_eq!(status, SiderustStatus::Ok);
    assert!(!p.is_null());
    unsafe { siderust_propagator_free(p) };
}

#[test]
fn propagator_custom_gm_create_free() {
    let mut p: *mut siderust_ffi::SiderustPropagator = std::ptr::null_mut();
    let status = siderust_propagator_two_body_new(MU, &mut p);
    assert_eq!(status, SiderustStatus::Ok);
    assert!(!p.is_null());
    unsafe { siderust_propagator_free(p) };
}

// ── Main regression: circular orbit closes within 1e-3 km ────────────────

/// Build a circular LEO state, propagate one full period with TwoBody DOP853,
/// and verify the position closure error is < 1e-3 km.
#[test]
fn two_body_circular_orbit_closes_within_1e3_km() {
    // Circular velocity
    let v_circ = (MU / R).sqrt();
    let period = orbital_period(R, MU);
    let jd = 2_451_545.0_f64;

    // 1. Create initial orbit state
    let mut s0: *mut siderust_ffi::SiderustOrbitState = std::ptr::null_mut();
    let st = siderust_orbit_state_new(jd, R, 0.0, 0.0, 0.0, v_circ, 0.0, &mut s0);
    assert_eq!(st, SiderustStatus::Ok, "orbit state creation failed");

    // 2. Create propagator (Earth GM, DOP853)
    let mut prop: *mut siderust_ffi::SiderustPropagator = std::ptr::null_mut();
    let pt = siderust_propagator_two_body_earth_new(&mut prop);
    assert_eq!(pt, SiderustStatus::Ok, "propagator creation failed");

    // 3. Propagate one orbital period (null ctx → empty context)
    let mut s_final: *mut siderust_ffi::SiderustOrbitState = std::ptr::null_mut();
    let ds = siderust_propagator_propagate(prop, s0, period, std::ptr::null(), &mut s_final);
    assert_eq!(
        ds,
        SiderustDynamicsStatus::Ok,
        "propagation failed with status {ds:?}"
    );
    assert!(!s_final.is_null());

    // 4. Read final position and check closure
    let mut xf = 0.0_f64;
    let mut yf = 0.0_f64;
    let mut zf = 0.0_f64;
    siderust_orbit_state_position(s_final, &mut xf, &mut yf, &mut zf);

    let dx = xf - R;
    let dy = yf - 0.0;
    let dz = zf - 0.0;
    let closure_km = (dx * dx + dy * dy + dz * dz).sqrt();

    assert!(
        closure_km < 1e-3,
        "orbit closure error {closure_km:.2e} km exceeds 1e-3 km threshold"
    );

    // 5. Cleanup
    unsafe {
        siderust_orbit_state_free(s_final);
        siderust_orbit_state_free(s0);
        siderust_propagator_free(prop);
    }
}

/// Same test but using an explicit dynamics context handle (not null).
#[test]
fn two_body_with_explicit_context_closes() {
    let v_circ = (MU / R).sqrt();
    let period = orbital_period(R, MU);

    let mut s0: *mut siderust_ffi::SiderustOrbitState = std::ptr::null_mut();
    siderust_orbit_state_new(2_451_545.0, R, 0.0, 0.0, 0.0, v_circ, 0.0, &mut s0);

    let mut prop: *mut siderust_ffi::SiderustPropagator = std::ptr::null_mut();
    siderust_propagator_two_body_earth_new(&mut prop);

    let mut ctx: *mut siderust_ffi::SiderustDynamicsContext = std::ptr::null_mut();
    siderust_dynamics_context_new(&mut ctx);

    let mut s_final: *mut siderust_ffi::SiderustOrbitState = std::ptr::null_mut();
    let ds = siderust_propagator_propagate(prop, s0, period, ctx, &mut s_final);
    assert_eq!(ds, SiderustDynamicsStatus::Ok, "propagation failed: {ds:?}");

    let mut xf = 0.0_f64;
    let mut yf = 0.0_f64;
    let mut zf = 0.0_f64;
    siderust_orbit_state_position(s_final, &mut xf, &mut yf, &mut zf);

    let dx = xf - R;
    let dy = yf;
    let dz = zf;
    let closure_km = (dx * dx + dy * dy + dz * dz).sqrt();
    assert!(
        closure_km < 1e-3,
        "orbit closure error {closure_km:.2e} km exceeds 1e-3 km"
    );

    unsafe {
        siderust_orbit_state_free(s_final);
        siderust_orbit_state_free(s0);
        siderust_propagator_free(prop);
        siderust_dynamics_context_free(ctx);
    }
}

/// Backward propagation: propagate −1 period and check closure.
#[test]
fn two_body_backward_propagation_closes() {
    let v_circ = (MU / R).sqrt();
    let period = orbital_period(R, MU);

    let mut s0: *mut siderust_ffi::SiderustOrbitState = std::ptr::null_mut();
    siderust_orbit_state_new(2_451_545.0, R, 0.0, 0.0, 0.0, v_circ, 0.0, &mut s0);

    let mut prop: *mut siderust_ffi::SiderustPropagator = std::ptr::null_mut();
    siderust_propagator_two_body_earth_new(&mut prop);

    let mut s_final: *mut siderust_ffi::SiderustOrbitState = std::ptr::null_mut();
    let ds = siderust_propagator_propagate(prop, s0, -period, std::ptr::null(), &mut s_final);
    assert_eq!(
        ds,
        SiderustDynamicsStatus::Ok,
        "backward propagation failed: {ds:?}"
    );

    let mut xf = 0.0_f64;
    let mut yf = 0.0_f64;
    let mut zf = 0.0_f64;
    siderust_orbit_state_position(s_final, &mut xf, &mut yf, &mut zf);

    let dx = xf - R;
    let dy = yf;
    let dz = zf;
    let closure_km = (dx * dx + dy * dy + dz * dz).sqrt();
    assert!(
        closure_km < 1e-3,
        "backward orbit closure error {closure_km:.2e} km exceeds 1e-3 km"
    );

    unsafe {
        siderust_orbit_state_free(s_final);
        siderust_orbit_state_free(s0);
        siderust_propagator_free(prop);
    }
}
