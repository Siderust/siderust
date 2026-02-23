// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Smoke tests for the siderust-ffi C API.
//!
//! These tests exercise the lifecycle (create → query → free) of every major
//! handle type and verify that the null-pointer guards and invalid-period
//! checks return the expected status codes.
//!
//! Run with:
//!   cargo test --manifest-path siderust/siderust-ffi/Cargo.toml
//!
//! And optionally with the address sanitizer to catch leaks:
//!   RUSTFLAGS="-Z sanitizer=address" cargo +nightly test --target x86_64-unknown-linux-gnu

// The query functions are `extern "C"` (not `unsafe extern "C"`) so calling
// them from Rust is safe.  We keep the `unsafe` blocks for documentation
// clarity (raw pointers are involved), hence we suppress this lint.
#![allow(unused_unsafe)]

use siderust_ffi::*;
use std::ffi::CString;
use std::ptr;

// ─────────────────────────────────────────────────────────────────────────────
// Opaque handle lifecycle
// ─────────────────────────────────────────────────────────────────────────────

#[test]
fn star_catalog_create_and_free() {
    let name = CString::new("VEGA").unwrap();
    let mut handle: *mut SiderustStar = ptr::null_mut();
    let status = unsafe { siderust_star_catalog(name.as_ptr(), &mut handle) };
    assert_eq!(status, SiderustStatus::Ok);
    assert!(!handle.is_null());
    unsafe { siderust_star_free(handle) };
}

#[test]
fn star_unknown_name_returns_error() {
    let name = CString::new("NOTASTAR").unwrap();
    let mut handle: *mut SiderustStar = ptr::null_mut();
    let status = unsafe { siderust_star_catalog(name.as_ptr(), &mut handle) };
    assert_eq!(status, SiderustStatus::UnknownStar);
    assert!(handle.is_null());
}

#[test]
fn star_custom_create_and_free() {
    let name = CString::new("ProximaCen").unwrap();
    let mut handle: *mut SiderustStar = ptr::null_mut();
    let status = unsafe {
        siderust_star_create(
            name.as_ptr(),
            4.244,     // distance_ly
            0.12,      // mass_solar
            0.15,      // radius_solar
            0.0017,    // luminosity_solar
            217.4291,  // ra_deg
            -62.679,   // dec_deg
            2451545.0, // epoch_jd (J2000.0)
            ptr::null(),
            &mut handle,
        )
    };
    assert_eq!(status, SiderustStatus::Ok);
    assert!(!handle.is_null());
    unsafe { siderust_star_free(handle) };
}

#[test]
fn target_create_and_free() {
    // Orion Nebula direction
    let mut handle: *mut SiderustTarget = ptr::null_mut();
    let status = unsafe { siderust_target_create(83.82, -5.39, 2451545.0, &mut handle) };
    assert_eq!(status, SiderustStatus::Ok);
    assert!(!handle.is_null());
    let mut ra = 0.0_f64;
    let mut dec = 0.0_f64;
    unsafe {
        assert_eq!(siderust_target_ra_deg(handle, &mut ra), SiderustStatus::Ok);
        assert_eq!(
            siderust_target_dec_deg(handle, &mut dec),
            SiderustStatus::Ok
        );
    }
    assert!((ra - 83.82).abs() < 1e-10);
    assert!((dec - (-5.39)).abs() < 1e-10);
    unsafe { siderust_target_free(handle) };
}

// ─────────────────────────────────────────────────────────────────────────────
// Vec-returning functions: allocate then free
// ─────────────────────────────────────────────────────────────────────────────

fn paris_observer() -> SiderustGeodetict {
    SiderustGeodetict {
        lon_deg: 2.35,
        lat_deg: 48.85,
        height_m: 35.0,
    }
}

fn one_day_window() -> TempochPeriodMjd {
    TempochPeriodMjd {
        start_mjd: 60000.0,
        end_mjd: 60001.0,
    }
}

fn default_opts() -> SiderustSearchOpts {
    SiderustSearchOpts {
        time_tolerance_days: 1e-9,
        scan_step_days: 0.0,
        has_scan_step: false,
    }
}

#[test]
fn sun_crossings_alloc_and_free() {
    let mut out: *mut SiderustCrossingEvent = ptr::null_mut();
    let mut count: usize = 0;
    let status = unsafe {
        siderust_sun_crossings(
            paris_observer(),
            one_day_window(),
            0.0,
            default_opts(),
            &mut out,
            &mut count,
        )
    };
    assert_eq!(status, SiderustStatus::Ok);
    unsafe { siderust_crossings_free(out, count) };
}

#[test]
fn sun_culminations_alloc_and_free() {
    let mut out: *mut SiderustCulminationEvent = ptr::null_mut();
    let mut count: usize = 0;
    let status = unsafe {
        siderust_sun_culminations(
            paris_observer(),
            one_day_window(),
            default_opts(),
            &mut out,
            &mut count,
        )
    };
    assert_eq!(status, SiderustStatus::Ok);
    unsafe { siderust_culminations_free(out, count) };
}

#[test]
fn sun_periods_above_threshold_alloc_and_free() {
    let mut out: *mut TempochPeriodMjd = ptr::null_mut();
    let mut count: usize = 0;
    let status = unsafe {
        siderust_sun_above_threshold(
            paris_observer(),
            one_day_window(),
            0.0,
            default_opts(),
            &mut out,
            &mut count,
        )
    };
    assert_eq!(status, SiderustStatus::Ok);
    unsafe { siderust_periods_free(out, count) };
}

#[test]
fn moon_azimuth_crossings_alloc_and_free() {
    let mut out: *mut SiderustAzimuthCrossingEvent = ptr::null_mut();
    let mut count: usize = 0;
    let status = unsafe {
        siderust_moon_azimuth_crossings(
            paris_observer(),
            one_day_window(),
            180.0, // bearing south
            default_opts(),
            &mut out,
            &mut count,
        )
    };
    assert_eq!(status, SiderustStatus::Ok);
    unsafe { siderust_azimuth_crossings_free(out, count) };
}

#[test]
fn phase_events_alloc_and_free() {
    let window = TempochPeriodMjd {
        start_mjd: 60000.0,
        end_mjd: 60030.0, // ~1 month to be sure we get events
    };
    let mut out: *mut SiderustPhaseEvent = ptr::null_mut();
    let mut count: usize = 0;
    let status =
        unsafe { siderust_find_phase_events(window, default_opts(), &mut out, &mut count) };
    assert_eq!(status, SiderustStatus::Ok);
    // One month should contain at least 2 principal phases
    assert!(
        count >= 2,
        "Expected >=2 phase events in a 30-day window, got {}",
        count
    );
    unsafe { siderust_phase_events_free(out, count) };
}

// ─────────────────────────────────────────────────────────────────────────────
// Null-pointer guards
// ─────────────────────────────────────────────────────────────────────────────

#[test]
fn null_ptr_guard_altitude_out() {
    let status = unsafe { siderust_sun_altitude_at(paris_observer(), 60000.0, ptr::null_mut()) };
    assert_eq!(status, SiderustStatus::NullPointer);
}

#[test]
fn null_ptr_guard_crossings_out() {
    let status = unsafe {
        siderust_sun_crossings(
            paris_observer(),
            one_day_window(),
            0.0,
            default_opts(),
            ptr::null_mut(),
            ptr::null_mut(),
        )
    };
    assert_eq!(status, SiderustStatus::NullPointer);
}

#[test]
fn null_ptr_guard_star_handle() {
    let mut out: f64 = 0.0;
    let status =
        unsafe { siderust_star_altitude_at(ptr::null(), paris_observer(), 60000.0, &mut out) };
    assert_eq!(status, SiderustStatus::NullPointer);
}

#[test]
fn null_ptr_guard_target_handle() {
    let mut out: f64 = 0.0;
    let status =
        unsafe { siderust_target_altitude_at(ptr::null(), paris_observer(), 60000.0, &mut out) };
    assert_eq!(status, SiderustStatus::NullPointer);
}

// ─────────────────────────────────────────────────────────────────────────────
// Validation errors
// ─────────────────────────────────────────────────────────────────────────────

#[test]
fn invalid_period_end_before_start() {
    let bad_window = TempochPeriodMjd {
        start_mjd: 60001.0,
        end_mjd: 60000.0, // end < start
    };
    let mut out: *mut TempochPeriodMjd = ptr::null_mut();
    let mut count: usize = 0;
    let status = unsafe {
        siderust_sun_above_threshold(
            paris_observer(),
            bad_window,
            0.0,
            default_opts(),
            &mut out,
            &mut count,
        )
    };
    assert_eq!(status, SiderustStatus::InvalidPeriod);
}

#[test]
fn icrs_azimuth_wrong_frame_returns_error() {
    let dir = SiderustSphericalDir {
        lon_deg: 83.82,
        lat_deg: -5.39,
        frame: SiderustFrame::EquatorialMeanJ2000, // not ICRS
    };
    let mut out: f64 = 0.0;
    let status = unsafe { siderust_icrs_azimuth_at(dir, paris_observer(), 60000.0, &mut out) };
    assert_eq!(status, SiderustStatus::InvalidFrame);
}
