// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Smoke tests for the reduced siderust-ffi C API.
//!
//! These tests focus on the forward-looking surface: `siderust_subject_t`,
//! `SiderustGenericTarget`, array lifetimes, and the remaining non-subject
//! helpers such as phase/event APIs.

#![allow(unused_unsafe)]

use siderust_ffi::*;
use std::ffi::CString;
use std::ptr;

fn mjd(value: f64) -> TempochMjd {
    TempochMjd::new(value)
}

fn paris_observer() -> SiderustGeodetict {
    SiderustGeodetict {
        lon_deg: 2.35,
        lat_deg: 48.85,
        height_m: 35.0,
    }
}

fn one_day_window() -> TempochPeriodMjd {
    TempochPeriodMjd {
        start_mjd: mjd(60000.0),
        end_mjd: mjd(60001.0),
    }
}

fn default_opts() -> SiderustSearchOpts {
    SiderustSearchOpts {
        time_tolerance_days: 1e-9,
        scan_step_days: 0.0,
        has_scan_step: false,
    }
}

fn sun_subject() -> SiderustSubject {
    SiderustSubject {
        kind: SiderustSubjectKind::Body,
        body: SiderustBody::Sun,
        star_handle: ptr::null(),
        icrs_dir: SiderustSphericalDir {
            polar_deg: 0.0,
            azimuth_deg: 0.0,
            frame: SiderustFrame::ICRS,
        },
        generic_target_handle: ptr::null(),
    }
}

fn moon_subject() -> SiderustSubject {
    SiderustSubject {
        kind: SiderustSubjectKind::Body,
        body: SiderustBody::Moon,
        star_handle: ptr::null(),
        icrs_dir: SiderustSphericalDir {
            polar_deg: 0.0,
            azimuth_deg: 0.0,
            frame: SiderustFrame::ICRS,
        },
        generic_target_handle: ptr::null(),
    }
}

fn icrs_subject() -> SiderustSubject {
    SiderustSubject {
        kind: SiderustSubjectKind::Icrs,
        body: SiderustBody::Sun,
        star_handle: ptr::null(),
        icrs_dir: SiderustSphericalDir {
            polar_deg: 38.78,
            azimuth_deg: 279.23,
            frame: SiderustFrame::ICRS,
        },
        generic_target_handle: ptr::null(),
    }
}

fn vega_subject() -> (*mut SiderustStar, SiderustSubject) {
    let name = CString::new("VEGA").unwrap();
    let mut handle: *mut SiderustStar = ptr::null_mut();
    let status = unsafe { siderust_star_catalog(name.as_ptr(), &mut handle) };
    assert_eq!(status, SiderustStatus::Ok);
    (
        handle,
        SiderustSubject {
            kind: SiderustSubjectKind::Star,
            body: SiderustBody::Sun,
            star_handle: handle,
            icrs_dir: SiderustSphericalDir {
                polar_deg: 0.0,
                azimuth_deg: 0.0,
                frame: SiderustFrame::ICRS,
            },
            generic_target_handle: ptr::null(),
        },
    )
}

fn generic_target_subject() -> (*mut SiderustGenericTarget, SiderustSubject) {
    let mut handle: *mut SiderustGenericTarget = ptr::null_mut();
    let status =
        unsafe { siderust_generic_target_create_icrs(83.82, -5.39, 2451545.0, &mut handle) };
    assert_eq!(status, SiderustStatus::Ok);
    (
        handle,
        SiderustSubject {
            kind: SiderustSubjectKind::GenericTarget,
            body: SiderustBody::Sun,
            star_handle: ptr::null(),
            icrs_dir: SiderustSphericalDir {
                polar_deg: 0.0,
                azimuth_deg: 0.0,
                frame: SiderustFrame::ICRS,
            },
            generic_target_handle: handle,
        },
    )
}

#[test]
fn star_lifecycle_still_works() {
    let (handle, _) = vega_subject();
    assert!(!handle.is_null());
    unsafe { siderust_star_free(handle) };
}

#[test]
fn generic_target_lifecycle_and_data_work() {
    let mut handle: *mut SiderustGenericTarget = ptr::null_mut();
    let status = unsafe {
        siderust_generic_target_create_icrs_with_pm(
            83.82,
            -5.39,
            2451545.0,
            0.0001,
            -0.0002,
            SiderustRaConvention::MuAlphaStar,
            &mut handle,
        )
    };
    assert_eq!(status, SiderustStatus::Ok);
    assert!(!handle.is_null());

    let mut data = SiderustGenericTargetData {
        kind: SiderustTargetCoordKind::CartesianPos,
        _pad1: [0; 4],
        coord: SiderustTargetCoordUnion {
            cartesian_pos: SiderustCartesianPos {
                x: 0.0,
                y: 0.0,
                z: 0.0,
                frame: SiderustFrame::ICRS,
                center: SiderustCenter::Geocentric,
                length_unit: SiderustLengthUnit::Meter,
            },
        },
        epoch_jd: TempochJd::new(0.0),
        has_proper_motion: false,
        _pad2: [0; 7],
        proper_motion: SiderustProperMotion {
            pm_ra_deg_yr: 0.0,
            pm_dec_deg_yr: 0.0,
            ra_convention: SiderustRaConvention::MuAlpha,
        },
    };
    let status = unsafe { siderust_generic_target_get_data(handle, &mut data) };
    assert_eq!(status, SiderustStatus::Ok);
    assert_eq!(data.kind, SiderustTargetCoordKind::SphericalDir);
    assert!(data.has_proper_motion);
    assert_eq!(data.epoch_jd.value, 2451545.0);

    unsafe { siderust_generic_target_free(handle) };
}

#[test]
fn subject_altitude_queries_cover_all_supported_kinds() {
    let mut out = 0.0_f64;

    let status =
        unsafe { siderust_altitude_at(sun_subject(), paris_observer(), mjd(60000.0), &mut out) };
    assert_eq!(status, SiderustStatus::Ok);
    assert!(out.is_finite());

    let (star_handle, star_subject) = vega_subject();
    let status =
        unsafe { siderust_altitude_at(star_subject, paris_observer(), mjd(60000.0), &mut out) };
    assert_eq!(status, SiderustStatus::Ok);
    assert!(out.is_finite());
    unsafe { siderust_star_free(star_handle) };

    let status =
        unsafe { siderust_altitude_at(icrs_subject(), paris_observer(), mjd(60000.0), &mut out) };
    assert_eq!(status, SiderustStatus::Ok);
    assert!(out.is_finite());

    let (target_handle, target_subject) = generic_target_subject();
    let status =
        unsafe { siderust_altitude_at(target_subject, paris_observer(), mjd(60000.0), &mut out) };
    assert_eq!(status, SiderustStatus::Ok);
    assert!(out.is_finite());
    unsafe { siderust_generic_target_free(target_handle) };
}

#[test]
fn subject_vector_queries_allocate_and_free() {
    let mut crossings: *mut SiderustCrossingEvent = ptr::null_mut();
    let mut count: usize = 0;
    let status = unsafe {
        siderust_crossings(
            sun_subject(),
            paris_observer(),
            one_day_window(),
            0.0,
            default_opts(),
            &mut crossings,
            &mut count,
        )
    };
    assert_eq!(status, SiderustStatus::Ok);
    unsafe { siderust_crossings_free(crossings, count) };

    let mut culminations: *mut SiderustCulminationEvent = ptr::null_mut();
    let status = unsafe {
        siderust_culminations(
            sun_subject(),
            paris_observer(),
            one_day_window(),
            default_opts(),
            &mut culminations,
            &mut count,
        )
    };
    assert_eq!(status, SiderustStatus::Ok);
    unsafe { siderust_culminations_free(culminations, count) };

    let (star_handle, star_subject) = vega_subject();
    let mut periods: *mut TempochPeriodMjd = ptr::null_mut();
    let status = unsafe {
        siderust_above_threshold(
            star_subject,
            paris_observer(),
            one_day_window(),
            0.0,
            default_opts(),
            &mut periods,
            &mut count,
        )
    };
    assert_eq!(status, SiderustStatus::Ok);
    unsafe {
        siderust_periods_free(periods, count);
        siderust_star_free(star_handle);
    }

    let (target_handle, target_subject) = generic_target_subject();
    let mut az_crossings: *mut SiderustAzimuthCrossingEvent = ptr::null_mut();
    let status = unsafe {
        siderust_azimuth_crossings(
            target_subject,
            paris_observer(),
            one_day_window(),
            180.0,
            default_opts(),
            &mut az_crossings,
            &mut count,
        )
    };
    assert_eq!(status, SiderustStatus::Ok);
    unsafe {
        siderust_azimuth_crossings_free(az_crossings, count);
        siderust_generic_target_free(target_handle);
    }

    let mut extrema: *mut SiderustAzimuthExtremum = ptr::null_mut();
    let status = unsafe {
        siderust_azimuth_extrema(
            moon_subject(),
            paris_observer(),
            one_day_window(),
            default_opts(),
            &mut extrema,
            &mut count,
        )
    };
    assert_eq!(status, SiderustStatus::Ok);
    unsafe { siderust_azimuth_extrema_free(extrema, count) };
}

#[test]
fn altitude_periods_stays_body_only() {
    let query = SiderustAltitudeQuery {
        observer: paris_observer(),
        start_mjd: mjd(60000.0),
        end_mjd: mjd(60001.0),
        min_altitude_deg: -90.0,
        max_altitude_deg: 90.0,
    };

    let mut out: *mut TempochPeriodMjd = ptr::null_mut();
    let mut count = 0usize;
    let status = unsafe { siderust_altitude_periods(sun_subject(), query, &mut out, &mut count) };
    assert_eq!(status, SiderustStatus::Ok);
    unsafe { siderust_periods_free(out, count) };

    let (target_handle, target_subject) = generic_target_subject();
    let status = unsafe { siderust_altitude_periods(target_subject, query, &mut out, &mut count) };
    assert_eq!(status, SiderustStatus::InvalidArgument);
    unsafe { siderust_generic_target_free(target_handle) };
}

#[test]
fn phase_events_alloc_and_free() {
    let window = TempochPeriodMjd {
        start_mjd: mjd(60000.0),
        end_mjd: mjd(60030.0),
    };
    let mut out: *mut SiderustPhaseEvent = ptr::null_mut();
    let mut count: usize = 0;
    let status =
        unsafe { siderust_find_phase_events(window, default_opts(), &mut out, &mut count) };
    assert_eq!(status, SiderustStatus::Ok);
    assert!(count >= 2);
    unsafe { siderust_phase_events_free(out, count) };
}

#[test]
fn null_pointer_guards_use_subject_handles() {
    let status = unsafe {
        siderust_altitude_at(
            sun_subject(),
            paris_observer(),
            mjd(60000.0),
            ptr::null_mut(),
        )
    };
    assert_eq!(status, SiderustStatus::NullPointer);

    let status = unsafe {
        siderust_crossings(
            sun_subject(),
            paris_observer(),
            one_day_window(),
            0.0,
            default_opts(),
            ptr::null_mut(),
            ptr::null_mut(),
        )
    };
    assert_eq!(status, SiderustStatus::NullPointer);

    let mut out: f64 = 0.0;
    let status = unsafe {
        siderust_altitude_at(
            SiderustSubject {
                kind: SiderustSubjectKind::Star,
                body: SiderustBody::Sun,
                star_handle: ptr::null(),
                icrs_dir: SiderustSphericalDir {
                    polar_deg: 0.0,
                    azimuth_deg: 0.0,
                    frame: SiderustFrame::ICRS,
                },
                generic_target_handle: ptr::null(),
            },
            paris_observer(),
            mjd(60000.0),
            &mut out,
        )
    };
    assert_eq!(status, SiderustStatus::NullPointer);

    let status = unsafe {
        siderust_altitude_at(
            SiderustSubject {
                kind: SiderustSubjectKind::GenericTarget,
                body: SiderustBody::Sun,
                star_handle: ptr::null(),
                icrs_dir: SiderustSphericalDir {
                    polar_deg: 0.0,
                    azimuth_deg: 0.0,
                    frame: SiderustFrame::ICRS,
                },
                generic_target_handle: ptr::null(),
            },
            paris_observer(),
            mjd(60000.0),
            &mut out,
        )
    };
    assert_eq!(status, SiderustStatus::NullPointer);
}

#[test]
fn invalid_period_is_still_reported() {
    let bad = TempochPeriodMjd {
        start_mjd: mjd(60001.0),
        end_mjd: mjd(60000.0),
    };
    let mut out: *mut TempochPeriodMjd = ptr::null_mut();
    let mut count = 0usize;
    let status = unsafe {
        siderust_above_threshold(
            sun_subject(),
            paris_observer(),
            bad,
            0.0,
            default_opts(),
            &mut out,
            &mut count,
        )
    };
    assert_eq!(status, SiderustStatus::InvalidPeriod);
}
