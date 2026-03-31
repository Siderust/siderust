// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! **Unified subject dispatch**, the **preferred FFI entry point** for altitude
//! and azimuth computations.
//!
//! # Overview
//!
//! Callers construct a [`SiderustSubject`] tagged value and pass it to these
//! generic entry-points.
//!
//! # Supported Subject Kinds
//!
//! - `Body`: Solar-system bodies (Sun, Moon, planets)
//! - `Star`: Catalog or custom stars via `SiderustStar` handle
//! - `Icrs`: Fixed ICRS direction (RA/Dec)
//! - `GenericTarget`: Full `CoordinateWithPM<T>` via [`SiderustGenericTarget`] handle

use crate::altitude::{crossings_to_c, culminations_to_c, periods_to_c, window_from_c};
use crate::azimuth::{vec_az_crossings_to_c, vec_az_extrema_to_c};
use crate::error::SiderustStatus;
use crate::types::*;
use qtty::*;
use siderust::calculus::azimuth::{azimuth_crossings, azimuth_extrema, in_azimuth_range};
use siderust::AltitudePeriodsProvider;
use siderust::AzimuthProvider;
use tempoch::ModifiedJulianDate;

// ═══════════════════════════════════════════════════════════════════════════
// Altitude, instantaneous
// ═══════════════════════════════════════════════════════════════════════════

/// Altitude of any subject at an instant (radians).
#[no_mangle]
pub extern "C" fn siderust_altitude_at(
    subject: SiderustSubject,
    observer: SiderustGeodetict,
    mjd: f64,
    out_rad: *mut f64,
) -> SiderustStatus {
    ffi_guard! {{
        if out_rad.is_null() {
            return SiderustStatus::NullPointer;
        }
        dispatch_subject!(subject, |p| {
            unsafe {
                *out_rad = p
                    .altitude_at(&observer.to_rust(), ModifiedJulianDate::new(mjd))
                    .value();
            }
            SiderustStatus::Ok
        })
    }}
}

// ═══════════════════════════════════════════════════════════════════════════
// Altitude, threshold periods
// ═══════════════════════════════════════════════════════════════════════════

/// Periods when a subject is above a threshold altitude.
#[no_mangle]
pub extern "C" fn siderust_above_threshold(
    subject: SiderustSubject,
    observer: SiderustGeodetict,
    window: TempochPeriodMjd,
    threshold_deg: f64,
    opts: SiderustSearchOpts,
    out: *mut *mut TempochPeriodMjd,
    count: *mut usize,
) -> SiderustStatus {
    ffi_guard! {{
        let window = match window_from_c(window) {
            Ok(w) => w,
            Err(e) => return e,
        };
        dispatch_subject!(subject, |p| {
            periods_to_c(
                siderust::above_threshold(
                    p,
                    &observer.to_rust(),
                    window,
                    Degrees::new(threshold_deg),
                    opts.to_rust(),
                ),
                out,
                count,
            )
        })
    }}
}

/// Periods when a subject is below a threshold altitude.
#[no_mangle]
pub extern "C" fn siderust_below_threshold(
    subject: SiderustSubject,
    observer: SiderustGeodetict,
    window: TempochPeriodMjd,
    threshold_deg: f64,
    opts: SiderustSearchOpts,
    out: *mut *mut TempochPeriodMjd,
    count: *mut usize,
) -> SiderustStatus {
    ffi_guard! {{
        let window = match window_from_c(window) {
            Ok(w) => w,
            Err(e) => return e,
        };
        dispatch_subject!(subject, |p| {
            periods_to_c(
                siderust::below_threshold(
                    p,
                    &observer.to_rust(),
                    window,
                    Degrees::new(threshold_deg),
                    opts.to_rust(),
                ),
                out,
                count,
            )
        })
    }}
}

// ═══════════════════════════════════════════════════════════════════════════
// Altitude, events
// ═══════════════════════════════════════════════════════════════════════════

/// Threshold-crossing events for a subject.
#[no_mangle]
pub extern "C" fn siderust_crossings(
    subject: SiderustSubject,
    observer: SiderustGeodetict,
    window: TempochPeriodMjd,
    threshold_deg: f64,
    opts: SiderustSearchOpts,
    out: *mut *mut SiderustCrossingEvent,
    count: *mut usize,
) -> SiderustStatus {
    ffi_guard! {{
        let window = match window_from_c(window) {
            Ok(w) => w,
            Err(e) => return e,
        };
        dispatch_subject!(subject, |p| {
            crossings_to_c(
                siderust::crossings(
                    p,
                    &observer.to_rust(),
                    window,
                    Degrees::new(threshold_deg),
                    opts.to_rust(),
                ),
                out,
                count,
            )
        })
    }}
}

/// Culmination (local extrema) events for a subject.
#[no_mangle]
pub extern "C" fn siderust_culminations(
    subject: SiderustSubject,
    observer: SiderustGeodetict,
    window: TempochPeriodMjd,
    opts: SiderustSearchOpts,
    out: *mut *mut SiderustCulminationEvent,
    count: *mut usize,
) -> SiderustStatus {
    ffi_guard! {{
        let window = match window_from_c(window) {
            Ok(w) => w,
            Err(e) => return e,
        };
        dispatch_subject!(subject, |p| {
            culminations_to_c(
                siderust::culminations(p, &observer.to_rust(), window, opts.to_rust()),
                out,
                count,
            )
        })
    }}
}

// ═══════════════════════════════════════════════════════════════════════════
// Altitude periods (body-only)
// ═══════════════════════════════════════════════════════════════════════════

/// Periods when a body's altitude is within [min, max].
///
/// Only `Body` subjects support this operation. For `Star`, `Icrs`, and
/// `GenericTarget`, `SIDERUST_STATUS_T_INVALID_ARGUMENT` is returned.
#[no_mangle]
pub extern "C" fn siderust_altitude_periods(
    subject: SiderustSubject,
    query: SiderustAltitudeQuery,
    out: *mut *mut TempochPeriodMjd,
    count: *mut usize,
) -> SiderustStatus {
    ffi_guard! {{
        match subject.kind {
            SiderustSubjectKind::Body => {
                let q = query.to_rust();
                dispatch_body!(subject.body, |b| {
                    periods_to_c(b.altitude_periods(&q), out, count)
                })
            }
            _ => SiderustStatus::InvalidArgument,
        }
    }}
}

// ═══════════════════════════════════════════════════════════════════════════
// Azimuth, instantaneous
// ═══════════════════════════════════════════════════════════════════════════

/// Azimuth of any subject at an instant (degrees, North-clockwise).
#[no_mangle]
pub extern "C" fn siderust_azimuth_at(
    subject: SiderustSubject,
    observer: SiderustGeodetict,
    mjd: f64,
    out_deg: *mut f64,
) -> SiderustStatus {
    ffi_guard! {{
        if out_deg.is_null() {
            return SiderustStatus::NullPointer;
        }
        dispatch_subject!(subject, |p| {
            unsafe {
                *out_deg = p
                    .azimuth_at(&observer.to_rust(), ModifiedJulianDate::new(mjd))
                    .value();
            }
            SiderustStatus::Ok
        })
    }}
}

// ═══════════════════════════════════════════════════════════════════════════
// Azimuth, events
// ═══════════════════════════════════════════════════════════════════════════

/// Azimuth bearing-crossing events for a subject.
#[no_mangle]
pub extern "C" fn siderust_azimuth_crossings(
    subject: SiderustSubject,
    observer: SiderustGeodetict,
    window: TempochPeriodMjd,
    bearing_deg: f64,
    opts: SiderustSearchOpts,
    out: *mut *mut SiderustAzimuthCrossingEvent,
    count: *mut usize,
) -> SiderustStatus {
    ffi_guard! {{
        let window = match window_from_c(window) {
            Ok(w) => w,
            Err(e) => return e,
        };
        dispatch_subject!(subject, |p| {
            vec_az_crossings_to_c(
                azimuth_crossings(
                    p,
                    &observer.to_rust(),
                    window,
                    Degrees::new(bearing_deg),
                    opts.to_rust(),
                ),
                out,
                count,
            )
        })
    }}
}

/// Azimuth extrema (northernmost / southernmost bearing) for a subject.
#[no_mangle]
pub extern "C" fn siderust_azimuth_extrema(
    subject: SiderustSubject,
    observer: SiderustGeodetict,
    window: TempochPeriodMjd,
    opts: SiderustSearchOpts,
    out: *mut *mut SiderustAzimuthExtremum,
    count: *mut usize,
) -> SiderustStatus {
    ffi_guard! {{
        let window = match window_from_c(window) {
            Ok(w) => w,
            Err(e) => return e,
        };
        dispatch_subject!(subject, |p| {
            vec_az_extrema_to_c(
                azimuth_extrema(p, &observer.to_rust(), window, opts.to_rust()),
                out,
                count,
            )
        })
    }}
}

/// Periods when a subject's azimuth is within [min_deg, max_deg].
#[no_mangle]
pub extern "C" fn siderust_in_azimuth_range(
    subject: SiderustSubject,
    observer: SiderustGeodetict,
    window: TempochPeriodMjd,
    min_deg: f64,
    max_deg: f64,
    opts: SiderustSearchOpts,
    out: *mut *mut TempochPeriodMjd,
    count: *mut usize,
) -> SiderustStatus {
    ffi_guard! {{
        let window = match window_from_c(window) {
            Ok(w) => w,
            Err(e) => return e,
        };
        dispatch_subject!(subject, |p| {
            periods_to_c(
                in_azimuth_range(
                    p,
                    &observer.to_rust(),
                    window,
                    Degrees::new(min_deg),
                    Degrees::new(max_deg),
                    opts.to_rust(),
                ),
                out,
                count,
            )
        })
    }}
}

// ═══════════════════════════════════════════════════════════════════════════
// Tests
// ═══════════════════════════════════════════════════════════════════════════

#[cfg(test)]
mod tests {
    use super::*;
    use crate::target::{
        siderust_generic_target_create_icrs, siderust_generic_target_create_icrs_with_pm,
        siderust_generic_target_free, SiderustGenericTarget,
    };
    use crate::test_helpers::*;
    use std::ffi::CString;
    use std::ptr;

    fn star_subject(name: &str) -> (*mut crate::bodies::SiderustStar, SiderustSubject) {
        let cname = CString::new(name).unwrap();
        let mut handle: *mut crate::bodies::SiderustStar = ptr::null_mut();
        let st = crate::bodies::siderust_star_catalog(cname.as_ptr(), &mut handle);
        assert_eq!(st, SiderustStatus::Ok);
        (handle, SiderustSubject::star(handle))
    }

    fn icrs_vega_subject() -> SiderustSubject {
        SiderustSubject::icrs(icrs_vega())
    }

    fn generic_target_subject() -> (*mut SiderustGenericTarget, SiderustSubject) {
        let mut handle: *mut SiderustGenericTarget = ptr::null_mut();
        let st = siderust_generic_target_create_icrs(279.23, 38.78, 2451545.0, &mut handle);
        assert_eq!(st, SiderustStatus::Ok);
        (handle, SiderustSubject::generic_target(handle))
    }

    #[test]
    fn altitude_at_body_star_icrs_and_generic_target() {
        let mut out = 0.0f64;
        let st = siderust_altitude_at(
            SiderustSubject::body(SiderustBody::Sun),
            paris(),
            mjd(60000.5),
            &mut out,
        );
        assert_eq!(st, SiderustStatus::Ok);
        assert!(out.is_finite());

        let (star_handle, star_subj) = star_subject("VEGA");
        let st = siderust_altitude_at(star_subj, paris(), mjd(60000.5), &mut out);
        assert_eq!(st, SiderustStatus::Ok);
        assert!(out.is_finite());

        let st = siderust_altitude_at(icrs_vega_subject(), paris(), mjd(60000.5), &mut out);
        assert_eq!(st, SiderustStatus::Ok);
        assert!(out.is_finite());

        let (target_handle, target_subj) = generic_target_subject();
        let st = siderust_altitude_at(target_subj, paris(), mjd(60000.5), &mut out);
        assert_eq!(st, SiderustStatus::Ok);
        assert!(out.is_finite());

        unsafe {
            crate::bodies::siderust_star_free(star_handle);
            siderust_generic_target_free(target_handle);
        }
    }

    #[test]
    fn altitude_at_null_out() {
        let st = siderust_altitude_at(
            SiderustSubject::body(SiderustBody::Sun),
            paris(),
            mjd(60000.5),
            ptr::null_mut(),
        );
        assert_eq!(st, SiderustStatus::NullPointer);
    }

    #[test]
    fn altitude_at_null_star_handle() {
        let subj = SiderustSubject::star(ptr::null());
        let mut out = 0.0f64;
        let st = siderust_altitude_at(subj, paris(), mjd(60000.5), &mut out);
        assert_eq!(st, SiderustStatus::NullPointer);
    }

    #[test]
    fn altitude_at_null_generic_target_handle() {
        let subj = SiderustSubject::generic_target(ptr::null());
        let mut out = 0.0f64;
        let st = siderust_altitude_at(subj, paris(), mjd(60000.5), &mut out);
        assert_eq!(st, SiderustStatus::NullPointer);
    }

    #[test]
    fn above_threshold_star_and_generic_target_with_pm() {
        let (star_handle, star_subj) = star_subject("VEGA");
        let mut out: *mut TempochPeriodMjd = ptr::null_mut();
        let mut count = 0usize;
        let st = siderust_above_threshold(
            star_subj,
            paris(),
            one_day_window(),
            0.0,
            default_opts(),
            &mut out,
            &mut count,
        );
        assert_eq!(st, SiderustStatus::Ok);
        unsafe { crate::altitude::siderust_periods_free(out, count) };

        let mut handle: *mut SiderustGenericTarget = ptr::null_mut();
        let st = siderust_generic_target_create_icrs_with_pm(
            279.23,
            38.78,
            2451545.0,
            0.0001,
            -0.0001,
            SiderustRaConvention::MuAlphaStar,
            &mut handle,
        );
        assert_eq!(st, SiderustStatus::Ok);
        let st = siderust_above_threshold(
            SiderustSubject::generic_target(handle),
            paris(),
            one_day_window(),
            0.0,
            default_opts(),
            &mut out,
            &mut count,
        );
        assert_eq!(st, SiderustStatus::Ok);
        unsafe {
            crate::altitude::siderust_periods_free(out, count);
            crate::bodies::siderust_star_free(star_handle);
            siderust_generic_target_free(handle);
        }
    }

    #[test]
    fn below_threshold_and_crossings_body() {
        let mut out: *mut TempochPeriodMjd = ptr::null_mut();
        let mut count = 0usize;
        let st = siderust_below_threshold(
            SiderustSubject::body(SiderustBody::Sun),
            paris(),
            one_day_window(),
            0.0,
            default_opts(),
            &mut out,
            &mut count,
        );
        assert_eq!(st, SiderustStatus::Ok);
        unsafe { crate::altitude::siderust_periods_free(out, count) };
        let mut crossings: *mut SiderustCrossingEvent = ptr::null_mut();
        let st = siderust_crossings(
            SiderustSubject::body(SiderustBody::Sun),
            paris(),
            one_day_window(),
            0.0,
            default_opts(),
            &mut crossings,
            &mut count,
        );
        assert_eq!(st, SiderustStatus::Ok);
        unsafe { crate::altitude::siderust_crossings_free(crossings, count) };
    }

    #[test]
    fn culminations_generic_target() {
        let (handle, subj) = generic_target_subject();
        let mut out: *mut SiderustCulminationEvent = ptr::null_mut();
        let mut count = 0usize;
        let st = siderust_culminations(
            subj,
            paris(),
            one_day_window(),
            default_opts(),
            &mut out,
            &mut count,
        );
        assert_eq!(st, SiderustStatus::Ok);
        unsafe {
            crate::altitude::siderust_culminations_free(out, count);
            siderust_generic_target_free(handle);
        }
    }

    #[test]
    fn altitude_periods_body_only() {
        let query = SiderustAltitudeQuery {
            observer: paris(),
            start_mjd: mjd(60000.0),
            end_mjd: mjd(60001.0),
            min_altitude_deg: -90.0,
            max_altitude_deg: 90.0,
        };
        let mut out: *mut TempochPeriodMjd = ptr::null_mut();
        let mut count = 0usize;
        let st = siderust_altitude_periods(
            SiderustSubject::body(SiderustBody::Sun),
            query,
            &mut out,
            &mut count,
        );
        assert_eq!(st, SiderustStatus::Ok);
        unsafe { crate::altitude::siderust_periods_free(out, count) };
        let (star_handle, star_subj) = star_subject("VEGA");
        let st = siderust_altitude_periods(star_subj, query, &mut out, &mut count);
        assert_eq!(st, SiderustStatus::InvalidArgument);
        unsafe { crate::bodies::siderust_star_free(star_handle) };

        let (target_handle, target_subj) = generic_target_subject();
        let st = siderust_altitude_periods(target_subj, query, &mut out, &mut count);
        assert_eq!(st, SiderustStatus::InvalidArgument);
        unsafe { siderust_generic_target_free(target_handle) };
    }

    #[test]
    fn azimuth_subject_queries() {
        let mut out = 0.0f64;
        let st = siderust_azimuth_at(icrs_vega_subject(), paris(), mjd(60000.5), &mut out);
        assert_eq!(st, SiderustStatus::Ok);
        assert!(out.is_finite());

        let mut crossings: *mut SiderustAzimuthCrossingEvent = ptr::null_mut();
        let mut count = 0usize;
        let st = siderust_azimuth_crossings(
            SiderustSubject::body(SiderustBody::Moon),
            paris(),
            one_day_window(),
            180.0,
            default_opts(),
            &mut crossings,
            &mut count,
        );
        assert_eq!(st, SiderustStatus::Ok);
        unsafe { crate::azimuth::siderust_azimuth_crossings_free(crossings, count) };

        let mut extrema: *mut SiderustAzimuthExtremum = ptr::null_mut();
        let st = siderust_azimuth_extrema(
            SiderustSubject::body(SiderustBody::Sun),
            paris(),
            one_day_window(),
            default_opts(),
            &mut extrema,
            &mut count,
        );
        assert_eq!(st, SiderustStatus::Ok);
        unsafe { crate::azimuth::siderust_azimuth_extrema_free(extrema, count) };

        let mut periods: *mut TempochPeriodMjd = ptr::null_mut();
        let st = siderust_in_azimuth_range(
            SiderustSubject::body(SiderustBody::Sun),
            paris(),
            one_day_window(),
            90.0,
            270.0,
            default_opts(),
            &mut periods,
            &mut count,
        );
        assert_eq!(st, SiderustStatus::Ok);
        unsafe { crate::altitude::siderust_periods_free(periods, count) };
    }
}
