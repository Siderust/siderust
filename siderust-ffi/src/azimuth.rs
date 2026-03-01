// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Azimuth FFI — bearing crossings, extrema, and range periods.
//!
//! Per-body functions are thin wrappers delegating to the unified
//! [`crate::subject`] functions.  Helpers and free functions are kept
//! here because `subject.rs` depends on them.

use crate::bodies::SiderustStar;
use crate::error::SiderustStatus;
use crate::ffi_utils::{free_boxed_slice, vec_to_c, FfiFrom};
use crate::types::*;
use siderust::calculus::azimuth::{AzimuthCrossingEvent, AzimuthExtremum};

// ═══════════════════════════════════════════════════════════════════════════
// Conversion helpers
// ═══════════════════════════════════════════════════════════════════════════

pub(crate) fn vec_az_crossings_to_c(
    events: Vec<AzimuthCrossingEvent>,
    out: *mut *mut SiderustAzimuthCrossingEvent,
    count: *mut usize,
) -> SiderustStatus {
    vec_to_c(events, SiderustAzimuthCrossingEvent::ffi_from, out, count)
}

pub(crate) fn vec_az_extrema_to_c(
    events: Vec<AzimuthExtremum>,
    out: *mut *mut SiderustAzimuthExtremum,
    count: *mut usize,
) -> SiderustStatus {
    vec_to_c(events, SiderustAzimuthExtremum::ffi_from, out, count)
}

// ═══════════════════════════════════════════════════════════════════════════
// Memory management
// ═══════════════════════════════════════════════════════════════════════════

/// Free an array of azimuth crossing events.
#[no_mangle]
pub unsafe extern "C" fn siderust_azimuth_crossings_free(
    ptr: *mut SiderustAzimuthCrossingEvent,
    count: usize,
) {
    free_boxed_slice(ptr, count);
}

/// Free an array of azimuth extrema.
#[no_mangle]
pub unsafe extern "C" fn siderust_azimuth_extrema_free(
    ptr: *mut SiderustAzimuthExtremum,
    count: usize,
) {
    free_boxed_slice(ptr, count);
}

// ═══════════════════════════════════════════════════════════════════════════
// Sun azimuth  — thin wrappers
// ═══════════════════════════════════════════════════════════════════════════

/// Azimuth of the Sun at an instant (degrees, North-clockwise).
#[no_mangle]
pub extern "C" fn siderust_sun_azimuth_at(
    observer: SiderustGeodetict,
    mjd: f64,
    out_deg: *mut f64,
) -> SiderustStatus {
    crate::subject::siderust_azimuth_at(
        SiderustSubject::body(SiderustBody::Sun),
        observer,
        mjd,
        out_deg,
    )
}

/// Azimuth bearing-crossing events for the Sun.
#[no_mangle]
pub extern "C" fn siderust_sun_azimuth_crossings(
    observer: SiderustGeodetict,
    window: TempochPeriodMjd,
    bearing_deg: f64,
    opts: SiderustSearchOpts,
    out: *mut *mut SiderustAzimuthCrossingEvent,
    count: *mut usize,
) -> SiderustStatus {
    crate::subject::siderust_azimuth_crossings(
        SiderustSubject::body(SiderustBody::Sun),
        observer,
        window,
        bearing_deg,
        opts,
        out,
        count,
    )
}

/// Azimuth extrema (northernmost and southernmost bearing) for the Sun.
#[no_mangle]
pub extern "C" fn siderust_sun_azimuth_extrema(
    observer: SiderustGeodetict,
    window: TempochPeriodMjd,
    opts: SiderustSearchOpts,
    out: *mut *mut SiderustAzimuthExtremum,
    count: *mut usize,
) -> SiderustStatus {
    crate::subject::siderust_azimuth_extrema(
        SiderustSubject::body(SiderustBody::Sun),
        observer,
        window,
        opts,
        out,
        count,
    )
}

/// Periods when the Sun's azimuth is within [min_deg, max_deg].
#[no_mangle]
pub extern "C" fn siderust_sun_in_azimuth_range(
    observer: SiderustGeodetict,
    window: TempochPeriodMjd,
    min_deg: f64,
    max_deg: f64,
    opts: SiderustSearchOpts,
    out: *mut *mut TempochPeriodMjd,
    count: *mut usize,
) -> SiderustStatus {
    crate::subject::siderust_in_azimuth_range(
        SiderustSubject::body(SiderustBody::Sun),
        observer,
        window,
        min_deg,
        max_deg,
        opts,
        out,
        count,
    )
}

// ═══════════════════════════════════════════════════════════════════════════
// Moon azimuth  — thin wrappers
// ═══════════════════════════════════════════════════════════════════════════

/// Azimuth of the Moon at an instant (degrees, North-clockwise).
#[no_mangle]
pub extern "C" fn siderust_moon_azimuth_at(
    observer: SiderustGeodetict,
    mjd: f64,
    out_deg: *mut f64,
) -> SiderustStatus {
    crate::subject::siderust_azimuth_at(
        SiderustSubject::body(SiderustBody::Moon),
        observer,
        mjd,
        out_deg,
    )
}

/// Azimuth bearing-crossing events for the Moon.
#[no_mangle]
pub extern "C" fn siderust_moon_azimuth_crossings(
    observer: SiderustGeodetict,
    window: TempochPeriodMjd,
    bearing_deg: f64,
    opts: SiderustSearchOpts,
    out: *mut *mut SiderustAzimuthCrossingEvent,
    count: *mut usize,
) -> SiderustStatus {
    crate::subject::siderust_azimuth_crossings(
        SiderustSubject::body(SiderustBody::Moon),
        observer,
        window,
        bearing_deg,
        opts,
        out,
        count,
    )
}

/// Azimuth extrema for the Moon.
#[no_mangle]
pub extern "C" fn siderust_moon_azimuth_extrema(
    observer: SiderustGeodetict,
    window: TempochPeriodMjd,
    opts: SiderustSearchOpts,
    out: *mut *mut SiderustAzimuthExtremum,
    count: *mut usize,
) -> SiderustStatus {
    crate::subject::siderust_azimuth_extrema(
        SiderustSubject::body(SiderustBody::Moon),
        observer,
        window,
        opts,
        out,
        count,
    )
}

/// Periods when the Moon's azimuth is within [min_deg, max_deg].
#[no_mangle]
pub extern "C" fn siderust_moon_in_azimuth_range(
    observer: SiderustGeodetict,
    window: TempochPeriodMjd,
    min_deg: f64,
    max_deg: f64,
    opts: SiderustSearchOpts,
    out: *mut *mut TempochPeriodMjd,
    count: *mut usize,
) -> SiderustStatus {
    crate::subject::siderust_in_azimuth_range(
        SiderustSubject::body(SiderustBody::Moon),
        observer,
        window,
        min_deg,
        max_deg,
        opts,
        out,
        count,
    )
}

// ═══════════════════════════════════════════════════════════════════════════
// Star azimuth  — thin wrappers
// ═══════════════════════════════════════════════════════════════════════════

/// Azimuth of a star at an instant (degrees, North-clockwise).
#[no_mangle]
pub extern "C" fn siderust_star_azimuth_at(
    handle: *const SiderustStar,
    observer: SiderustGeodetict,
    mjd: f64,
    out_deg: *mut f64,
) -> SiderustStatus {
    crate::subject::siderust_azimuth_at(SiderustSubject::star(handle), observer, mjd, out_deg)
}

/// Azimuth bearing-crossing events for a star.
#[no_mangle]
pub extern "C" fn siderust_star_azimuth_crossings(
    handle: *const SiderustStar,
    observer: SiderustGeodetict,
    window: TempochPeriodMjd,
    bearing_deg: f64,
    opts: SiderustSearchOpts,
    out: *mut *mut SiderustAzimuthCrossingEvent,
    count: *mut usize,
) -> SiderustStatus {
    crate::subject::siderust_azimuth_crossings(
        SiderustSubject::star(handle),
        observer,
        window,
        bearing_deg,
        opts,
        out,
        count,
    )
}

/// Periods when a star's azimuth is within [min_deg, max_deg].
#[no_mangle]
pub extern "C" fn siderust_star_in_azimuth_range(
    handle: *const SiderustStar,
    observer: SiderustGeodetict,
    window: TempochPeriodMjd,
    min_deg: f64,
    max_deg: f64,
    opts: SiderustSearchOpts,
    out: *mut *mut TempochPeriodMjd,
    count: *mut usize,
) -> SiderustStatus {
    crate::subject::siderust_in_azimuth_range(
        SiderustSubject::star(handle),
        observer,
        window,
        min_deg,
        max_deg,
        opts,
        out,
        count,
    )
}

// ═══════════════════════════════════════════════════════════════════════════
// ICRS direction azimuth  — thin wrappers
// ═══════════════════════════════════════════════════════════════════════════

/// Azimuth of an ICRS direction at an instant.
#[no_mangle]
pub extern "C" fn siderust_icrs_azimuth_at(
    dir: SiderustSphericalDir,
    observer: SiderustGeodetict,
    mjd: f64,
    out_deg: *mut f64,
) -> SiderustStatus {
    if dir.frame != SiderustFrame::ICRS {
        return SiderustStatus::InvalidFrame;
    }
    crate::subject::siderust_azimuth_at(SiderustSubject::icrs(dir), observer, mjd, out_deg)
}

/// Azimuth bearing-crossing events for an ICRS direction.
#[no_mangle]
pub extern "C" fn siderust_icrs_azimuth_crossings(
    dir: SiderustSphericalDir,
    observer: SiderustGeodetict,
    window: TempochPeriodMjd,
    bearing_deg: f64,
    opts: SiderustSearchOpts,
    out: *mut *mut SiderustAzimuthCrossingEvent,
    count: *mut usize,
) -> SiderustStatus {
    if dir.frame != SiderustFrame::ICRS {
        return SiderustStatus::InvalidFrame;
    }
    crate::subject::siderust_azimuth_crossings(
        SiderustSubject::icrs(dir),
        observer,
        window,
        bearing_deg,
        opts,
        out,
        count,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bodies::SiderustStar;
    use crate::test_helpers::*;
    use std::ffi::CString;
    use std::ptr;

    fn get_vega() -> *mut SiderustStar {
        let cname = CString::new("VEGA").unwrap();
        let mut h: *mut SiderustStar = ptr::null_mut();
        crate::bodies::siderust_star_catalog(cname.as_ptr(), &mut h);
        h
    }

    // ── Sun azimuth ───────────────────────────────────────────────────────

    #[test]
    fn sun_azimuth_at_is_finite() {
        let mut out = 0.0f64;
        assert_eq!(
            siderust_sun_azimuth_at(paris(), 60000.0, &mut out),
            SiderustStatus::Ok
        );
        assert!(out.is_finite() && out >= 0.0 && out < 360.0);
    }

    #[test]
    fn sun_azimuth_at_null_out() {
        assert_eq!(
            siderust_sun_azimuth_at(paris(), 60000.0, ptr::null_mut()),
            SiderustStatus::NullPointer
        );
    }

    #[test]
    fn sun_azimuth_crossings_ok() {
        let mut out: *mut SiderustAzimuthCrossingEvent = ptr::null_mut();
        let mut count = 0usize;
        let s = siderust_sun_azimuth_crossings(
            paris(),
            one_day_window(),
            180.0,
            default_opts(),
            &mut out,
            &mut count,
        );
        assert_eq!(s, SiderustStatus::Ok);
        unsafe { siderust_azimuth_crossings_free(out, count) };
    }

    #[test]
    fn sun_azimuth_extrema_ok() {
        let mut out: *mut SiderustAzimuthExtremum = ptr::null_mut();
        let mut count = 0usize;
        let s = siderust_sun_azimuth_extrema(
            paris(),
            one_day_window(),
            default_opts(),
            &mut out,
            &mut count,
        );
        assert_eq!(s, SiderustStatus::Ok);
        unsafe { siderust_azimuth_extrema_free(out, count) };
    }

    #[test]
    fn sun_in_azimuth_range_ok() {
        let mut out: *mut TempochPeriodMjd = ptr::null_mut();
        let mut count = 0usize;
        let s = siderust_sun_in_azimuth_range(
            paris(),
            one_day_window(),
            90.0,
            270.0,
            default_opts(),
            &mut out,
            &mut count,
        );
        assert_eq!(s, SiderustStatus::Ok);
        unsafe { crate::altitude::siderust_periods_free(out, count) };
    }

    // ── Moon azimuth ──────────────────────────────────────────────────────

    #[test]
    fn moon_azimuth_at_is_finite() {
        let mut out = 0.0f64;
        assert_eq!(
            siderust_moon_azimuth_at(paris(), 60000.0, &mut out),
            SiderustStatus::Ok
        );
        assert!(out.is_finite());
    }

    #[test]
    fn moon_azimuth_at_null_out() {
        assert_eq!(
            siderust_moon_azimuth_at(paris(), 60000.0, ptr::null_mut()),
            SiderustStatus::NullPointer
        );
    }

    #[test]
    fn moon_azimuth_crossings_ok() {
        let mut out: *mut SiderustAzimuthCrossingEvent = ptr::null_mut();
        let mut count = 0usize;
        let s = siderust_moon_azimuth_crossings(
            paris(),
            one_day_window(),
            180.0,
            default_opts(),
            &mut out,
            &mut count,
        );
        assert_eq!(s, SiderustStatus::Ok);
        unsafe { siderust_azimuth_crossings_free(out, count) };
    }

    #[test]
    fn moon_azimuth_extrema_ok() {
        let mut out: *mut SiderustAzimuthExtremum = ptr::null_mut();
        let mut count = 0usize;
        let s = siderust_moon_azimuth_extrema(
            paris(),
            one_day_window(),
            default_opts(),
            &mut out,
            &mut count,
        );
        assert_eq!(s, SiderustStatus::Ok);
        unsafe { siderust_azimuth_extrema_free(out, count) };
    }

    #[test]
    fn moon_in_azimuth_range_ok() {
        let mut out: *mut TempochPeriodMjd = ptr::null_mut();
        let mut count = 0usize;
        let s = siderust_moon_in_azimuth_range(
            paris(),
            one_day_window(),
            90.0,
            270.0,
            default_opts(),
            &mut out,
            &mut count,
        );
        assert_eq!(s, SiderustStatus::Ok);
        unsafe { crate::altitude::siderust_periods_free(out, count) };
    }

    // ── Star azimuth ──────────────────────────────────────────────────────

    #[test]
    fn star_azimuth_at_ok() {
        let h = get_vega();
        let mut out = 0.0f64;
        assert_eq!(
            siderust_star_azimuth_at(h, paris(), 60000.0, &mut out),
            SiderustStatus::Ok
        );
        assert!(out.is_finite());
        unsafe { crate::bodies::siderust_star_free(h) };
    }

    #[test]
    fn star_azimuth_at_null_handle() {
        let mut out = 0.0f64;
        assert_eq!(
            siderust_star_azimuth_at(ptr::null(), paris(), 60000.0, &mut out),
            SiderustStatus::NullPointer
        );
    }

    #[test]
    fn star_azimuth_at_null_out() {
        let h = get_vega();
        assert_eq!(
            siderust_star_azimuth_at(h, paris(), 60000.0, ptr::null_mut()),
            SiderustStatus::NullPointer
        );
        unsafe { crate::bodies::siderust_star_free(h) };
    }

    #[test]
    fn star_azimuth_crossings_ok() {
        let h = get_vega();
        let mut out: *mut SiderustAzimuthCrossingEvent = ptr::null_mut();
        let mut count = 0usize;
        let s = siderust_star_azimuth_crossings(
            h,
            paris(),
            one_day_window(),
            180.0,
            default_opts(),
            &mut out,
            &mut count,
        );
        assert_eq!(s, SiderustStatus::Ok);
        unsafe {
            siderust_azimuth_crossings_free(out, count);
            crate::bodies::siderust_star_free(h);
        }
    }

    #[test]
    fn star_azimuth_crossings_null_handle() {
        let mut out: *mut SiderustAzimuthCrossingEvent = ptr::null_mut();
        let mut count = 0usize;
        assert_eq!(
            siderust_star_azimuth_crossings(
                ptr::null(),
                paris(),
                one_day_window(),
                180.0,
                default_opts(),
                &mut out,
                &mut count
            ),
            SiderustStatus::NullPointer
        );
    }

    // ── ICRS azimuth ──────────────────────────────────────────────────────

    #[test]
    fn icrs_azimuth_at_ok() {
        let mut out = 0.0f64;
        assert_eq!(
            siderust_icrs_azimuth_at(icrs_vega(), paris(), 60000.0, &mut out),
            SiderustStatus::Ok
        );
        assert!(out.is_finite());
    }

    #[test]
    fn icrs_azimuth_at_null_out() {
        assert_eq!(
            siderust_icrs_azimuth_at(icrs_vega(), paris(), 60000.0, ptr::null_mut()),
            SiderustStatus::NullPointer
        );
    }

    #[test]
    fn icrs_azimuth_at_wrong_frame() {
        let dir = SiderustSphericalDir {
            polar_deg: 0.0,
            azimuth_deg: 0.0,
            frame: SiderustFrame::EquatorialMeanJ2000,
        };
        let mut out = 0.0f64;
        assert_eq!(
            siderust_icrs_azimuth_at(dir, paris(), 60000.0, &mut out),
            SiderustStatus::InvalidFrame
        );
    }

    #[test]
    fn icrs_azimuth_crossings_ok() {
        let mut out: *mut SiderustAzimuthCrossingEvent = ptr::null_mut();
        let mut count = 0usize;
        let s = siderust_icrs_azimuth_crossings(
            icrs_vega(),
            paris(),
            one_day_window(),
            180.0,
            default_opts(),
            &mut out,
            &mut count,
        );
        assert_eq!(s, SiderustStatus::Ok);
        unsafe { siderust_azimuth_crossings_free(out, count) };
    }

    #[test]
    fn icrs_azimuth_crossings_wrong_frame() {
        let dir = SiderustSphericalDir {
            polar_deg: 0.0,
            azimuth_deg: 0.0,
            frame: SiderustFrame::Galactic,
        };
        let mut out: *mut SiderustAzimuthCrossingEvent = ptr::null_mut();
        let mut count = 0usize;
        assert_eq!(
            siderust_icrs_azimuth_crossings(
                dir,
                paris(),
                one_day_window(),
                180.0,
                default_opts(),
                &mut out,
                &mut count
            ),
            SiderustStatus::InvalidFrame
        );
    }

    // ── Invalid window ────────────────────────────────────────────────────

    #[test]
    fn sun_crossings_invalid_window() {
        let bad = TempochPeriodMjd {
            start_mjd: 60001.0,
            end_mjd: 60000.0,
        };
        let mut out: *mut SiderustAzimuthCrossingEvent = ptr::null_mut();
        let mut count = 0usize;
        assert_eq!(
            siderust_sun_azimuth_crossings(paris(), bad, 0.0, default_opts(), &mut out, &mut count),
            SiderustStatus::InvalidPeriod
        );
    }

    #[test]
    fn sun_extrema_invalid_window() {
        let bad = TempochPeriodMjd {
            start_mjd: 60001.0,
            end_mjd: 60000.0,
        };
        let mut out: *mut SiderustAzimuthExtremum = ptr::null_mut();
        let mut count = 0usize;
        assert_eq!(
            siderust_sun_azimuth_extrema(paris(), bad, default_opts(), &mut out, &mut count),
            SiderustStatus::InvalidPeriod
        );
    }

    #[test]
    fn sun_in_range_invalid_window() {
        let bad = TempochPeriodMjd {
            start_mjd: 60001.0,
            end_mjd: 60000.0,
        };
        let mut out: *mut TempochPeriodMjd = ptr::null_mut();
        let mut count = 0usize;
        assert_eq!(
            siderust_sun_in_azimuth_range(
                paris(),
                bad,
                0.0,
                180.0,
                default_opts(),
                &mut out,
                &mut count
            ),
            SiderustStatus::InvalidPeriod
        );
    }

    #[test]
    fn free_null_pointers_safe() {
        unsafe {
            siderust_azimuth_crossings_free(ptr::null_mut(), 0);
            siderust_azimuth_extrema_free(ptr::null_mut(), 0);
        }
    }
}
