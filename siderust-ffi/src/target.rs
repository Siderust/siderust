// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Target FFI — opaque handle for an ICRS celestial target.
//!
//! A `SiderustTarget` wraps an ICRS direction (RA/Dec) and exposes altitude
//! and azimuth queries by delegating to the unified [`crate::subject`] module.

use crate::error::SiderustStatus;
use crate::types::*;
use qtty::*;
use siderust::coordinates::spherical;

// ═══════════════════════════════════════════════════════════════════════════
// Opaque handle
// ═══════════════════════════════════════════════════════════════════════════

/// Opaque handle representing an ICRS celestial target (RA/Dec).
pub struct SiderustTarget {
    pub(crate) dir: spherical::direction::ICRS,
    pub(crate) ra_deg: f64,
    pub(crate) dec_deg: f64,
    pub(crate) epoch_jd: f64,
}

// ═══════════════════════════════════════════════════════════════════════════
// Lifecycle
// ═══════════════════════════════════════════════════════════════════════════

/// Create a new target from right ascension and declination (degrees) and epoch
/// (Julian Date).
///
/// On success, writes a newly-allocated pointer to `*out`.
/// The caller must free it with [`siderust_target_free`].
#[no_mangle]
pub extern "C" fn siderust_target_create(
    ra_deg: f64,
    dec_deg: f64,
    epoch_jd: f64,
    out: *mut *mut SiderustTarget,
) -> SiderustStatus {
    ffi_guard! {{
        if out.is_null() {
            return SiderustStatus::NullPointer;
        }
        let dir = spherical::direction::ICRS::new(Degrees::new(ra_deg), Degrees::new(dec_deg));
        let handle = Box::new(SiderustTarget {
            dir,
            ra_deg,
            dec_deg,
            epoch_jd,
        });
        unsafe { *out = Box::into_raw(handle) };
        SiderustStatus::Ok
    }}
}

/// Free a target handle created by [`siderust_target_create`].
///
/// # Safety
/// The handle must not be used after this call.
#[no_mangle]
pub unsafe extern "C" fn siderust_target_free(handle: *mut SiderustTarget) {
    if !handle.is_null() {
        drop(Box::from_raw(handle));
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// Property accessors
// ═══════════════════════════════════════════════════════════════════════════

/// Write the right ascension (degrees) to `*out`.
#[no_mangle]
pub extern "C" fn siderust_target_ra_deg(
    handle: *const SiderustTarget,
    out: *mut f64,
) -> SiderustStatus {
    ffi_guard! {{
        if handle.is_null() || out.is_null() {
            return SiderustStatus::NullPointer;
        }
        unsafe { *out = (*handle).ra_deg };
        SiderustStatus::Ok
    }}
}

/// Write the declination (degrees) to `*out`.
#[no_mangle]
pub extern "C" fn siderust_target_dec_deg(
    handle: *const SiderustTarget,
    out: *mut f64,
) -> SiderustStatus {
    ffi_guard! {{
        if handle.is_null() || out.is_null() {
            return SiderustStatus::NullPointer;
        }
        unsafe { *out = (*handle).dec_deg };
        SiderustStatus::Ok
    }}
}

/// Write the reference epoch (Julian Date) to `*out`.
#[no_mangle]
pub extern "C" fn siderust_target_epoch_jd(
    handle: *const SiderustTarget,
    out: *mut f64,
) -> SiderustStatus {
    ffi_guard! {{
        if handle.is_null() || out.is_null() {
            return SiderustStatus::NullPointer;
        }
        unsafe { *out = (*handle).epoch_jd };
        SiderustStatus::Ok
    }}
}

// ═══════════════════════════════════════════════════════════════════════════
// Altitude  — thin wrappers
// ═══════════════════════════════════════════════════════════════════════════

/// Altitude of the target (radians).
#[no_mangle]
pub extern "C" fn siderust_target_altitude_at(
    handle: *const SiderustTarget,
    observer: SiderustGeodetict,
    mjd: f64,
    out_rad: *mut f64,
) -> SiderustStatus {
    crate::subject::siderust_altitude_at(SiderustSubject::target(handle), observer, mjd, out_rad)
}

/// Periods when the target is above `threshold_deg`.
#[no_mangle]
pub extern "C" fn siderust_target_above_threshold(
    handle: *const SiderustTarget,
    observer: SiderustGeodetict,
    window: TempochPeriodMjd,
    threshold_deg: f64,
    opts: SiderustSearchOpts,
    out: *mut *mut TempochPeriodMjd,
    count: *mut usize,
) -> SiderustStatus {
    crate::subject::siderust_above_threshold(
        SiderustSubject::target(handle),
        observer,
        window,
        threshold_deg,
        opts,
        out,
        count,
    )
}

/// Periods when the target is below `threshold_deg`.
#[no_mangle]
pub extern "C" fn siderust_target_below_threshold(
    handle: *const SiderustTarget,
    observer: SiderustGeodetict,
    window: TempochPeriodMjd,
    threshold_deg: f64,
    opts: SiderustSearchOpts,
    out: *mut *mut TempochPeriodMjd,
    count: *mut usize,
) -> SiderustStatus {
    crate::subject::siderust_below_threshold(
        SiderustSubject::target(handle),
        observer,
        window,
        threshold_deg,
        opts,
        out,
        count,
    )
}

/// Altitude crossing events for the target.
#[no_mangle]
pub extern "C" fn siderust_target_crossings(
    handle: *const SiderustTarget,
    observer: SiderustGeodetict,
    window: TempochPeriodMjd,
    threshold_deg: f64,
    opts: SiderustSearchOpts,
    out: *mut *mut SiderustCrossingEvent,
    count: *mut usize,
) -> SiderustStatus {
    crate::subject::siderust_crossings(
        SiderustSubject::target(handle),
        observer,
        window,
        threshold_deg,
        opts,
        out,
        count,
    )
}

/// Culmination events for the target.
#[no_mangle]
pub extern "C" fn siderust_target_culminations(
    handle: *const SiderustTarget,
    observer: SiderustGeodetict,
    window: TempochPeriodMjd,
    opts: SiderustSearchOpts,
    out: *mut *mut SiderustCulminationEvent,
    count: *mut usize,
) -> SiderustStatus {
    crate::subject::siderust_culminations(
        SiderustSubject::target(handle),
        observer,
        window,
        opts,
        out,
        count,
    )
}

// ═══════════════════════════════════════════════════════════════════════════
// Azimuth  — thin wrappers
// ═══════════════════════════════════════════════════════════════════════════

/// Azimuth of the target (degrees, North-clockwise).
#[no_mangle]
pub extern "C" fn siderust_target_azimuth_at(
    handle: *const SiderustTarget,
    observer: SiderustGeodetict,
    mjd: f64,
    out_deg: *mut f64,
) -> SiderustStatus {
    crate::subject::siderust_azimuth_at(SiderustSubject::target(handle), observer, mjd, out_deg)
}

/// Azimuth bearing-crossing events for the target.
#[no_mangle]
pub extern "C" fn siderust_target_azimuth_crossings(
    handle: *const SiderustTarget,
    observer: SiderustGeodetict,
    window: TempochPeriodMjd,
    bearing_deg: f64,
    opts: SiderustSearchOpts,
    out: *mut *mut SiderustAzimuthCrossingEvent,
    count: *mut usize,
) -> SiderustStatus {
    crate::subject::siderust_azimuth_crossings(
        SiderustSubject::target(handle),
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
    use crate::test_helpers::*;
    use std::ptr;

    // Vega: RA = 279.23°, Dec = 38.78° (approximate)
    fn create_vega() -> *mut SiderustTarget {
        let mut h: *mut SiderustTarget = ptr::null_mut();
        assert_eq!(
            siderust_target_create(279.23, 38.78, 2451545.0, &mut h),
            SiderustStatus::Ok
        );
        assert!(!h.is_null());
        h
    }

    // ── Lifecycle ─────────────────────────────────────────────────────────

    #[test]
    fn create_null_out() {
        assert_eq!(
            siderust_target_create(0.0, 0.0, 0.0, ptr::null_mut()),
            SiderustStatus::NullPointer
        );
    }

    #[test]
    fn create_and_free() {
        let h = create_vega();
        unsafe { siderust_target_free(h) };
    }

    #[test]
    fn free_null_safe() {
        unsafe { siderust_target_free(ptr::null_mut()) };
    }

    // ── Accessors ─────────────────────────────────────────────────────────

    #[test]
    fn ra_deg_roundtrip() {
        let h = create_vega();
        let mut out = 0.0f64;
        assert_eq!(siderust_target_ra_deg(h, &mut out), SiderustStatus::Ok);
        assert!((out - 279.23).abs() < 1e-10);
        unsafe { siderust_target_free(h) };
    }

    #[test]
    fn dec_deg_roundtrip() {
        let h = create_vega();
        let mut out = 0.0f64;
        assert_eq!(siderust_target_dec_deg(h, &mut out), SiderustStatus::Ok);
        assert!((out - 38.78).abs() < 1e-10);
        unsafe { siderust_target_free(h) };
    }

    #[test]
    fn epoch_jd_roundtrip() {
        let h = create_vega();
        let mut out = 0.0f64;
        assert_eq!(siderust_target_epoch_jd(h, &mut out), SiderustStatus::Ok);
        assert!((out - 2451545.0).abs() < 1e-6);
        unsafe { siderust_target_free(h) };
    }

    #[test]
    fn ra_null_handle() {
        let mut out = 0.0f64;
        assert_eq!(
            siderust_target_ra_deg(ptr::null(), &mut out),
            SiderustStatus::NullPointer
        );
    }

    #[test]
    fn ra_null_out() {
        let h = create_vega();
        assert_eq!(
            siderust_target_ra_deg(h, ptr::null_mut()),
            SiderustStatus::NullPointer
        );
        unsafe { siderust_target_free(h) };
    }

    #[test]
    fn dec_null_out() {
        let h = create_vega();
        assert_eq!(
            siderust_target_dec_deg(h, ptr::null_mut()),
            SiderustStatus::NullPointer
        );
        unsafe { siderust_target_free(h) };
    }

    #[test]
    fn epoch_null_out() {
        let h = create_vega();
        assert_eq!(
            siderust_target_epoch_jd(h, ptr::null_mut()),
            SiderustStatus::NullPointer
        );
        unsafe { siderust_target_free(h) };
    }

    // ── Altitude ──────────────────────────────────────────────────────────

    #[test]
    fn altitude_at_is_finite() {
        let h = create_vega();
        let mut out = 0.0f64;
        assert_eq!(
            siderust_target_altitude_at(h, paris(), 60000.0, &mut out),
            SiderustStatus::Ok
        );
        assert!(out.is_finite());
        unsafe { siderust_target_free(h) };
    }

    #[test]
    fn altitude_at_null_handle() {
        let mut out = 0.0f64;
        assert_eq!(
            siderust_target_altitude_at(ptr::null(), paris(), 60000.0, &mut out),
            SiderustStatus::NullPointer
        );
    }

    #[test]
    fn altitude_at_null_out() {
        let h = create_vega();
        assert_eq!(
            siderust_target_altitude_at(h, paris(), 60000.0, ptr::null_mut()),
            SiderustStatus::NullPointer
        );
        unsafe { siderust_target_free(h) };
    }

    #[test]
    fn above_threshold_ok() {
        let h = create_vega();
        let mut out: *mut TempochPeriodMjd = ptr::null_mut();
        let mut count = 0usize;
        let s = siderust_target_above_threshold(
            h,
            paris(),
            one_day_window(),
            0.0,
            default_opts(),
            &mut out,
            &mut count,
        );
        assert_eq!(s, SiderustStatus::Ok);
        unsafe {
            crate::altitude::siderust_periods_free(out, count);
            siderust_target_free(h);
        }
    }

    #[test]
    fn above_threshold_null_handle() {
        let mut out: *mut TempochPeriodMjd = ptr::null_mut();
        let mut count = 0usize;
        assert_eq!(
            siderust_target_above_threshold(
                ptr::null(),
                paris(),
                one_day_window(),
                0.0,
                default_opts(),
                &mut out,
                &mut count
            ),
            SiderustStatus::NullPointer
        );
    }

    #[test]
    fn above_threshold_invalid_window() {
        let h = create_vega();
        let bad = TempochPeriodMjd {
            start_mjd: 60001.0,
            end_mjd: 60000.0,
        };
        let mut out: *mut TempochPeriodMjd = ptr::null_mut();
        let mut count = 0usize;
        assert_eq!(
            siderust_target_above_threshold(
                h,
                paris(),
                bad,
                0.0,
                default_opts(),
                &mut out,
                &mut count
            ),
            SiderustStatus::InvalidPeriod
        );
        unsafe { siderust_target_free(h) };
    }

    #[test]
    fn crossings_ok() {
        let h = create_vega();
        let mut out: *mut SiderustCrossingEvent = ptr::null_mut();
        let mut count = 0usize;
        let s = siderust_target_crossings(
            h,
            paris(),
            one_day_window(),
            0.0,
            default_opts(),
            &mut out,
            &mut count,
        );
        assert_eq!(s, SiderustStatus::Ok);
        unsafe {
            crate::altitude::siderust_crossings_free(out, count);
            siderust_target_free(h);
        }
    }

    #[test]
    fn crossings_null_handle() {
        let mut out: *mut SiderustCrossingEvent = ptr::null_mut();
        let mut count = 0usize;
        assert_eq!(
            siderust_target_crossings(
                ptr::null(),
                paris(),
                one_day_window(),
                0.0,
                default_opts(),
                &mut out,
                &mut count
            ),
            SiderustStatus::NullPointer
        );
    }

    #[test]
    fn culminations_ok() {
        let h = create_vega();
        let mut out: *mut SiderustCulminationEvent = ptr::null_mut();
        let mut count = 0usize;
        let s = siderust_target_culminations(
            h,
            paris(),
            one_day_window(),
            default_opts(),
            &mut out,
            &mut count,
        );
        assert_eq!(s, SiderustStatus::Ok);
        unsafe {
            crate::altitude::siderust_culminations_free(out, count);
            siderust_target_free(h);
        }
    }

    #[test]
    fn culminations_null_handle() {
        let mut out: *mut SiderustCulminationEvent = ptr::null_mut();
        let mut count = 0usize;
        assert_eq!(
            siderust_target_culminations(
                ptr::null(),
                paris(),
                one_day_window(),
                default_opts(),
                &mut out,
                &mut count
            ),
            SiderustStatus::NullPointer
        );
    }

    // ── Azimuth ───────────────────────────────────────────────────────────

    #[test]
    fn azimuth_at_is_finite() {
        let h = create_vega();
        let mut out = 0.0f64;
        assert_eq!(
            siderust_target_azimuth_at(h, paris(), 60000.0, &mut out),
            SiderustStatus::Ok
        );
        assert!(out.is_finite() && out >= 0.0 && out < 360.0);
        unsafe { siderust_target_free(h) };
    }

    #[test]
    fn azimuth_at_null_handle() {
        let mut out = 0.0f64;
        assert_eq!(
            siderust_target_azimuth_at(ptr::null(), paris(), 60000.0, &mut out),
            SiderustStatus::NullPointer
        );
    }

    #[test]
    fn azimuth_at_null_out() {
        let h = create_vega();
        assert_eq!(
            siderust_target_azimuth_at(h, paris(), 60000.0, ptr::null_mut()),
            SiderustStatus::NullPointer
        );
        unsafe { siderust_target_free(h) };
    }

    #[test]
    fn azimuth_crossings_ok() {
        let h = create_vega();
        let mut out: *mut SiderustAzimuthCrossingEvent = ptr::null_mut();
        let mut count = 0usize;
        let s = siderust_target_azimuth_crossings(
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
            crate::azimuth::siderust_azimuth_crossings_free(out, count);
            siderust_target_free(h);
        }
    }

    #[test]
    fn azimuth_crossings_null_handle() {
        let mut out: *mut SiderustAzimuthCrossingEvent = ptr::null_mut();
        let mut count = 0usize;
        assert_eq!(
            siderust_target_azimuth_crossings(
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

    #[test]
    fn azimuth_crossings_invalid_window() {
        let h = create_vega();
        let bad = TempochPeriodMjd {
            start_mjd: 60001.0,
            end_mjd: 60000.0,
        };
        let mut out: *mut SiderustAzimuthCrossingEvent = ptr::null_mut();
        let mut count = 0usize;
        assert_eq!(
            siderust_target_azimuth_crossings(
                h,
                paris(),
                bad,
                180.0,
                default_opts(),
                &mut out,
                &mut count
            ),
            SiderustStatus::InvalidPeriod
        );
        unsafe { siderust_target_free(h) };
    }
}
