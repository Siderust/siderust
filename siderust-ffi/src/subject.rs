// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! **Unified subject dispatch** — one set of FFI functions that operate on any
//! [`SiderustSubject`].
//!
//! Instead of separate `siderust_sun_*`, `siderust_moon_*`, `siderust_body_*`,
//! `siderust_star_*`, `siderust_icrs_*`, and `siderust_target_*` families, the
//! caller constructs a [`SiderustSubject`] tagged value and passes it to these
//! generic entry-points.
//!
//! The old per-type functions remain for backwards compatibility but new code
//! should prefer these.

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
// Altitude — instantaneous
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
// Altitude — threshold periods
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
// Altitude — events
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
/// Only `Body` subjects support this operation.  For `Star`, `Icrs`, and
/// `Target`, `SIDERUST_STATUS_T_INVALID_ARGUMENT` is returned.
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
// Azimuth — instantaneous
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
// Azimuth — events
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
    use crate::bodies::SiderustStar;
    use crate::target::SiderustTarget;
    use crate::test_helpers::*;
    use std::ffi::CString;
    use std::ptr;

    // ── Subject constructors ──────────────────────────────────────────────

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
            target_handle: ptr::null(),
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
            target_handle: ptr::null(),
        }
    }

    fn mars_subject() -> SiderustSubject {
        SiderustSubject {
            kind: SiderustSubjectKind::Body,
            body: SiderustBody::Mars,
            star_handle: ptr::null(),
            icrs_dir: SiderustSphericalDir {
                polar_deg: 0.0,
                azimuth_deg: 0.0,
                frame: SiderustFrame::ICRS,
            },
            target_handle: ptr::null(),
        }
    }

    fn star_subject(name: &str) -> (*mut SiderustStar, SiderustSubject) {
        let cname = CString::new(name).unwrap();
        let mut handle: *mut SiderustStar = ptr::null_mut();
        let st = crate::bodies::siderust_star_catalog(cname.as_ptr(), &mut handle);
        assert_eq!(st, SiderustStatus::Ok);
        let subject = SiderustSubject {
            kind: SiderustSubjectKind::Star,
            body: SiderustBody::Sun, // unused
            star_handle: handle,
            icrs_dir: SiderustSphericalDir {
                polar_deg: 0.0,
                azimuth_deg: 0.0,
                frame: SiderustFrame::ICRS,
            },
            target_handle: ptr::null(),
        };
        (handle, subject)
    }

    fn icrs_vega_subject() -> SiderustSubject {
        SiderustSubject {
            kind: SiderustSubjectKind::Icrs,
            body: SiderustBody::Sun, // unused
            star_handle: ptr::null(),
            icrs_dir: SiderustSphericalDir {
                polar_deg: 38.78,
                azimuth_deg: 279.23,
                frame: SiderustFrame::ICRS,
            },
            target_handle: ptr::null(),
        }
    }

    fn target_subject(ra: f64, dec: f64) -> (*mut SiderustTarget, SiderustSubject) {
        let mut handle: *mut SiderustTarget = ptr::null_mut();
        let st = crate::target::siderust_target_create(ra, dec, 2451545.0, &mut handle);
        assert_eq!(st, SiderustStatus::Ok);
        let subject = SiderustSubject {
            kind: SiderustSubjectKind::Target,
            body: SiderustBody::Sun, // unused
            star_handle: ptr::null(),
            icrs_dir: SiderustSphericalDir {
                polar_deg: 0.0,
                azimuth_deg: 0.0,
                frame: SiderustFrame::ICRS,
            },
            target_handle: handle,
        };
        (handle, subject)
    }

    // ── altitude_at ───────────────────────────────────────────────────────

    #[test]
    fn altitude_at_sun() {
        let mut out = 0.0f64;
        let st = siderust_altitude_at(sun_subject(), paris(), 60000.5, &mut out);
        assert_eq!(st, SiderustStatus::Ok);
        assert!(out.is_finite());
    }

    #[test]
    fn altitude_at_moon() {
        let mut out = 0.0f64;
        let st = siderust_altitude_at(moon_subject(), paris(), 60000.5, &mut out);
        assert_eq!(st, SiderustStatus::Ok);
        assert!(out.is_finite());
    }

    #[test]
    fn altitude_at_mars() {
        let mut out = 0.0f64;
        let st = siderust_altitude_at(mars_subject(), paris(), 60000.5, &mut out);
        assert_eq!(st, SiderustStatus::Ok);
        assert!(out.is_finite());
    }

    #[test]
    fn altitude_at_star() {
        let (handle, subj) = star_subject("VEGA");
        let mut out = 0.0f64;
        let st = siderust_altitude_at(subj, paris(), 60000.5, &mut out);
        assert_eq!(st, SiderustStatus::Ok);
        assert!(out.is_finite());
        unsafe { crate::bodies::siderust_star_free(handle) };
    }

    #[test]
    fn altitude_at_icrs() {
        let mut out = 0.0f64;
        let st = siderust_altitude_at(icrs_vega_subject(), paris(), 60000.5, &mut out);
        assert_eq!(st, SiderustStatus::Ok);
        assert!(out.is_finite());
    }

    #[test]
    fn altitude_at_target() {
        let (handle, subj) = target_subject(279.23, 38.78);
        let mut out = 0.0f64;
        let st = siderust_altitude_at(subj, paris(), 60000.5, &mut out);
        assert_eq!(st, SiderustStatus::Ok);
        assert!(out.is_finite());
        unsafe { crate::target::siderust_target_free(handle) };
    }

    #[test]
    fn altitude_at_null_out() {
        let st = siderust_altitude_at(sun_subject(), paris(), 60000.5, ptr::null_mut());
        assert_eq!(st, SiderustStatus::NullPointer);
    }

    #[test]
    fn altitude_at_null_star_handle() {
        let subj = SiderustSubject {
            kind: SiderustSubjectKind::Star,
            body: SiderustBody::Sun,
            star_handle: ptr::null(),
            icrs_dir: SiderustSphericalDir {
                polar_deg: 0.0,
                azimuth_deg: 0.0,
                frame: SiderustFrame::ICRS,
            },
            target_handle: ptr::null(),
        };
        let mut out = 0.0f64;
        let st = siderust_altitude_at(subj, paris(), 60000.5, &mut out);
        assert_eq!(st, SiderustStatus::NullPointer);
    }

    // ── above_threshold ───────────────────────────────────────────────────

    #[test]
    fn above_threshold_sun() {
        let mut out: *mut TempochPeriodMjd = ptr::null_mut();
        let mut count = 0usize;
        let st = siderust_above_threshold(
            sun_subject(),
            paris(),
            one_day_window(),
            0.0,
            default_opts(),
            &mut out,
            &mut count,
        );
        assert_eq!(st, SiderustStatus::Ok);
        assert!(count > 0);
        unsafe { crate::altitude::siderust_periods_free(out, count) };
    }

    #[test]
    fn above_threshold_star() {
        let (handle, subj) = star_subject("VEGA");
        let mut out: *mut TempochPeriodMjd = ptr::null_mut();
        let mut count = 0usize;
        let st = siderust_above_threshold(
            subj,
            paris(),
            one_day_window(),
            0.0,
            default_opts(),
            &mut out,
            &mut count,
        );
        assert_eq!(st, SiderustStatus::Ok);
        unsafe { crate::altitude::siderust_periods_free(out, count) };
        unsafe { crate::bodies::siderust_star_free(handle) };
    }

    #[test]
    fn above_threshold_icrs() {
        let mut out: *mut TempochPeriodMjd = ptr::null_mut();
        let mut count = 0usize;
        let st = siderust_above_threshold(
            icrs_vega_subject(),
            paris(),
            one_day_window(),
            0.0,
            default_opts(),
            &mut out,
            &mut count,
        );
        assert_eq!(st, SiderustStatus::Ok);
        unsafe { crate::altitude::siderust_periods_free(out, count) };
    }

    #[test]
    fn above_threshold_target() {
        let (handle, subj) = target_subject(279.23, 38.78);
        let mut out: *mut TempochPeriodMjd = ptr::null_mut();
        let mut count = 0usize;
        let st = siderust_above_threshold(
            subj,
            paris(),
            one_day_window(),
            0.0,
            default_opts(),
            &mut out,
            &mut count,
        );
        assert_eq!(st, SiderustStatus::Ok);
        unsafe { crate::altitude::siderust_periods_free(out, count) };
        unsafe { crate::target::siderust_target_free(handle) };
    }

    // ── below_threshold ───────────────────────────────────────────────────

    #[test]
    fn below_threshold_sun() {
        let mut out: *mut TempochPeriodMjd = ptr::null_mut();
        let mut count = 0usize;
        let st = siderust_below_threshold(
            sun_subject(),
            paris(),
            one_day_window(),
            0.0,
            default_opts(),
            &mut out,
            &mut count,
        );
        assert_eq!(st, SiderustStatus::Ok);
        assert!(count > 0);
        unsafe { crate::altitude::siderust_periods_free(out, count) };
    }

    // ── crossings ─────────────────────────────────────────────────────────

    #[test]
    fn crossings_sun() {
        let mut out: *mut SiderustCrossingEvent = ptr::null_mut();
        let mut count = 0usize;
        let st = siderust_crossings(
            sun_subject(),
            paris(),
            one_day_window(),
            0.0,
            default_opts(),
            &mut out,
            &mut count,
        );
        assert_eq!(st, SiderustStatus::Ok);
        unsafe { crate::altitude::siderust_crossings_free(out, count) };
    }

    #[test]
    fn crossings_star() {
        let (handle, subj) = star_subject("VEGA");
        let mut out: *mut SiderustCrossingEvent = ptr::null_mut();
        let mut count = 0usize;
        let st = siderust_crossings(
            subj,
            paris(),
            one_day_window(),
            0.0,
            default_opts(),
            &mut out,
            &mut count,
        );
        assert_eq!(st, SiderustStatus::Ok);
        unsafe { crate::altitude::siderust_crossings_free(out, count) };
        unsafe { crate::bodies::siderust_star_free(handle) };
    }

    // ── culminations ──────────────────────────────────────────────────────

    #[test]
    fn culminations_sun() {
        let mut out: *mut SiderustCulminationEvent = ptr::null_mut();
        let mut count = 0usize;
        let st = siderust_culminations(
            sun_subject(),
            paris(),
            one_day_window(),
            default_opts(),
            &mut out,
            &mut count,
        );
        assert_eq!(st, SiderustStatus::Ok);
        unsafe { crate::altitude::siderust_culminations_free(out, count) };
    }

    #[test]
    fn culminations_target() {
        let (handle, subj) = target_subject(279.23, 38.78);
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
        unsafe { crate::altitude::siderust_culminations_free(out, count) };
        unsafe { crate::target::siderust_target_free(handle) };
    }

    // ── altitude_periods (body-only) ──────────────────────────────────────

    #[test]
    fn altitude_periods_sun() {
        let query = SiderustAltitudeQuery {
            observer: paris(),
            start_mjd: 60000.0,
            end_mjd: 60001.0,
            min_altitude_deg: -90.0,
            max_altitude_deg: 90.0,
        };
        let mut out: *mut TempochPeriodMjd = ptr::null_mut();
        let mut count = 0usize;
        let st = siderust_altitude_periods(sun_subject(), query, &mut out, &mut count);
        assert_eq!(st, SiderustStatus::Ok);
        unsafe { crate::altitude::siderust_periods_free(out, count) };
    }

    #[test]
    fn altitude_periods_star_unsupported() {
        let (handle, subj) = star_subject("VEGA");
        let query = SiderustAltitudeQuery {
            observer: paris(),
            start_mjd: 60000.0,
            end_mjd: 60001.0,
            min_altitude_deg: -90.0,
            max_altitude_deg: 90.0,
        };
        let mut out: *mut TempochPeriodMjd = ptr::null_mut();
        let mut count = 0usize;
        let st = siderust_altitude_periods(subj, query, &mut out, &mut count);
        assert_eq!(st, SiderustStatus::InvalidArgument);
        unsafe { crate::bodies::siderust_star_free(handle) };
    }

    // ── azimuth_at ────────────────────────────────────────────────────────

    #[test]
    fn azimuth_at_sun() {
        let mut out = 0.0f64;
        let st = siderust_azimuth_at(sun_subject(), paris(), 60000.5, &mut out);
        assert_eq!(st, SiderustStatus::Ok);
        assert!(out.is_finite());
    }

    #[test]
    fn azimuth_at_star() {
        let (handle, subj) = star_subject("VEGA");
        let mut out = 0.0f64;
        let st = siderust_azimuth_at(subj, paris(), 60000.5, &mut out);
        assert_eq!(st, SiderustStatus::Ok);
        assert!(out.is_finite());
        unsafe { crate::bodies::siderust_star_free(handle) };
    }

    #[test]
    fn azimuth_at_icrs() {
        let mut out = 0.0f64;
        let st = siderust_azimuth_at(icrs_vega_subject(), paris(), 60000.5, &mut out);
        assert_eq!(st, SiderustStatus::Ok);
        assert!(out.is_finite());
    }

    #[test]
    fn azimuth_at_target() {
        let (handle, subj) = target_subject(279.23, 38.78);
        let mut out = 0.0f64;
        let st = siderust_azimuth_at(subj, paris(), 60000.5, &mut out);
        assert_eq!(st, SiderustStatus::Ok);
        assert!(out.is_finite());
        unsafe { crate::target::siderust_target_free(handle) };
    }

    // ── azimuth_crossings ─────────────────────────────────────────────────

    #[test]
    fn azimuth_crossings_sun() {
        let mut out: *mut SiderustAzimuthCrossingEvent = ptr::null_mut();
        let mut count = 0usize;
        let st = siderust_azimuth_crossings(
            sun_subject(),
            paris(),
            one_day_window(),
            180.0,
            default_opts(),
            &mut out,
            &mut count,
        );
        assert_eq!(st, SiderustStatus::Ok);
        unsafe { crate::azimuth::siderust_azimuth_crossings_free(out, count) };
    }

    // ── azimuth_extrema ───────────────────────────────────────────────────

    #[test]
    fn azimuth_extrema_sun() {
        let mut out: *mut SiderustAzimuthExtremum = ptr::null_mut();
        let mut count = 0usize;
        let st = siderust_azimuth_extrema(
            sun_subject(),
            paris(),
            one_day_window(),
            default_opts(),
            &mut out,
            &mut count,
        );
        assert_eq!(st, SiderustStatus::Ok);
        unsafe { crate::azimuth::siderust_azimuth_extrema_free(out, count) };
    }

    // ── in_azimuth_range ──────────────────────────────────────────────────

    #[test]
    fn in_azimuth_range_sun() {
        let mut out: *mut TempochPeriodMjd = ptr::null_mut();
        let mut count = 0usize;
        let st = siderust_in_azimuth_range(
            sun_subject(),
            paris(),
            one_day_window(),
            90.0,
            270.0,
            default_opts(),
            &mut out,
            &mut count,
        );
        assert_eq!(st, SiderustStatus::Ok);
        unsafe { crate::altitude::siderust_periods_free(out, count) };
    }
}
