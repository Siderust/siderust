// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! FFI bindings for the altitude / observation API.
//!
//! Provides per-body-type functions for altitude_at, above_threshold,
//! below_threshold, crossings, culminations, and altitude_periods.

use crate::bodies::SiderustStar;
use crate::error::SiderustStatus;
use crate::ffi_utils::{free_boxed_slice, vec_to_c, FfiFrom};
use crate::types::*;
use qtty::*;
use siderust::bodies::{Moon, Sun};
use siderust::coordinates::spherical;
use siderust::AltitudePeriodsProvider;
use tempoch::{Interval, ModifiedJulianDate, Period, MJD};

// ═══════════════════════════════════════════════════════════════════════════
// Helpers
// ═══════════════════════════════════════════════════════════════════════════

pub(crate) fn window_from_c(w: TempochPeriodMjd) -> Result<Period<MJD>, SiderustStatus> {
    if w.start_mjd > w.end_mjd {
        return Err(SiderustStatus::InvalidPeriod);
    }
    Ok(Interval::new(
        ModifiedJulianDate::new(w.start_mjd),
        ModifiedJulianDate::new(w.end_mjd),
    ))
}

pub(crate) fn periods_to_c(
    periods: Vec<Period<MJD>>,
    out: *mut *mut TempochPeriodMjd,
    count: *mut usize,
) -> SiderustStatus {
    vec_to_c(periods, TempochPeriodMjd::ffi_from, out, count)
}

pub(crate) fn crossings_to_c(
    events: Vec<siderust::CrossingEvent>,
    out: *mut *mut SiderustCrossingEvent,
    count: *mut usize,
) -> SiderustStatus {
    vec_to_c(events, SiderustCrossingEvent::ffi_from, out, count)
}

pub(crate) fn culminations_to_c(
    events: Vec<siderust::CulminationEvent>,
    out: *mut *mut SiderustCulminationEvent,
    count: *mut usize,
) -> SiderustStatus {
    vec_to_c(events, SiderustCulminationEvent::ffi_from, out, count)
}

pub(crate) fn icrs_from_c(
    dir: SiderustSphericalDir,
) -> Result<spherical::direction::ICRS, SiderustStatus> {
    if dir.frame != SiderustFrame::ICRS {
        return Err(SiderustStatus::InvalidFrame);
    }
    Ok(spherical::direction::ICRS::new(
        Degrees::new(dir.azimuth_deg),
        Degrees::new(dir.polar_deg),
    ))
}

// ═══════════════════════════════════════════════════════════════════════════
// Free functions for dynamic arrays
// ═══════════════════════════════════════════════════════════════════════════

/// Free an array of MJD periods.
#[no_mangle]
pub unsafe extern "C" fn siderust_periods_free(ptr: *mut TempochPeriodMjd, count: usize) {
    free_boxed_slice(ptr, count);
}

/// Free an array of crossing events.
#[no_mangle]
pub unsafe extern "C" fn siderust_crossings_free(ptr: *mut SiderustCrossingEvent, count: usize) {
    free_boxed_slice(ptr, count);
}

/// Free an array of culmination events.
#[no_mangle]
pub unsafe extern "C" fn siderust_culminations_free(
    ptr: *mut SiderustCulminationEvent,
    count: usize,
) {
    free_boxed_slice(ptr, count);
}

// ═══════════════════════════════════════════════════════════════════════════
// Sun
// ═══════════════════════════════════════════════════════════════════════════

/// Altitude of the Sun at an instant (radians).
#[no_mangle]
pub extern "C" fn siderust_sun_altitude_at(
    observer: SiderustGeodetict,
    mjd: f64,
    out_rad: *mut f64,
) -> SiderustStatus {
    if out_rad.is_null() {
        return SiderustStatus::NullPointer;
    }
    unsafe {
        *out_rad = Sun
            .altitude_at(&observer.to_rust(), ModifiedJulianDate::new(mjd))
            .value();
    }
    SiderustStatus::Ok
}

/// Periods when the Sun is above a threshold altitude.
#[no_mangle]
pub extern "C" fn siderust_sun_above_threshold(
    observer: SiderustGeodetict,
    window: TempochPeriodMjd,
    threshold_deg: f64,
    opts: SiderustSearchOpts,
    out: *mut *mut TempochPeriodMjd,
    count: *mut usize,
) -> SiderustStatus {
    let window = match window_from_c(window) {
        Ok(w) => w,
        Err(e) => return e,
    };
    periods_to_c(
        siderust::above_threshold(
            &Sun,
            &observer.to_rust(),
            window,
            Degrees::new(threshold_deg),
            opts.to_rust(),
        ),
        out,
        count,
    )
}

/// Periods when the Sun is below a threshold altitude.
#[no_mangle]
pub extern "C" fn siderust_sun_below_threshold(
    observer: SiderustGeodetict,
    window: TempochPeriodMjd,
    threshold_deg: f64,
    opts: SiderustSearchOpts,
    out: *mut *mut TempochPeriodMjd,
    count: *mut usize,
) -> SiderustStatus {
    let window = match window_from_c(window) {
        Ok(w) => w,
        Err(e) => return e,
    };
    periods_to_c(
        siderust::below_threshold(
            &Sun,
            &observer.to_rust(),
            window,
            Degrees::new(threshold_deg),
            opts.to_rust(),
        ),
        out,
        count,
    )
}

/// Threshold-crossing events for the Sun.
#[no_mangle]
pub extern "C" fn siderust_sun_crossings(
    observer: SiderustGeodetict,
    window: TempochPeriodMjd,
    threshold_deg: f64,
    opts: SiderustSearchOpts,
    out: *mut *mut SiderustCrossingEvent,
    count: *mut usize,
) -> SiderustStatus {
    let window = match window_from_c(window) {
        Ok(w) => w,
        Err(e) => return e,
    };
    crossings_to_c(
        siderust::crossings(
            &Sun,
            &observer.to_rust(),
            window,
            Degrees::new(threshold_deg),
            opts.to_rust(),
        ),
        out,
        count,
    )
}

/// Culmination (local extrema) events for the Sun.
#[no_mangle]
pub extern "C" fn siderust_sun_culminations(
    observer: SiderustGeodetict,
    window: TempochPeriodMjd,
    opts: SiderustSearchOpts,
    out: *mut *mut SiderustCulminationEvent,
    count: *mut usize,
) -> SiderustStatus {
    let window = match window_from_c(window) {
        Ok(w) => w,
        Err(e) => return e,
    };
    culminations_to_c(
        siderust::culminations(&Sun, &observer.to_rust(), window, opts.to_rust()),
        out,
        count,
    )
}

/// Periods when the Sun's altitude is within [min, max].
#[no_mangle]
pub extern "C" fn siderust_sun_altitude_periods(
    query: SiderustAltitudeQuery,
    out: *mut *mut TempochPeriodMjd,
    count: *mut usize,
) -> SiderustStatus {
    periods_to_c(Sun.altitude_periods(&query.to_rust()), out, count)
}

// ═══════════════════════════════════════════════════════════════════════════
// Moon
// ═══════════════════════════════════════════════════════════════════════════

/// Altitude of the Moon at an instant (radians).
#[no_mangle]
pub extern "C" fn siderust_moon_altitude_at(
    observer: SiderustGeodetict,
    mjd: f64,
    out_rad: *mut f64,
) -> SiderustStatus {
    if out_rad.is_null() {
        return SiderustStatus::NullPointer;
    }
    unsafe {
        *out_rad = Moon
            .altitude_at(&observer.to_rust(), ModifiedJulianDate::new(mjd))
            .value();
    }
    SiderustStatus::Ok
}

/// Periods when the Moon is above a threshold altitude.
#[no_mangle]
pub extern "C" fn siderust_moon_above_threshold(
    observer: SiderustGeodetict,
    window: TempochPeriodMjd,
    threshold_deg: f64,
    opts: SiderustSearchOpts,
    out: *mut *mut TempochPeriodMjd,
    count: *mut usize,
) -> SiderustStatus {
    let window = match window_from_c(window) {
        Ok(w) => w,
        Err(e) => return e,
    };
    periods_to_c(
        siderust::above_threshold(
            &Moon,
            &observer.to_rust(),
            window,
            Degrees::new(threshold_deg),
            opts.to_rust(),
        ),
        out,
        count,
    )
}

/// Periods when the Moon is below a threshold altitude.
#[no_mangle]
pub extern "C" fn siderust_moon_below_threshold(
    observer: SiderustGeodetict,
    window: TempochPeriodMjd,
    threshold_deg: f64,
    opts: SiderustSearchOpts,
    out: *mut *mut TempochPeriodMjd,
    count: *mut usize,
) -> SiderustStatus {
    let window = match window_from_c(window) {
        Ok(w) => w,
        Err(e) => return e,
    };
    periods_to_c(
        siderust::below_threshold(
            &Moon,
            &observer.to_rust(),
            window,
            Degrees::new(threshold_deg),
            opts.to_rust(),
        ),
        out,
        count,
    )
}

/// Threshold-crossing events for the Moon.
#[no_mangle]
pub extern "C" fn siderust_moon_crossings(
    observer: SiderustGeodetict,
    window: TempochPeriodMjd,
    threshold_deg: f64,
    opts: SiderustSearchOpts,
    out: *mut *mut SiderustCrossingEvent,
    count: *mut usize,
) -> SiderustStatus {
    let window = match window_from_c(window) {
        Ok(w) => w,
        Err(e) => return e,
    };
    crossings_to_c(
        siderust::crossings(
            &Moon,
            &observer.to_rust(),
            window,
            Degrees::new(threshold_deg),
            opts.to_rust(),
        ),
        out,
        count,
    )
}

/// Culmination (local extrema) events for the Moon.
#[no_mangle]
pub extern "C" fn siderust_moon_culminations(
    observer: SiderustGeodetict,
    window: TempochPeriodMjd,
    opts: SiderustSearchOpts,
    out: *mut *mut SiderustCulminationEvent,
    count: *mut usize,
) -> SiderustStatus {
    let window = match window_from_c(window) {
        Ok(w) => w,
        Err(e) => return e,
    };
    culminations_to_c(
        siderust::culminations(&Moon, &observer.to_rust(), window, opts.to_rust()),
        out,
        count,
    )
}

/// Periods when the Moon's altitude is within [min, max].
#[no_mangle]
pub extern "C" fn siderust_moon_altitude_periods(
    query: SiderustAltitudeQuery,
    out: *mut *mut TempochPeriodMjd,
    count: *mut usize,
) -> SiderustStatus {
    periods_to_c(Moon.altitude_periods(&query.to_rust()), out, count)
}

// ═══════════════════════════════════════════════════════════════════════════
// Star (opaque handle)
// ═══════════════════════════════════════════════════════════════════════════

/// Altitude of a star at an instant (radians).
#[no_mangle]
pub extern "C" fn siderust_star_altitude_at(
    handle: *const SiderustStar,
    observer: SiderustGeodetict,
    mjd: f64,
    out_rad: *mut f64,
) -> SiderustStatus {
    if handle.is_null() || out_rad.is_null() {
        return SiderustStatus::NullPointer;
    }
    let star = unsafe { &(*handle).inner };
    unsafe {
        *out_rad = star
            .altitude_at(&observer.to_rust(), ModifiedJulianDate::new(mjd))
            .value();
    }
    SiderustStatus::Ok
}

/// Periods when a star is above a threshold altitude.
#[no_mangle]
pub extern "C" fn siderust_star_above_threshold(
    handle: *const SiderustStar,
    observer: SiderustGeodetict,
    window: TempochPeriodMjd,
    threshold_deg: f64,
    opts: SiderustSearchOpts,
    out: *mut *mut TempochPeriodMjd,
    count: *mut usize,
) -> SiderustStatus {
    if handle.is_null() {
        return SiderustStatus::NullPointer;
    }
    let star = unsafe { &(*handle).inner };
    let window = match window_from_c(window) {
        Ok(w) => w,
        Err(e) => return e,
    };
    periods_to_c(
        siderust::above_threshold(
            star,
            &observer.to_rust(),
            window,
            Degrees::new(threshold_deg),
            opts.to_rust(),
        ),
        out,
        count,
    )
}

/// Periods when a star is below a threshold altitude.
#[no_mangle]
pub extern "C" fn siderust_star_below_threshold(
    handle: *const SiderustStar,
    observer: SiderustGeodetict,
    window: TempochPeriodMjd,
    threshold_deg: f64,
    opts: SiderustSearchOpts,
    out: *mut *mut TempochPeriodMjd,
    count: *mut usize,
) -> SiderustStatus {
    if handle.is_null() {
        return SiderustStatus::NullPointer;
    }
    let star = unsafe { &(*handle).inner };
    let window = match window_from_c(window) {
        Ok(w) => w,
        Err(e) => return e,
    };
    periods_to_c(
        siderust::below_threshold(
            star,
            &observer.to_rust(),
            window,
            Degrees::new(threshold_deg),
            opts.to_rust(),
        ),
        out,
        count,
    )
}

/// Threshold-crossing events for a star.
#[no_mangle]
pub extern "C" fn siderust_star_crossings(
    handle: *const SiderustStar,
    observer: SiderustGeodetict,
    window: TempochPeriodMjd,
    threshold_deg: f64,
    opts: SiderustSearchOpts,
    out: *mut *mut SiderustCrossingEvent,
    count: *mut usize,
) -> SiderustStatus {
    if handle.is_null() {
        return SiderustStatus::NullPointer;
    }
    let star = unsafe { &(*handle).inner };
    let window = match window_from_c(window) {
        Ok(w) => w,
        Err(e) => return e,
    };
    crossings_to_c(
        siderust::crossings(
            star,
            &observer.to_rust(),
            window,
            Degrees::new(threshold_deg),
            opts.to_rust(),
        ),
        out,
        count,
    )
}

/// Culmination events for a star.
#[no_mangle]
pub extern "C" fn siderust_star_culminations(
    handle: *const SiderustStar,
    observer: SiderustGeodetict,
    window: TempochPeriodMjd,
    opts: SiderustSearchOpts,
    out: *mut *mut SiderustCulminationEvent,
    count: *mut usize,
) -> SiderustStatus {
    if handle.is_null() {
        return SiderustStatus::NullPointer;
    }
    let star = unsafe { &(*handle).inner };
    let window = match window_from_c(window) {
        Ok(w) => w,
        Err(e) => return e,
    };
    culminations_to_c(
        siderust::culminations(star, &observer.to_rust(), window, opts.to_rust()),
        out,
        count,
    )
}

// ═══════════════════════════════════════════════════════════════════════════
// ICRS direction
// ═══════════════════════════════════════════════════════════════════════════

/// Altitude of an ICRS direction at an instant (radians).
#[no_mangle]
pub extern "C" fn siderust_icrs_altitude_at(
    dir: SiderustSphericalDir,
    observer: SiderustGeodetict,
    mjd: f64,
    out_rad: *mut f64,
) -> SiderustStatus {
    if out_rad.is_null() {
        return SiderustStatus::NullPointer;
    }
    let dir = match icrs_from_c(dir) {
        Ok(d) => d,
        Err(e) => return e,
    };
    unsafe {
        *out_rad = dir
            .altitude_at(&observer.to_rust(), ModifiedJulianDate::new(mjd))
            .value();
    }
    SiderustStatus::Ok
}

/// Periods when an ICRS direction is above a threshold altitude.
#[no_mangle]
pub extern "C" fn siderust_icrs_above_threshold(
    dir: SiderustSphericalDir,
    observer: SiderustGeodetict,
    window: TempochPeriodMjd,
    threshold_deg: f64,
    opts: SiderustSearchOpts,
    out: *mut *mut TempochPeriodMjd,
    count: *mut usize,
) -> SiderustStatus {
    let dir = match icrs_from_c(dir) {
        Ok(d) => d,
        Err(e) => return e,
    };
    let window = match window_from_c(window) {
        Ok(w) => w,
        Err(e) => return e,
    };
    periods_to_c(
        siderust::above_threshold(
            &dir,
            &observer.to_rust(),
            window,
            Degrees::new(threshold_deg),
            opts.to_rust(),
        ),
        out,
        count,
    )
}

/// Periods when an ICRS direction is below a threshold altitude.
#[no_mangle]
pub extern "C" fn siderust_icrs_below_threshold(
    dir: SiderustSphericalDir,
    observer: SiderustGeodetict,
    window: TempochPeriodMjd,
    threshold_deg: f64,
    opts: SiderustSearchOpts,
    out: *mut *mut TempochPeriodMjd,
    count: *mut usize,
) -> SiderustStatus {
    let dir = match icrs_from_c(dir) {
        Ok(d) => d,
        Err(e) => return e,
    };
    let window = match window_from_c(window) {
        Ok(w) => w,
        Err(e) => return e,
    };
    periods_to_c(
        siderust::below_threshold(
            &dir,
            &observer.to_rust(),
            window,
            Degrees::new(threshold_deg),
            opts.to_rust(),
        ),
        out,
        count,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bodies::SiderustStar;
    use std::ffi::CString;
    use std::ptr;

    fn paris() -> SiderustGeodetict {
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

    fn icrs_vega() -> SiderustSphericalDir {
        SiderustSphericalDir {
            polar_deg: 38.78,
            azimuth_deg: 279.23,
            frame: SiderustFrame::ICRS,
        }
    }

    fn get_vega_handle() -> *mut SiderustStar {
        let cname = CString::new("VEGA").unwrap();
        let mut handle: *mut SiderustStar = ptr::null_mut();
        unsafe { crate::bodies::siderust_star_catalog(cname.as_ptr(), &mut handle) };
        handle
    }

    // ── window_from_c ─────────────────────────────────────────────────────

    #[test]
    fn window_invalid_period() {
        let bad = TempochPeriodMjd {
            start_mjd: 60001.0,
            end_mjd: 60000.0,
        };
        assert!(window_from_c(bad).is_err());
    }

    #[test]
    fn window_valid_period() {
        let w = TempochPeriodMjd {
            start_mjd: 60000.0,
            end_mjd: 60001.0,
        };
        assert!(window_from_c(w).is_ok());
    }

    // ── icrs_from_c ───────────────────────────────────────────────────────

    #[test]
    fn icrs_from_c_wrong_frame() {
        let dir = SiderustSphericalDir {
            polar_deg: 0.0,
            azimuth_deg: 0.0,
            frame: SiderustFrame::EquatorialMeanJ2000,
        };
        assert!(icrs_from_c(dir).is_err());
    }

    #[test]
    fn icrs_from_c_correct_frame() {
        assert!(icrs_from_c(icrs_vega()).is_ok());
    }

    // ── Sun altitude ──────────────────────────────────────────────────────

    #[test]
    fn sun_altitude_at_is_finite() {
        let mut out = 0.0f64;
        assert_eq!(
            siderust_sun_altitude_at(paris(), 60000.0, &mut out),
            SiderustStatus::Ok
        );
        assert!(out.is_finite());
    }

    #[test]
    fn sun_altitude_at_null_ptr() {
        assert_eq!(
            siderust_sun_altitude_at(paris(), 60000.0, ptr::null_mut()),
            SiderustStatus::NullPointer
        );
    }

    #[test]
    fn sun_above_threshold_ok() {
        let mut out: *mut TempochPeriodMjd = ptr::null_mut();
        let mut count = 0usize;
        let s = siderust_sun_above_threshold(
            paris(),
            one_day_window(),
            0.0,
            default_opts(),
            &mut out,
            &mut count,
        );
        assert_eq!(s, SiderustStatus::Ok);
        unsafe { siderust_periods_free(out, count) };
    }

    #[test]
    fn sun_below_threshold_ok() {
        let mut out: *mut TempochPeriodMjd = ptr::null_mut();
        let mut count = 0usize;
        let s = siderust_sun_below_threshold(
            paris(),
            one_day_window(),
            0.0,
            default_opts(),
            &mut out,
            &mut count,
        );
        assert_eq!(s, SiderustStatus::Ok);
        unsafe { siderust_periods_free(out, count) };
    }

    #[test]
    fn sun_crossings_ok() {
        let mut out: *mut SiderustCrossingEvent = ptr::null_mut();
        let mut count = 0usize;
        let s = siderust_sun_crossings(
            paris(),
            one_day_window(),
            0.0,
            default_opts(),
            &mut out,
            &mut count,
        );
        assert_eq!(s, SiderustStatus::Ok);
        unsafe { siderust_crossings_free(out, count) };
    }

    #[test]
    fn sun_culminations_ok() {
        let mut out: *mut SiderustCulminationEvent = ptr::null_mut();
        let mut count = 0usize;
        let s = siderust_sun_culminations(
            paris(),
            one_day_window(),
            default_opts(),
            &mut out,
            &mut count,
        );
        assert_eq!(s, SiderustStatus::Ok);
        unsafe { siderust_culminations_free(out, count) };
    }

    #[test]
    fn sun_altitude_periods_ok() {
        let query = SiderustAltitudeQuery {
            observer: paris(),
            start_mjd: 60000.0,
            end_mjd: 60001.0,
            min_altitude_deg: -90.0,
            max_altitude_deg: 90.0,
        };
        let mut out: *mut TempochPeriodMjd = ptr::null_mut();
        let mut count = 0usize;
        let s = siderust_sun_altitude_periods(query, &mut out, &mut count);
        assert_eq!(s, SiderustStatus::Ok);
        unsafe { siderust_periods_free(out, count) };
    }

    // ── Moon altitude ─────────────────────────────────────────────────────

    #[test]
    fn moon_altitude_at_is_finite() {
        let mut out = 0.0f64;
        assert_eq!(
            siderust_moon_altitude_at(paris(), 60000.0, &mut out),
            SiderustStatus::Ok
        );
        assert!(out.is_finite());
    }

    #[test]
    fn moon_altitude_at_null_ptr() {
        assert_eq!(
            siderust_moon_altitude_at(paris(), 60000.0, ptr::null_mut()),
            SiderustStatus::NullPointer
        );
    }

    #[test]
    fn moon_above_threshold_ok() {
        let mut out: *mut TempochPeriodMjd = ptr::null_mut();
        let mut count = 0usize;
        let s = siderust_moon_above_threshold(
            paris(),
            one_day_window(),
            0.0,
            default_opts(),
            &mut out,
            &mut count,
        );
        assert_eq!(s, SiderustStatus::Ok);
        unsafe { siderust_periods_free(out, count) };
    }

    #[test]
    fn moon_below_threshold_ok() {
        let mut out: *mut TempochPeriodMjd = ptr::null_mut();
        let mut count = 0usize;
        let s = siderust_moon_below_threshold(
            paris(),
            one_day_window(),
            0.0,
            default_opts(),
            &mut out,
            &mut count,
        );
        assert_eq!(s, SiderustStatus::Ok);
        unsafe { siderust_periods_free(out, count) };
    }

    #[test]
    fn moon_crossings_ok() {
        let mut out: *mut SiderustCrossingEvent = ptr::null_mut();
        let mut count = 0usize;
        let s = siderust_moon_crossings(
            paris(),
            one_day_window(),
            0.0,
            default_opts(),
            &mut out,
            &mut count,
        );
        assert_eq!(s, SiderustStatus::Ok);
        unsafe { siderust_crossings_free(out, count) };
    }

    #[test]
    fn moon_culminations_ok() {
        let mut out: *mut SiderustCulminationEvent = ptr::null_mut();
        let mut count = 0usize;
        let s = siderust_moon_culminations(
            paris(),
            one_day_window(),
            default_opts(),
            &mut out,
            &mut count,
        );
        assert_eq!(s, SiderustStatus::Ok);
        unsafe { siderust_culminations_free(out, count) };
    }

    #[test]
    fn moon_altitude_periods_ok() {
        let query = SiderustAltitudeQuery {
            observer: paris(),
            start_mjd: 60000.0,
            end_mjd: 60001.0,
            min_altitude_deg: -90.0,
            max_altitude_deg: 90.0,
        };
        let mut out: *mut TempochPeriodMjd = ptr::null_mut();
        let mut count = 0usize;
        let s = siderust_moon_altitude_periods(query, &mut out, &mut count);
        assert_eq!(s, SiderustStatus::Ok);
        unsafe { siderust_periods_free(out, count) };
    }

    // ── Star altitude ─────────────────────────────────────────────────────

    #[test]
    fn star_altitude_at_ok() {
        let h = get_vega_handle();
        let mut out = 0.0f64;
        assert_eq!(
            siderust_star_altitude_at(h, paris(), 60000.0, &mut out),
            SiderustStatus::Ok
        );
        assert!(out.is_finite());
        unsafe { crate::bodies::siderust_star_free(h) };
    }

    #[test]
    fn star_altitude_null_handle() {
        let mut out = 0.0f64;
        assert_eq!(
            siderust_star_altitude_at(ptr::null(), paris(), 60000.0, &mut out),
            SiderustStatus::NullPointer
        );
    }

    #[test]
    fn star_altitude_null_out() {
        let h = get_vega_handle();
        assert_eq!(
            siderust_star_altitude_at(h, paris(), 60000.0, ptr::null_mut()),
            SiderustStatus::NullPointer
        );
        unsafe { crate::bodies::siderust_star_free(h) };
    }

    #[test]
    fn star_above_threshold_ok() {
        let h = get_vega_handle();
        let mut out: *mut TempochPeriodMjd = ptr::null_mut();
        let mut count = 0usize;
        let s = siderust_star_above_threshold(
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
            siderust_periods_free(out, count);
            crate::bodies::siderust_star_free(h);
        }
    }

    #[test]
    fn star_above_threshold_null_handle() {
        let mut out: *mut TempochPeriodMjd = ptr::null_mut();
        let mut count = 0usize;
        assert_eq!(
            siderust_star_above_threshold(
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
    fn star_below_threshold_ok() {
        let h = get_vega_handle();
        let mut out: *mut TempochPeriodMjd = ptr::null_mut();
        let mut count = 0usize;
        let s = siderust_star_below_threshold(
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
            siderust_periods_free(out, count);
            crate::bodies::siderust_star_free(h);
        }
    }

    #[test]
    fn star_below_threshold_null_handle() {
        let mut out: *mut TempochPeriodMjd = ptr::null_mut();
        let mut count = 0usize;
        assert_eq!(
            siderust_star_below_threshold(
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
    fn star_crossings_ok() {
        let h = get_vega_handle();
        let mut out: *mut SiderustCrossingEvent = ptr::null_mut();
        let mut count = 0usize;
        let s = siderust_star_crossings(
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
            siderust_crossings_free(out, count);
            crate::bodies::siderust_star_free(h);
        }
    }

    #[test]
    fn star_crossings_null_handle() {
        let mut out: *mut SiderustCrossingEvent = ptr::null_mut();
        let mut count = 0usize;
        assert_eq!(
            siderust_star_crossings(
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
    fn star_culminations_ok() {
        let h = get_vega_handle();
        let mut out: *mut SiderustCulminationEvent = ptr::null_mut();
        let mut count = 0usize;
        let s = siderust_star_culminations(
            h,
            paris(),
            one_day_window(),
            default_opts(),
            &mut out,
            &mut count,
        );
        assert_eq!(s, SiderustStatus::Ok);
        unsafe {
            siderust_culminations_free(out, count);
            crate::bodies::siderust_star_free(h);
        }
    }

    #[test]
    fn star_culminations_null_handle() {
        let mut out: *mut SiderustCulminationEvent = ptr::null_mut();
        let mut count = 0usize;
        assert_eq!(
            siderust_star_culminations(
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

    // ── ICRS altitude ─────────────────────────────────────────────────────

    #[test]
    fn icrs_altitude_at_ok() {
        let mut out = 0.0f64;
        let s = siderust_icrs_altitude_at(icrs_vega(), paris(), 60000.0, &mut out);
        assert_eq!(s, SiderustStatus::Ok);
        assert!(out.is_finite());
    }

    #[test]
    fn icrs_altitude_null_out() {
        assert_eq!(
            siderust_icrs_altitude_at(icrs_vega(), paris(), 60000.0, ptr::null_mut()),
            SiderustStatus::NullPointer
        );
    }

    #[test]
    fn icrs_altitude_wrong_frame() {
        let dir = SiderustSphericalDir {
            polar_deg: 0.0,
            azimuth_deg: 0.0,
            frame: SiderustFrame::EquatorialMeanJ2000,
        };
        let mut out = 0.0f64;
        assert_eq!(
            siderust_icrs_altitude_at(dir, paris(), 60000.0, &mut out),
            SiderustStatus::InvalidFrame
        );
    }

    #[test]
    fn icrs_above_threshold_ok() {
        let mut out: *mut TempochPeriodMjd = ptr::null_mut();
        let mut count = 0usize;
        let s = siderust_icrs_above_threshold(
            icrs_vega(),
            paris(),
            one_day_window(),
            0.0,
            default_opts(),
            &mut out,
            &mut count,
        );
        assert_eq!(s, SiderustStatus::Ok);
        unsafe { siderust_periods_free(out, count) };
    }

    #[test]
    fn icrs_above_threshold_wrong_frame() {
        let dir = SiderustSphericalDir {
            polar_deg: 0.0,
            azimuth_deg: 0.0,
            frame: SiderustFrame::Galactic,
        };
        let mut out: *mut TempochPeriodMjd = ptr::null_mut();
        let mut count = 0usize;
        assert_eq!(
            siderust_icrs_above_threshold(
                dir,
                paris(),
                one_day_window(),
                0.0,
                default_opts(),
                &mut out,
                &mut count
            ),
            SiderustStatus::InvalidFrame
        );
    }

    #[test]
    fn icrs_below_threshold_ok() {
        let mut out: *mut TempochPeriodMjd = ptr::null_mut();
        let mut count = 0usize;
        let s = siderust_icrs_below_threshold(
            icrs_vega(),
            paris(),
            one_day_window(),
            0.0,
            default_opts(),
            &mut out,
            &mut count,
        );
        assert_eq!(s, SiderustStatus::Ok);
        unsafe { siderust_periods_free(out, count) };
    }

    #[test]
    fn icrs_below_threshold_wrong_frame() {
        let dir = SiderustSphericalDir {
            polar_deg: 0.0,
            azimuth_deg: 0.0,
            frame: SiderustFrame::Galactic,
        };
        let mut out: *mut TempochPeriodMjd = ptr::null_mut();
        let mut count = 0usize;
        assert_eq!(
            siderust_icrs_below_threshold(
                dir,
                paris(),
                one_day_window(),
                0.0,
                default_opts(),
                &mut out,
                &mut count
            ),
            SiderustStatus::InvalidFrame
        );
    }

    // ── Invalid window ────────────────────────────────────────────────────

    #[test]
    fn sun_above_threshold_invalid_window() {
        let bad = TempochPeriodMjd {
            start_mjd: 60001.0,
            end_mjd: 60000.0,
        };
        let mut out: *mut TempochPeriodMjd = ptr::null_mut();
        let mut count = 0usize;
        assert_eq!(
            siderust_sun_above_threshold(paris(), bad, 0.0, default_opts(), &mut out, &mut count),
            SiderustStatus::InvalidPeriod
        );
    }

    #[test]
    fn free_null_pointers_safe() {
        unsafe {
            siderust_periods_free(ptr::null_mut(), 0);
            siderust_crossings_free(ptr::null_mut(), 0);
            siderust_culminations_free(ptr::null_mut(), 0);
        }
    }
}
