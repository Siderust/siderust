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
        Degrees::new(dir.lon_deg),
        Degrees::new(dir.lat_deg),
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
