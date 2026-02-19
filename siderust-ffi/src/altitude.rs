// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! FFI bindings for the altitude / observation API.
//!
//! Provides per-body-type functions for altitude_at, above_threshold,
//! below_threshold, crossings, culminations, and altitude_periods.

use crate::bodies::SiderustStar;
use crate::error::SiderustStatus;
use crate::types::*;
use qtty::*;
use siderust::bodies::{Moon, Sun};
use siderust::coordinates::spherical;
use siderust::AltitudePeriodsProvider;
use tempoch::{Interval, ModifiedJulianDate, Period, MJD};

// Re-use the TempochPeriodMjd from our types module (re-exported from tempoch_ffi).

// ═══════════════════════════════════════════════════════════════════════════
// Helper: convert Vec<Period<MJD>> to allocated C array
// ═══════════════════════════════════════════════════════════════════════════

fn periods_to_c(
    periods: Vec<Period<MJD>>,
    out: *mut *mut TempochPeriodMjd,
    count: *mut usize,
) -> SiderustStatus {
    if out.is_null() || count.is_null() {
        return SiderustStatus::NullPointer;
    }
    let c_periods: Vec<TempochPeriodMjd> = periods
        .iter()
        .map(|p| TempochPeriodMjd {
            start_mjd: p.start.value(),
            end_mjd: p.end.value(),
        })
        .collect();
    let len = c_periods.len();
    let boxed = c_periods.into_boxed_slice();
    let ptr = Box::into_raw(boxed) as *mut TempochPeriodMjd;
    unsafe {
        *out = ptr;
        *count = len;
    }
    SiderustStatus::Ok
}

fn crossings_to_c(
    events: Vec<siderust::CrossingEvent>,
    out: *mut *mut SiderustCrossingEvent,
    count: *mut usize,
) -> SiderustStatus {
    if out.is_null() || count.is_null() {
        return SiderustStatus::NullPointer;
    }
    let c_events: Vec<SiderustCrossingEvent> = events
        .iter()
        .map(SiderustCrossingEvent::from_rust)
        .collect();
    let len = c_events.len();
    let boxed = c_events.into_boxed_slice();
    let ptr = Box::into_raw(boxed) as *mut SiderustCrossingEvent;
    unsafe {
        *out = ptr;
        *count = len;
    }
    SiderustStatus::Ok
}

fn culminations_to_c(
    events: Vec<siderust::CulminationEvent>,
    out: *mut *mut SiderustCulminationEvent,
    count: *mut usize,
) -> SiderustStatus {
    if out.is_null() || count.is_null() {
        return SiderustStatus::NullPointer;
    }
    let c_events: Vec<SiderustCulminationEvent> = events
        .iter()
        .map(SiderustCulminationEvent::from_rust)
        .collect();
    let len = c_events.len();
    let boxed = c_events.into_boxed_slice();
    let ptr = Box::into_raw(boxed) as *mut SiderustCulminationEvent;
    unsafe {
        *out = ptr;
        *count = len;
    }
    SiderustStatus::Ok
}

// ═══════════════════════════════════════════════════════════════════════════
// Free functions for dynamic arrays
// ═══════════════════════════════════════════════════════════════════════════

/// Free an array of crossing events.
#[no_mangle]
pub unsafe extern "C" fn siderust_crossings_free(
    ptr: *mut SiderustCrossingEvent,
    count: usize,
) {
    if !ptr.is_null() && count > 0 {
        let _ = Box::from_raw(std::slice::from_raw_parts_mut(ptr, count));
    }
}

/// Free an array of culmination events.
#[no_mangle]
pub unsafe extern "C" fn siderust_culminations_free(
    ptr: *mut SiderustCulminationEvent,
    count: usize,
) {
    if !ptr.is_null() && count > 0 {
        let _ = Box::from_raw(std::slice::from_raw_parts_mut(ptr, count));
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// Sun altitude functions
// ═══════════════════════════════════════════════════════════════════════════

/// Compute the altitude of the Sun at a given instant (radians).
#[no_mangle]
pub extern "C" fn siderust_sun_altitude_at(
    observer: SiderustGeodetict,
    mjd: f64,
    out_rad: *mut f64,
) -> SiderustStatus {
    if out_rad.is_null() { return SiderustStatus::NullPointer; }
    let obs = observer.to_rust();
    let t = ModifiedJulianDate::new(mjd);
    let alt = Sun.altitude_at(&obs, t);
    unsafe { *out_rad = alt.value() };
    SiderustStatus::Ok
}

/// Find periods when the Sun is above a threshold altitude.
#[no_mangle]
pub extern "C" fn siderust_sun_above_threshold(
    observer: SiderustGeodetict,
    start_mjd: f64,
    end_mjd: f64,
    threshold_deg: f64,
    opts: SiderustSearchOpts,
    out: *mut *mut TempochPeriodMjd,
    count: *mut usize,
) -> SiderustStatus {
    let obs = observer.to_rust();
    let window = Interval::new(ModifiedJulianDate::new(start_mjd), ModifiedJulianDate::new(end_mjd));
    let periods = siderust::above_threshold(&Sun, &obs, window, Degrees::new(threshold_deg), opts.to_rust());
    periods_to_c(periods, out, count)
}

/// Find periods when the Sun is below a threshold altitude.
#[no_mangle]
pub extern "C" fn siderust_sun_below_threshold(
    observer: SiderustGeodetict,
    start_mjd: f64,
    end_mjd: f64,
    threshold_deg: f64,
    opts: SiderustSearchOpts,
    out: *mut *mut TempochPeriodMjd,
    count: *mut usize,
) -> SiderustStatus {
    let obs = observer.to_rust();
    let window = Interval::new(ModifiedJulianDate::new(start_mjd), ModifiedJulianDate::new(end_mjd));
    let periods = siderust::below_threshold(&Sun, &obs, window, Degrees::new(threshold_deg), opts.to_rust());
    periods_to_c(periods, out, count)
}

/// Find threshold-crossing events for the Sun.
#[no_mangle]
pub extern "C" fn siderust_sun_crossings(
    observer: SiderustGeodetict,
    start_mjd: f64,
    end_mjd: f64,
    threshold_deg: f64,
    opts: SiderustSearchOpts,
    out: *mut *mut SiderustCrossingEvent,
    count: *mut usize,
) -> SiderustStatus {
    let obs = observer.to_rust();
    let window = Interval::new(ModifiedJulianDate::new(start_mjd), ModifiedJulianDate::new(end_mjd));
    let events = siderust::crossings(&Sun, &obs, window, Degrees::new(threshold_deg), opts.to_rust());
    crossings_to_c(events, out, count)
}

/// Find culmination (local extrema) events for the Sun.
#[no_mangle]
pub extern "C" fn siderust_sun_culminations(
    observer: SiderustGeodetict,
    start_mjd: f64,
    end_mjd: f64,
    opts: SiderustSearchOpts,
    out: *mut *mut SiderustCulminationEvent,
    count: *mut usize,
) -> SiderustStatus {
    let obs = observer.to_rust();
    let window = Interval::new(ModifiedJulianDate::new(start_mjd), ModifiedJulianDate::new(end_mjd));
    let events = siderust::culminations(&Sun, &obs, window, opts.to_rust());
    culminations_to_c(events, out, count)
}

/// Find periods when the Sun's altitude is within [min, max].
#[no_mangle]
pub extern "C" fn siderust_sun_altitude_periods(
    query: SiderustAltitudeQuery,
    out: *mut *mut TempochPeriodMjd,
    count: *mut usize,
) -> SiderustStatus {
    let q = query.to_rust();
    let periods = Sun.altitude_periods(&q);
    periods_to_c(periods, out, count)
}

// ═══════════════════════════════════════════════════════════════════════════
// Moon altitude functions
// ═══════════════════════════════════════════════════════════════════════════

/// Compute the altitude of the Moon at a given instant (radians).
#[no_mangle]
pub extern "C" fn siderust_moon_altitude_at(
    observer: SiderustGeodetict,
    mjd: f64,
    out_rad: *mut f64,
) -> SiderustStatus {
    if out_rad.is_null() { return SiderustStatus::NullPointer; }
    let obs = observer.to_rust();
    let t = ModifiedJulianDate::new(mjd);
    let alt = Moon.altitude_at(&obs, t);
    unsafe { *out_rad = alt.value() };
    SiderustStatus::Ok
}

/// Find periods when the Moon is above a threshold altitude.
#[no_mangle]
pub extern "C" fn siderust_moon_above_threshold(
    observer: SiderustGeodetict,
    start_mjd: f64,
    end_mjd: f64,
    threshold_deg: f64,
    opts: SiderustSearchOpts,
    out: *mut *mut TempochPeriodMjd,
    count: *mut usize,
) -> SiderustStatus {
    let obs = observer.to_rust();
    let window = Interval::new(ModifiedJulianDate::new(start_mjd), ModifiedJulianDate::new(end_mjd));
    let periods = siderust::above_threshold(&Moon, &obs, window, Degrees::new(threshold_deg), opts.to_rust());
    periods_to_c(periods, out, count)
}

/// Find periods when the Moon is below a threshold altitude.
#[no_mangle]
pub extern "C" fn siderust_moon_below_threshold(
    observer: SiderustGeodetict,
    start_mjd: f64,
    end_mjd: f64,
    threshold_deg: f64,
    opts: SiderustSearchOpts,
    out: *mut *mut TempochPeriodMjd,
    count: *mut usize,
) -> SiderustStatus {
    let obs = observer.to_rust();
    let window = Interval::new(ModifiedJulianDate::new(start_mjd), ModifiedJulianDate::new(end_mjd));
    let periods = siderust::below_threshold(&Moon, &obs, window, Degrees::new(threshold_deg), opts.to_rust());
    periods_to_c(periods, out, count)
}

/// Find threshold-crossing events for the Moon.
#[no_mangle]
pub extern "C" fn siderust_moon_crossings(
    observer: SiderustGeodetict,
    start_mjd: f64,
    end_mjd: f64,
    threshold_deg: f64,
    opts: SiderustSearchOpts,
    out: *mut *mut SiderustCrossingEvent,
    count: *mut usize,
) -> SiderustStatus {
    let obs = observer.to_rust();
    let window = Interval::new(ModifiedJulianDate::new(start_mjd), ModifiedJulianDate::new(end_mjd));
    let events = siderust::crossings(&Moon, &obs, window, Degrees::new(threshold_deg), opts.to_rust());
    crossings_to_c(events, out, count)
}

/// Find culmination (local extrema) events for the Moon.
#[no_mangle]
pub extern "C" fn siderust_moon_culminations(
    observer: SiderustGeodetict,
    start_mjd: f64,
    end_mjd: f64,
    opts: SiderustSearchOpts,
    out: *mut *mut SiderustCulminationEvent,
    count: *mut usize,
) -> SiderustStatus {
    let obs = observer.to_rust();
    let window = Interval::new(ModifiedJulianDate::new(start_mjd), ModifiedJulianDate::new(end_mjd));
    let events = siderust::culminations(&Moon, &obs, window, opts.to_rust());
    culminations_to_c(events, out, count)
}

/// Find periods when the Moon's altitude is within [min, max].
#[no_mangle]
pub extern "C" fn siderust_moon_altitude_periods(
    query: SiderustAltitudeQuery,
    out: *mut *mut TempochPeriodMjd,
    count: *mut usize,
) -> SiderustStatus {
    let q = query.to_rust();
    let periods = Moon.altitude_periods(&q);
    periods_to_c(periods, out, count)
}

// ═══════════════════════════════════════════════════════════════════════════
// Star altitude functions (take an opaque handle)
// ═══════════════════════════════════════════════════════════════════════════

/// Compute the altitude of a star at a given instant (radians).
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
    let obs = observer.to_rust();
    let t = ModifiedJulianDate::new(mjd);
    let alt = star.altitude_at(&obs, t);
    unsafe { *out_rad = alt.value() };
    SiderustStatus::Ok
}

/// Find periods when a star is above a threshold altitude.
#[no_mangle]
pub extern "C" fn siderust_star_above_threshold(
    handle: *const SiderustStar,
    observer: SiderustGeodetict,
    start_mjd: f64,
    end_mjd: f64,
    threshold_deg: f64,
    opts: SiderustSearchOpts,
    out: *mut *mut TempochPeriodMjd,
    count: *mut usize,
) -> SiderustStatus {
    if handle.is_null() {
        return SiderustStatus::NullPointer;
    }
    let star = unsafe { &(*handle).inner };
    let obs = observer.to_rust();
    let window = Interval::new(
        ModifiedJulianDate::new(start_mjd),
        ModifiedJulianDate::new(end_mjd),
    );
    let periods = siderust::above_threshold(
        star, &obs, window, Degrees::new(threshold_deg), opts.to_rust(),
    );
    periods_to_c(periods, out, count)
}

/// Find periods when a star is below a threshold altitude.
#[no_mangle]
pub extern "C" fn siderust_star_below_threshold(
    handle: *const SiderustStar,
    observer: SiderustGeodetict,
    start_mjd: f64,
    end_mjd: f64,
    threshold_deg: f64,
    opts: SiderustSearchOpts,
    out: *mut *mut TempochPeriodMjd,
    count: *mut usize,
) -> SiderustStatus {
    if handle.is_null() {
        return SiderustStatus::NullPointer;
    }
    let star = unsafe { &(*handle).inner };
    let obs = observer.to_rust();
    let window = Interval::new(
        ModifiedJulianDate::new(start_mjd),
        ModifiedJulianDate::new(end_mjd),
    );
    let periods = siderust::below_threshold(
        star, &obs, window, Degrees::new(threshold_deg), opts.to_rust(),
    );
    periods_to_c(periods, out, count)
}

/// Find crossing events for a star.
#[no_mangle]
pub extern "C" fn siderust_star_crossings(
    handle: *const SiderustStar,
    observer: SiderustGeodetict,
    start_mjd: f64,
    end_mjd: f64,
    threshold_deg: f64,
    opts: SiderustSearchOpts,
    out: *mut *mut SiderustCrossingEvent,
    count: *mut usize,
) -> SiderustStatus {
    if handle.is_null() {
        return SiderustStatus::NullPointer;
    }
    let star = unsafe { &(*handle).inner };
    let obs = observer.to_rust();
    let window = Interval::new(
        ModifiedJulianDate::new(start_mjd),
        ModifiedJulianDate::new(end_mjd),
    );
    let events = siderust::crossings(
        star, &obs, window, Degrees::new(threshold_deg), opts.to_rust(),
    );
    crossings_to_c(events, out, count)
}

/// Find culmination events for a star.
#[no_mangle]
pub extern "C" fn siderust_star_culminations(
    handle: *const SiderustStar,
    observer: SiderustGeodetict,
    start_mjd: f64,
    end_mjd: f64,
    opts: SiderustSearchOpts,
    out: *mut *mut SiderustCulminationEvent,
    count: *mut usize,
) -> SiderustStatus {
    if handle.is_null() {
        return SiderustStatus::NullPointer;
    }
    let star = unsafe { &(*handle).inner };
    let obs = observer.to_rust();
    let window = Interval::new(
        ModifiedJulianDate::new(start_mjd),
        ModifiedJulianDate::new(end_mjd),
    );
    let events = siderust::culminations(star, &obs, window, opts.to_rust());
    culminations_to_c(events, out, count)
}

// ═══════════════════════════════════════════════════════════════════════════
// ICRS direction altitude functions (lightweight — no body object)
// ═══════════════════════════════════════════════════════════════════════════

/// Compute the altitude of a fixed ICRS direction at a given instant (radians).
#[no_mangle]
pub extern "C" fn siderust_icrs_dir_altitude_at(
    ra_deg: f64,
    dec_deg: f64,
    observer: SiderustGeodetict,
    mjd: f64,
    out_rad: *mut f64,
) -> SiderustStatus {
    if out_rad.is_null() {
        return SiderustStatus::NullPointer;
    }
    let dir = spherical::direction::ICRS::new(Degrees::new(ra_deg), Degrees::new(dec_deg));
    let obs = observer.to_rust();
    let t = ModifiedJulianDate::new(mjd);
    let alt = dir.altitude_at(&obs, t);
    unsafe { *out_rad = alt.value() };
    SiderustStatus::Ok
}

/// Find periods when an ICRS direction is above a threshold altitude.
#[no_mangle]
pub extern "C" fn siderust_icrs_dir_above_threshold(
    ra_deg: f64,
    dec_deg: f64,
    observer: SiderustGeodetict,
    start_mjd: f64,
    end_mjd: f64,
    threshold_deg: f64,
    opts: SiderustSearchOpts,
    out: *mut *mut TempochPeriodMjd,
    count: *mut usize,
) -> SiderustStatus {
    let dir = spherical::direction::ICRS::new(Degrees::new(ra_deg), Degrees::new(dec_deg));
    let obs = observer.to_rust();
    let window = Interval::new(
        ModifiedJulianDate::new(start_mjd),
        ModifiedJulianDate::new(end_mjd),
    );
    let periods = siderust::above_threshold(
        &dir, &obs, window, Degrees::new(threshold_deg), opts.to_rust(),
    );
    periods_to_c(periods, out, count)
}

/// Find periods when an ICRS direction is below a threshold altitude.
#[no_mangle]
pub extern "C" fn siderust_icrs_dir_below_threshold(
    ra_deg: f64,
    dec_deg: f64,
    observer: SiderustGeodetict,
    start_mjd: f64,
    end_mjd: f64,
    threshold_deg: f64,
    opts: SiderustSearchOpts,
    out: *mut *mut TempochPeriodMjd,
    count: *mut usize,
) -> SiderustStatus {
    let dir = spherical::direction::ICRS::new(Degrees::new(ra_deg), Degrees::new(dec_deg));
    let obs = observer.to_rust();
    let window = Interval::new(
        ModifiedJulianDate::new(start_mjd),
        ModifiedJulianDate::new(end_mjd),
    );
    let periods = siderust::below_threshold(
        &dir, &obs, window, Degrees::new(threshold_deg), opts.to_rust(),
    );
    periods_to_c(periods, out, count)
}
