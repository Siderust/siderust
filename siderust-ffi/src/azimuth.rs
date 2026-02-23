// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Azimuth FFI — bearing crossings, extrema, and range periods.

use crate::altitude::{periods_to_c, window_from_c};
use crate::bodies::SiderustStar;
use crate::error::SiderustStatus;
use crate::ffi_utils::{free_boxed_slice, vec_to_c, FfiFrom};
use crate::types::*;
use qtty::*;
use siderust::bodies::{Moon, Sun};
use siderust::calculus::azimuth::{
    azimuth_crossings, azimuth_extrema, in_azimuth_range, AzimuthCrossingEvent, AzimuthExtremum,
};
use siderust::coordinates::spherical;
use siderust::time::ModifiedJulianDate;
use siderust::AzimuthProvider;

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
// Sun azimuth
// ═══════════════════════════════════════════════════════════════════════════

/// Azimuth of the Sun at an instant (degrees, North-clockwise).
#[no_mangle]
pub extern "C" fn siderust_sun_azimuth_at(
    observer: SiderustGeodetict,
    mjd: f64,
    out_deg: *mut f64,
) -> SiderustStatus {
    if out_deg.is_null() {
        return SiderustStatus::NullPointer;
    }
    unsafe {
        *out_deg = Sun
            .azimuth_at(&observer.to_rust(), ModifiedJulianDate::new(mjd))
            .value();
    }
    SiderustStatus::Ok
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
    let window = match window_from_c(window) {
        Ok(w) => w,
        Err(e) => return e,
    };
    vec_az_crossings_to_c(
        azimuth_crossings(
            &Sun,
            &observer.to_rust(),
            window,
            Degrees::new(bearing_deg),
            opts.to_rust(),
        ),
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
    let window = match window_from_c(window) {
        Ok(w) => w,
        Err(e) => return e,
    };
    vec_az_extrema_to_c(
        azimuth_extrema(&Sun, &observer.to_rust(), window, opts.to_rust()),
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
    let window = match window_from_c(window) {
        Ok(w) => w,
        Err(e) => return e,
    };
    periods_to_c(
        in_azimuth_range(
            &Sun,
            &observer.to_rust(),
            window,
            Degrees::new(min_deg),
            Degrees::new(max_deg),
            opts.to_rust(),
        ),
        out,
        count,
    )
}

// ═══════════════════════════════════════════════════════════════════════════
// Moon azimuth
// ═══════════════════════════════════════════════════════════════════════════

/// Azimuth of the Moon at an instant (degrees, North-clockwise).
#[no_mangle]
pub extern "C" fn siderust_moon_azimuth_at(
    observer: SiderustGeodetict,
    mjd: f64,
    out_deg: *mut f64,
) -> SiderustStatus {
    if out_deg.is_null() {
        return SiderustStatus::NullPointer;
    }
    unsafe {
        *out_deg = Moon
            .azimuth_at(&observer.to_rust(), ModifiedJulianDate::new(mjd))
            .value();
    }
    SiderustStatus::Ok
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
    let window = match window_from_c(window) {
        Ok(w) => w,
        Err(e) => return e,
    };
    vec_az_crossings_to_c(
        azimuth_crossings(
            &Moon,
            &observer.to_rust(),
            window,
            Degrees::new(bearing_deg),
            opts.to_rust(),
        ),
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
    let window = match window_from_c(window) {
        Ok(w) => w,
        Err(e) => return e,
    };
    vec_az_extrema_to_c(
        azimuth_extrema(&Moon, &observer.to_rust(), window, opts.to_rust()),
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
    let window = match window_from_c(window) {
        Ok(w) => w,
        Err(e) => return e,
    };
    periods_to_c(
        in_azimuth_range(
            &Moon,
            &observer.to_rust(),
            window,
            Degrees::new(min_deg),
            Degrees::new(max_deg),
            opts.to_rust(),
        ),
        out,
        count,
    )
}

// ═══════════════════════════════════════════════════════════════════════════
// Star azimuth
// ═══════════════════════════════════════════════════════════════════════════

/// Azimuth of a star at an instant (degrees, North-clockwise).
#[no_mangle]
pub extern "C" fn siderust_star_azimuth_at(
    handle: *const SiderustStar,
    observer: SiderustGeodetict,
    mjd: f64,
    out_deg: *mut f64,
) -> SiderustStatus {
    if handle.is_null() || out_deg.is_null() {
        return SiderustStatus::NullPointer;
    }
    let star = unsafe { &(*handle).inner };
    unsafe {
        *out_deg = star
            .azimuth_at(&observer.to_rust(), ModifiedJulianDate::new(mjd))
            .value();
    }
    SiderustStatus::Ok
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
    if handle.is_null() {
        return SiderustStatus::NullPointer;
    }
    let window = match window_from_c(window) {
        Ok(w) => w,
        Err(e) => return e,
    };
    let star = unsafe { &(*handle).inner };
    vec_az_crossings_to_c(
        azimuth_crossings(
            star,
            &observer.to_rust(),
            window,
            Degrees::new(bearing_deg),
            opts.to_rust(),
        ),
        out,
        count,
    )
}

// ═══════════════════════════════════════════════════════════════════════════
// ICRS direction azimuth
// ═══════════════════════════════════════════════════════════════════════════

/// Azimuth of an ICRS direction at an instant.
///
/// `dir.frame` must equal `SIDERUST_FRAME_T_ICRS`; otherwise
/// `SIDERUST_STATUS_T_INVALID_FRAME` is returned.
#[no_mangle]
pub extern "C" fn siderust_icrs_azimuth_at(
    dir: SiderustSphericalDir,
    observer: SiderustGeodetict,
    mjd: f64,
    out_deg: *mut f64,
) -> SiderustStatus {
    if out_deg.is_null() {
        return SiderustStatus::NullPointer;
    }
    if dir.frame != SiderustFrame::ICRS {
        return SiderustStatus::InvalidFrame;
    }
    let icrs =
        spherical::direction::ICRS::new(Degrees::new(dir.azimuth_deg), Degrees::new(dir.polar_deg));
    unsafe {
        *out_deg = icrs
            .azimuth_at(&observer.to_rust(), ModifiedJulianDate::new(mjd))
            .value();
    }
    SiderustStatus::Ok
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
    let window = match window_from_c(window) {
        Ok(w) => w,
        Err(e) => return e,
    };
    let icrs =
        spherical::direction::ICRS::new(Degrees::new(dir.azimuth_deg), Degrees::new(dir.polar_deg));
    vec_az_crossings_to_c(
        azimuth_crossings(
            &icrs,
            &observer.to_rust(),
            window,
            Degrees::new(bearing_deg),
            opts.to_rust(),
        ),
        out,
        count,
    )
}
