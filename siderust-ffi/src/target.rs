// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Target FFI — opaque handle for an ICRS celestial target.
//!
//! A `SiderustTarget` wraps an ICRS direction (RA/Dec) and exposes altitude
//! and azimuth queries using the same algorithms as the `icrs_*` family.

use crate::altitude::{crossings_to_c, culminations_to_c, periods_to_c, window_from_c};
use crate::azimuth::vec_az_crossings_to_c;
use crate::error::SiderustStatus;
use crate::types::*;
use qtty::*;
use siderust::calculus::altitude::{
    above_threshold, crossings, culminations, AltitudePeriodsProvider,
};
use siderust::calculus::azimuth::azimuth_crossings;
use siderust::coordinates::spherical;
use siderust::time::ModifiedJulianDate;
use siderust::AzimuthProvider;

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
    if handle.is_null() || out.is_null() {
        return SiderustStatus::NullPointer;
    }
    unsafe { *out = (*handle).ra_deg };
    SiderustStatus::Ok
}

/// Write the declination (degrees) to `*out`.
#[no_mangle]
pub extern "C" fn siderust_target_dec_deg(
    handle: *const SiderustTarget,
    out: *mut f64,
) -> SiderustStatus {
    if handle.is_null() || out.is_null() {
        return SiderustStatus::NullPointer;
    }
    unsafe { *out = (*handle).dec_deg };
    SiderustStatus::Ok
}

/// Write the reference epoch (Julian Date) to `*out`.
#[no_mangle]
pub extern "C" fn siderust_target_epoch_jd(
    handle: *const SiderustTarget,
    out: *mut f64,
) -> SiderustStatus {
    if handle.is_null() || out.is_null() {
        return SiderustStatus::NullPointer;
    }
    unsafe { *out = (*handle).epoch_jd };
    SiderustStatus::Ok
}

// ═══════════════════════════════════════════════════════════════════════════
// Altitude
// ═══════════════════════════════════════════════════════════════════════════

/// Altitude of the target (radians).
#[no_mangle]
pub extern "C" fn siderust_target_altitude_at(
    handle: *const SiderustTarget,
    observer: SiderustGeodetict,
    mjd: f64,
    out_rad: *mut f64,
) -> SiderustStatus {
    if handle.is_null() || out_rad.is_null() {
        return SiderustStatus::NullPointer;
    }
    let dir = unsafe { &(*handle).dir };
    unsafe {
        *out_rad = dir
            .altitude_at(&observer.to_rust(), ModifiedJulianDate::new(mjd))
            .value();
    }
    SiderustStatus::Ok
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
    if handle.is_null() {
        return SiderustStatus::NullPointer;
    }
    let window = match window_from_c(window) {
        Ok(w) => w,
        Err(e) => return e,
    };
    let dir = unsafe { &(*handle).dir };
    periods_to_c(
        above_threshold(
            dir,
            &observer.to_rust(),
            window,
            Degrees::new(threshold_deg),
            opts.to_rust(),
        ),
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
    if handle.is_null() {
        return SiderustStatus::NullPointer;
    }
    let window = match window_from_c(window) {
        Ok(w) => w,
        Err(e) => return e,
    };
    let dir = unsafe { &(*handle).dir };
    crossings_to_c(
        crossings(
            dir,
            &observer.to_rust(),
            window,
            Degrees::new(threshold_deg),
            opts.to_rust(),
        ),
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
    if handle.is_null() {
        return SiderustStatus::NullPointer;
    }
    let window = match window_from_c(window) {
        Ok(w) => w,
        Err(e) => return e,
    };
    let dir = unsafe { &(*handle).dir };
    culminations_to_c(
        culminations(dir, &observer.to_rust(), window, opts.to_rust()),
        out,
        count,
    )
}

// ═══════════════════════════════════════════════════════════════════════════
// Azimuth
// ═══════════════════════════════════════════════════════════════════════════

/// Azimuth of the target (degrees, North-clockwise).
#[no_mangle]
pub extern "C" fn siderust_target_azimuth_at(
    handle: *const SiderustTarget,
    observer: SiderustGeodetict,
    mjd: f64,
    out_deg: *mut f64,
) -> SiderustStatus {
    if handle.is_null() || out_deg.is_null() {
        return SiderustStatus::NullPointer;
    }
    let dir = unsafe { &(*handle).dir };
    unsafe {
        *out_deg = dir
            .azimuth_at(&observer.to_rust(), ModifiedJulianDate::new(mjd))
            .value();
    }
    SiderustStatus::Ok
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
    if handle.is_null() {
        return SiderustStatus::NullPointer;
    }
    let window = match window_from_c(window) {
        Ok(w) => w,
        Err(e) => return e,
    };
    let dir = unsafe { &(*handle).dir };
    vec_az_crossings_to_c(
        azimuth_crossings(
            dir,
            &observer.to_rust(),
            window,
            Degrees::new(bearing_deg),
            opts.to_rust(),
        ),
        out,
        count,
    )
}
