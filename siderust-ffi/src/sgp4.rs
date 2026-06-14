// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! C FFI for TLE parsing and the SGP4 propagator.
//!
//! ## Lifecycle
//!
//! ### TLE handle
//!
//! 1. Parse a two-line element set: `siderust_tle_parse(line1, line2, out)`.
//! 2. Inspect the NORAD catalog number: `siderust_tle_norad_id(tle, out)`.
//! 3. Pass the handle to [`siderust_sgp4_new`].
//! 4. Free when done: `siderust_tle_free(tle)`.
//!
//! ### SGP4 propagator handle
//!
//! 1. Create: `siderust_sgp4_new(tle, gravity_model, out)`.
//! 2. Query the TLE epoch: `siderust_sgp4_epoch_jd_utc(sgp4, out_jd)`.
//! 3. Propagate to a Julian date: `siderust_sgp4_propagate_at(sgp4, jd, out_pos, out_vel)`.
//! 4. Free when done: `siderust_sgp4_free(sgp4)`.
//!
//! ## Status codes
//!
//! | Value | Meaning |
//! |-------|---------|
//! | 0 | Success |
//! | 1 | A required pointer was null |
//! | 9 | TLE parse error or propagation failure |
//! | 10 | Internal panic (should never happen) |

use std::ffi::CStr;
use std::os::raw::c_char;

use qtty::time::Day;
use qtty::Quantity;
use siderust::astro::sgp4::{GravityModel, Sgp4Propagator};
use siderust::formats::tle::parse_tle;
use tempoch::{JulianDate, UTC};

use crate::error::SiderustStatus;

// ─── Gravity model selector ───────────────────────────────────────────────────

/// Gravity model selector for the SGP4 propagator.
///
/// | Value | Meaning |
/// |-------|---------|
/// | 0 | WGS-72 with AFSPC sidereal-time convention (default) |
/// | 1 | WGS-72 with IAU sidereal-time formula |
/// | 2 | WGS-84 with IAU sidereal-time formula |
#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SiderustGravityModel {
    /// WGS-72 / AFSPC — matches the TLE mean-element conventions.
    Wgs72 = 0,
    /// WGS-72 / IAU sidereal time.
    Wgs72Iau = 1,
    /// WGS-84 / IAU sidereal time — modern default for new pipelines.
    Wgs84 = 2,
}

impl From<SiderustGravityModel> for GravityModel {
    fn from(m: SiderustGravityModel) -> Self {
        match m {
            SiderustGravityModel::Wgs72 => GravityModel::Wgs72,
            SiderustGravityModel::Wgs72Iau => GravityModel::Wgs72Iau,
            SiderustGravityModel::Wgs84 => GravityModel::Wgs84,
        }
    }
}

// ─── Opaque TLE handle ────────────────────────────────────────────────────────

/// Opaque handle to a parsed Two-Line Element set.
///
/// Owned by the caller; must be freed with [`siderust_tle_free`].
pub struct SiderustTle {
    inner: siderust::formats::tle::TLE,
}

/// Parse a two-line element set from two NUL-terminated C strings.
///
/// # Parameters
///
/// * `line1` — TLE line 1 (NUL-terminated).
/// * `line2` — TLE line 2 (NUL-terminated).
/// * `out`   — receives a `*mut SiderustTle` on success; must **not** be null.
///
/// # Returns
///
/// [`SiderustStatus::Ok`] on success.
/// [`SiderustStatus::NullPointer`] if any pointer is null.
/// [`SiderustStatus::InvalidArgument`] if parsing fails.
#[no_mangle]
pub extern "C" fn siderust_tle_parse(
    line1: *const c_char,
    line2: *const c_char,
    out: *mut *mut SiderustTle,
) -> SiderustStatus {
    crate::ffi_guard!({
        crate::check_out!(out);
        if line1.is_null() || line2.is_null() {
            return SiderustStatus::NullPointer;
        }
        // SAFETY: caller guarantees NUL-terminated strings.
        let l1 = unsafe { CStr::from_ptr(line1) }
            .to_str()
            .unwrap_or("")
            .trim();
        let l2 = unsafe { CStr::from_ptr(line2) }
            .to_str()
            .unwrap_or("")
            .trim();

        match parse_tle(l1, l2) {
            Ok(tle) => {
                let boxed = Box::new(SiderustTle { inner: tle });
                // SAFETY: out checked non-null above.
                unsafe { *out = Box::into_raw(boxed) };
                SiderustStatus::Ok
            }
            Err(_) => SiderustStatus::InvalidArgument,
        }
    })
}

/// Return the NORAD catalog number of a TLE handle.
///
/// # Returns
///
/// [`SiderustStatus::Ok`] on success.
/// [`SiderustStatus::NullPointer`] if either pointer is null.
#[no_mangle]
pub extern "C" fn siderust_tle_norad_id(
    tle: *const SiderustTle,
    out_id: *mut u32,
) -> SiderustStatus {
    crate::ffi_guard!({
        crate::check_out!(out_id);
        if tle.is_null() {
            return SiderustStatus::NullPointer;
        }
        let tle = unsafe { &*tle };
        // SAFETY: out_id checked non-null above.
        unsafe { *out_id = tle.inner.norad_id.0 };
        SiderustStatus::Ok
    })
}

/// Free a TLE handle obtained from [`siderust_tle_parse`].
///
/// Passing null is a no-op.
#[no_mangle]
pub extern "C" fn siderust_tle_free(tle: *mut SiderustTle) {
    if !tle.is_null() {
        // SAFETY: tle was produced by Box::into_raw in siderust_tle_parse.
        drop(unsafe { Box::from_raw(tle) });
    }
}

// ─── Opaque SGP4 propagator handle ───────────────────────────────────────────

/// Opaque handle to an initialised SGP4 propagator.
///
/// Owned by the caller; must be freed with [`siderust_sgp4_free`].
pub struct SiderustSgp4 {
    inner: Sgp4Propagator,
}

/// Create an SGP4 propagator from a TLE handle.
///
/// # Parameters
///
/// * `tle`           — TLE handle (not consumed; caller still owns it).
/// * `gravity_model` — 0 = WGS-72 (default), 1 = WGS-72/IAU, 2 = WGS-84.
/// * `out`           — receives a `*mut SiderustSgp4` on success.
///
/// # Returns
///
/// [`SiderustStatus::Ok`] on success.
/// [`SiderustStatus::NullPointer`] if any pointer is null.
/// [`SiderustStatus::InvalidArgument`] if propagator initialisation fails.
#[no_mangle]
pub extern "C" fn siderust_sgp4_new(
    tle: *const SiderustTle,
    gravity_model: i32,
    out: *mut *mut SiderustSgp4,
) -> SiderustStatus {
    crate::ffi_guard!({
        crate::check_out!(out);
        if tle.is_null() {
            return SiderustStatus::NullPointer;
        }
        let tle = unsafe { &*tle };
        let model = match gravity_model {
            1 => GravityModel::Wgs72Iau,
            2 => GravityModel::Wgs84,
            _ => GravityModel::Wgs72,
        };
        match Sgp4Propagator::from_tle_with_model(&tle.inner, model) {
            Ok(prop) => {
                let boxed = Box::new(SiderustSgp4 { inner: prop });
                unsafe { *out = Box::into_raw(boxed) };
                SiderustStatus::Ok
            }
            Err(_) => SiderustStatus::InvalidArgument,
        }
    })
}

/// Return the gravity model used by an SGP4 propagator handle.
///
/// * `out_model` receives 0 (WGS-72), 1 (WGS-72/IAU), or 2 (WGS-84).
#[no_mangle]
pub extern "C" fn siderust_sgp4_gravity_model(
    sgp4: *const SiderustSgp4,
    out_model: *mut i32,
) -> SiderustStatus {
    crate::ffi_guard!({
        crate::check_out!(out_model);
        if sgp4.is_null() {
            return SiderustStatus::NullPointer;
        }
        let sgp4 = unsafe { &*sgp4 };
        let code = match sgp4.inner.gravity_model() {
            GravityModel::Wgs72 => 0,
            GravityModel::Wgs72Iau => 1,
            GravityModel::Wgs84 => 2,
        };
        unsafe { *out_model = code };
        SiderustStatus::Ok
    })
}

/// Return the UTC Julian date of the TLE epoch stored in an SGP4 handle.
///
/// `out_jd` receives the Julian date in days.
#[no_mangle]
pub extern "C" fn siderust_sgp4_epoch_jd_utc(
    sgp4: *const SiderustSgp4,
    out_jd: *mut f64,
) -> SiderustStatus {
    crate::ffi_guard!({
        crate::check_out!(out_jd);
        if sgp4.is_null() {
            return SiderustStatus::NullPointer;
        }
        let sgp4 = unsafe { &*sgp4 };
        unsafe { *out_jd = sgp4.inner.epoch_jd_utc().raw().value() };
        SiderustStatus::Ok
    })
}

/// Propagate an SGP4 handle to a UTC Julian date.
///
/// # Parameters
///
/// * `sgp4`        — propagator handle (must not be null).
/// * `jd_utc`      — target epoch as a UTC Julian date (days).
/// * `out_pos_km`  — receives position [x, y, z] in km (TEME frame, 3 elements).
/// * `out_vel_kms` — receives velocity [vx, vy, vz] in km/s (TEME frame, 3 elements).
///
/// # Returns
///
/// [`SiderustStatus::Ok`] on success.
/// [`SiderustStatus::NullPointer`] if any required pointer is null.
/// [`SiderustStatus::InvalidArgument`] if propagation fails (e.g. singular orbit).
#[no_mangle]
pub extern "C" fn siderust_sgp4_propagate_at(
    sgp4: *const SiderustSgp4,
    jd_utc: f64,
    out_pos_km: *mut f64,
    out_vel_kms: *mut f64,
) -> SiderustStatus {
    crate::ffi_guard!({
        crate::check_out!(out_pos_km);
        crate::check_out!(out_vel_kms);
        if sgp4.is_null() {
            return SiderustStatus::NullPointer;
        }
        let sgp4 = unsafe { &*sgp4 };

        let target = match JulianDate::<UTC>::try_new(Quantity::<Day>::new(jd_utc)) {
            Ok(jd) => jd,
            Err(_) => return SiderustStatus::InvalidArgument,
        };

        match sgp4.inner.propagate_at(target) {
            Ok(state) => {
                let p = state.position().as_array();
                let v = state.velocity().as_array();
                // SAFETY: out_pos_km / out_vel_kms point to at least 3 f64 each.
                unsafe {
                    *out_pos_km.add(0) = p[0].value();
                    *out_pos_km.add(1) = p[1].value();
                    *out_pos_km.add(2) = p[2].value();
                    *out_vel_kms.add(0) = v[0].value();
                    *out_vel_kms.add(1) = v[1].value();
                    *out_vel_kms.add(2) = v[2].value();
                }
                SiderustStatus::Ok
            }
            Err(_) => SiderustStatus::InvalidArgument,
        }
    })
}

/// Free an SGP4 handle obtained from [`siderust_sgp4_new`].
///
/// Passing null is a no-op.
#[no_mangle]
pub extern "C" fn siderust_sgp4_free(sgp4: *mut SiderustSgp4) {
    if !sgp4.is_null() {
        // SAFETY: sgp4 was produced by Box::into_raw in siderust_sgp4_new.
        drop(unsafe { Box::from_raw(sgp4) });
    }
}
