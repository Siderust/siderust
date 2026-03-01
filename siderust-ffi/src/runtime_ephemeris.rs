// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! FFI bindings for the runtime-loaded ephemeris backend.
//!
//! All functions in this module work with an opaque [`SiderustRuntimeEphemeris`]
//! handle that wraps a [`RuntimeEphemeris`](siderust::calculus::ephemeris::RuntimeEphemeris).
//!
//! ## Lifecycle
//!
//! 1. Create a handle via `siderust_runtime_ephemeris_load_bsp` (from file)
//!    or `siderust_runtime_ephemeris_load_bytes` (from memory).
//! 2. Query positions via `siderust_runtime_ephemeris_sun_barycentric`, etc.
//! 3. Free the handle via `siderust_runtime_ephemeris_free`.
//!
//! With the `runtime-data` feature, `siderust_runtime_ephemeris_ensure` can
//! download a BSP file and load it in one step.

use std::ffi::CStr;
use std::os::raw::c_char;

use crate::error::SiderustStatus;
use crate::types::*;
use siderust::calculus::ephemeris::{DynEphemeris, RuntimeEphemeris};
use siderust::time::JulianDate;

// ═══════════════════════════════════════════════════════════════════════════
// Opaque handle
// ═══════════════════════════════════════════════════════════════════════════

/// Opaque handle to a `RuntimeEphemeris`.
///
/// Created via `siderust_runtime_ephemeris_load_bsp`, `_load_bytes`, or
/// `_ensure`. Must be freed with `siderust_runtime_ephemeris_free`.
pub struct SiderustRuntimeEphemeris {
    pub(crate) inner: RuntimeEphemeris,
}

// ═══════════════════════════════════════════════════════════════════════════
// Constructors
// ═══════════════════════════════════════════════════════════════════════════

/// Load a runtime ephemeris from a BSP file on disk.
///
/// `path` must be a valid, NUL-terminated UTF-8 path to a JPL DE4xx BSP file.
/// On success, `*out` receives a newly allocated handle that must be freed
/// with `siderust_runtime_ephemeris_free`.
#[no_mangle]
pub extern "C" fn siderust_runtime_ephemeris_load_bsp(
    path: *const c_char,
    out: *mut *mut SiderustRuntimeEphemeris,
) -> SiderustStatus {
    ffi_guard! {{
        if path.is_null() || out.is_null() {
            return SiderustStatus::NullPointer;
        }
        let path_str = match unsafe { CStr::from_ptr(path) }.to_str() {
            Ok(s) => s,
            Err(_) => return SiderustStatus::InvalidArgument,
        };
        match RuntimeEphemeris::from_bsp(path_str) {
            Ok(eph) => {
                let handle = Box::new(SiderustRuntimeEphemeris { inner: eph });
                unsafe { *out = Box::into_raw(handle) };
                SiderustStatus::Ok
            }
            Err(_) => SiderustStatus::DataError,
        }
    }}
}

/// Load a runtime ephemeris from raw BSP bytes already in memory.
///
/// `data` must point to `len` bytes of a valid JPL DE4xx BSP file.
/// On success, `*out` receives a newly allocated handle.
#[no_mangle]
pub extern "C" fn siderust_runtime_ephemeris_load_bytes(
    data: *const u8,
    len: usize,
    out: *mut *mut SiderustRuntimeEphemeris,
) -> SiderustStatus {
    ffi_guard! {{
        if data.is_null() || out.is_null() {
            return SiderustStatus::NullPointer;
        }
        let slice = unsafe { std::slice::from_raw_parts(data, len) };
        match RuntimeEphemeris::from_bytes(slice) {
            Ok(eph) => {
                let handle = Box::new(SiderustRuntimeEphemeris { inner: eph });
                unsafe { *out = Box::into_raw(handle) };
                SiderustStatus::Ok
            }
            Err(_) => SiderustStatus::DataError,
        }
    }}
}

/// Download (if necessary) and load a runtime ephemeris for the given dataset.
///
/// `dataset_id`:
///   - 0 = DE440 (~120 MB)
///   - 1 = DE441 (~1.65 GB)
///
/// On first call this will download the BSP from JPL's servers into the
/// local cache (`~/.siderust/data/` or `$SIDERUST_DATA_DIR`). Subsequent
/// calls re-use the cached file.
///
/// On success, `*out` receives a newly allocated handle.
///
/// **Requires feature `runtime-data`.**
#[cfg(feature = "runtime-data")]
#[no_mangle]
pub extern "C" fn siderust_runtime_ephemeris_ensure(
    dataset_id: u32,
    out: *mut *mut SiderustRuntimeEphemeris,
) -> SiderustStatus {
    ffi_guard! {{
        if out.is_null() {
            return SiderustStatus::NullPointer;
        }
        let id = match dataset_id {
            0 => siderust::data::DatasetId::De440,
            1 => siderust::data::DatasetId::De441,
            _ => return SiderustStatus::InvalidArgument,
        };
        let dm = match siderust::data::DataManager::new() {
            Ok(dm) => dm,
            Err(_) => return SiderustStatus::DataError,
        };
        match RuntimeEphemeris::from_data_manager(&dm, id) {
            Ok(eph) => {
                let handle = Box::new(SiderustRuntimeEphemeris { inner: eph });
                unsafe { *out = Box::into_raw(handle) };
                SiderustStatus::Ok
            }
            Err(_) => SiderustStatus::DataError,
        }
    }}
}

// ═══════════════════════════════════════════════════════════════════════════
// Destructor
// ═══════════════════════════════════════════════════════════════════════════

/// Free a `SiderustRuntimeEphemeris` handle.
///
/// # Safety
///
/// The handle must have been allocated by one of the `siderust_runtime_ephemeris_*`
/// constructors, and must not be used after this call.
#[no_mangle]
pub unsafe extern "C" fn siderust_runtime_ephemeris_free(
    handle: *mut SiderustRuntimeEphemeris,
) {
    if !handle.is_null() {
        drop(Box::from_raw(handle));
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// Query functions
// ═══════════════════════════════════════════════════════════════════════════

/// Sun's barycentric position (EclipticMeanJ2000, AU) via the runtime ephemeris.
#[no_mangle]
pub extern "C" fn siderust_runtime_ephemeris_sun_barycentric(
    handle: *const SiderustRuntimeEphemeris,
    jd: f64,
    out: *mut SiderustCartesianPos,
) -> SiderustStatus {
    ffi_guard! {{
        if handle.is_null() || out.is_null() {
            return SiderustStatus::NullPointer;
        }
        let eph = unsafe { &*handle };
        let t = JulianDate::new(jd);
        let pos = eph.inner.sun_barycentric(t);
        unsafe {
            *out = SiderustCartesianPos {
                x: pos.x().value(),
                y: pos.y().value(),
                z: pos.z().value(),
                frame: SiderustFrame::EclipticMeanJ2000,
                center: SiderustCenter::Barycentric,
            };
        }
        SiderustStatus::Ok
    }}
}

/// Earth's barycentric position (EclipticMeanJ2000, AU) via the runtime ephemeris.
#[no_mangle]
pub extern "C" fn siderust_runtime_ephemeris_earth_barycentric(
    handle: *const SiderustRuntimeEphemeris,
    jd: f64,
    out: *mut SiderustCartesianPos,
) -> SiderustStatus {
    ffi_guard! {{
        if handle.is_null() || out.is_null() {
            return SiderustStatus::NullPointer;
        }
        let eph = unsafe { &*handle };
        let t = JulianDate::new(jd);
        let pos = eph.inner.earth_barycentric(t);
        unsafe {
            *out = SiderustCartesianPos {
                x: pos.x().value(),
                y: pos.y().value(),
                z: pos.z().value(),
                frame: SiderustFrame::EclipticMeanJ2000,
                center: SiderustCenter::Barycentric,
            };
        }
        SiderustStatus::Ok
    }}
}

/// Earth's heliocentric position (EclipticMeanJ2000, AU) via the runtime ephemeris.
#[no_mangle]
pub extern "C" fn siderust_runtime_ephemeris_earth_heliocentric(
    handle: *const SiderustRuntimeEphemeris,
    jd: f64,
    out: *mut SiderustCartesianPos,
) -> SiderustStatus {
    ffi_guard! {{
        if handle.is_null() || out.is_null() {
            return SiderustStatus::NullPointer;
        }
        let eph = unsafe { &*handle };
        let t = JulianDate::new(jd);
        let pos = eph.inner.earth_heliocentric(t);
        unsafe {
            *out = SiderustCartesianPos {
                x: pos.x().value(),
                y: pos.y().value(),
                z: pos.z().value(),
                frame: SiderustFrame::EclipticMeanJ2000,
                center: SiderustCenter::Heliocentric,
            };
        }
        SiderustStatus::Ok
    }}
}

/// Earth's barycentric velocity (EclipticMeanJ2000, AU/day) via the runtime ephemeris.
#[no_mangle]
pub extern "C" fn siderust_runtime_ephemeris_earth_barycentric_velocity(
    handle: *const SiderustRuntimeEphemeris,
    jd: f64,
    out: *mut SiderustCartesianVel,
) -> SiderustStatus {
    ffi_guard! {{
        if handle.is_null() || out.is_null() {
            return SiderustStatus::NullPointer;
        }
        let eph = unsafe { &*handle };
        let t = JulianDate::new(jd);
        let vel = eph.inner.earth_barycentric_velocity(t);
        unsafe {
            *out = SiderustCartesianVel {
                vx: vel.x().value(),
                vy: vel.y().value(),
                vz: vel.z().value(),
                frame: SiderustFrame::EclipticMeanJ2000,
            };
        }
        SiderustStatus::Ok
    }}
}

/// Moon's geocentric position (EclipticMeanJ2000, km) via the runtime ephemeris.
#[no_mangle]
pub extern "C" fn siderust_runtime_ephemeris_moon_geocentric(
    handle: *const SiderustRuntimeEphemeris,
    jd: f64,
    out: *mut SiderustCartesianPos,
) -> SiderustStatus {
    ffi_guard! {{
        if handle.is_null() || out.is_null() {
            return SiderustStatus::NullPointer;
        }
        let eph = unsafe { &*handle };
        let t = JulianDate::new(jd);
        let pos = eph.inner.moon_geocentric(t);
        unsafe {
            *out = SiderustCartesianPos {
                x: pos.x().value(),
                y: pos.y().value(),
                z: pos.z().value(),
                frame: SiderustFrame::EclipticMeanJ2000,
                center: SiderustCenter::Geocentric,
            };
        }
        SiderustStatus::Ok
    }}
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::ptr;

    #[test]
    fn load_bsp_null_path() {
        let mut out: *mut SiderustRuntimeEphemeris = ptr::null_mut();
        assert_eq!(
            siderust_runtime_ephemeris_load_bsp(ptr::null(), &mut out),
            SiderustStatus::NullPointer
        );
    }

    #[test]
    fn load_bsp_null_out() {
        let path = b"dummy.bsp\0";
        assert_eq!(
            siderust_runtime_ephemeris_load_bsp(
                path.as_ptr() as *const c_char,
                ptr::null_mut()
            ),
            SiderustStatus::NullPointer
        );
    }

    #[test]
    fn load_bsp_nonexistent_file() {
        let path = b"/nonexistent/path/de440.bsp\0";
        let mut out: *mut SiderustRuntimeEphemeris = ptr::null_mut();
        assert_eq!(
            siderust_runtime_ephemeris_load_bsp(
                path.as_ptr() as *const c_char,
                &mut out
            ),
            SiderustStatus::DataError
        );
        assert!(out.is_null());
    }

    #[test]
    fn load_bytes_null_data() {
        let mut out: *mut SiderustRuntimeEphemeris = ptr::null_mut();
        assert_eq!(
            siderust_runtime_ephemeris_load_bytes(ptr::null(), 0, &mut out),
            SiderustStatus::NullPointer
        );
    }

    #[test]
    fn load_bytes_invalid_data() {
        let bad_data = [0u8; 64];
        let mut out: *mut SiderustRuntimeEphemeris = ptr::null_mut();
        assert_eq!(
            siderust_runtime_ephemeris_load_bytes(
                bad_data.as_ptr(),
                bad_data.len(),
                &mut out
            ),
            SiderustStatus::DataError
        );
    }

    #[test]
    fn query_null_handle() {
        let mut pos = SiderustCartesianPos {
            x: 0.0,
            y: 0.0,
            z: 0.0,
            frame: SiderustFrame::ICRS,
            center: SiderustCenter::Barycentric,
        };
        assert_eq!(
            siderust_runtime_ephemeris_sun_barycentric(ptr::null(), 2451545.0, &mut pos),
            SiderustStatus::NullPointer
        );
        assert_eq!(
            siderust_runtime_ephemeris_earth_barycentric(ptr::null(), 2451545.0, &mut pos),
            SiderustStatus::NullPointer
        );
        assert_eq!(
            siderust_runtime_ephemeris_earth_heliocentric(ptr::null(), 2451545.0, &mut pos),
            SiderustStatus::NullPointer
        );
        assert_eq!(
            siderust_runtime_ephemeris_moon_geocentric(ptr::null(), 2451545.0, &mut pos),
            SiderustStatus::NullPointer
        );
        let mut vel = SiderustCartesianVel {
            vx: 0.0,
            vy: 0.0,
            vz: 0.0,
            frame: SiderustFrame::ICRS,
        };
        assert_eq!(
            siderust_runtime_ephemeris_earth_barycentric_velocity(ptr::null(), 2451545.0, &mut vel),
            SiderustStatus::NullPointer
        );
    }

    #[test]
    fn free_null_is_safe() {
        unsafe { siderust_runtime_ephemeris_free(ptr::null_mut()) };
    }
}
