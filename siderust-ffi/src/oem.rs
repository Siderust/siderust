// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! C FFI for CCSDS Orbit Ephemeris Message (OEM) parsing.
//!
//! Provides a simple flat-array interface: parse an OEM text document and
//! receive a C array of [`SiderustOemState`] values (epoch + position + velocity).
//! All segments from the OEM file are concatenated into a single flat array.
//!
//! ## Lifecycle
//!
//! 1. Call [`siderust_oem_parse_str`] with a NUL-terminated OEM text.
//! 2. Iterate over the returned `out_states` array (`out_count` elements).
//! 3. Free the array with [`siderust_oem_states_free`] when done.
//!
//! ## Status codes
//!
//! | Value | Meaning |
//! |-------|---------|
//! | 0 | Success |
//! | 1 | A required pointer was null |
//! | 9 | OEM parse error |

use std::ffi::CStr;
use std::os::raw::{c_char, c_ulong};

use siderust::formats::ccsds::oem::read_oem;

use crate::error::SiderustStatus;

// ─── Public types ─────────────────────────────────────────────────────────────

/// A single spacecraft state vector as parsed from a CCSDS OEM file.
///
/// All numeric fields use the conventional OEM units: position in km,
/// velocity in km/s, epoch as a Julian Date.
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct SiderustOemState {
    /// Epoch as a Julian Date (any time system declared in the OEM metadata).
    pub epoch_jd: f64,
    /// Position components [x, y, z] in km.
    pub pos_km: [f64; 3],
    /// Velocity components [vx, vy, vz] in km/s.
    pub vel_kms: [f64; 3],
}

// ─── Entry points ─────────────────────────────────────────────────────────────

/// Parse a CCSDS OEM (KVN) document from a NUL-terminated C string.
///
/// All state vectors from all segments are flattened into a single heap-
/// allocated C array returned through `out_states`.  The caller must free it
/// with [`siderust_oem_states_free`] when no longer needed.
///
/// # Parameters
///
/// * `text`       — NUL-terminated OEM document (UTF-8 or ASCII).
/// * `out_states` — receives a `*mut SiderustOemState` pointing to the array.
/// * `out_count`  — receives the number of elements in the array.
///
/// # Returns
///
/// [`SiderustStatus::Ok`] on success (including zero states).
/// [`SiderustStatus::NullPointer`] if any required pointer is null.
/// [`SiderustStatus::InvalidArgument`] if the OEM document cannot be parsed.
#[no_mangle]
pub extern "C" fn siderust_oem_parse_str(
    text: *const c_char,
    out_states: *mut *mut SiderustOemState,
    out_count: *mut c_ulong,
) -> SiderustStatus {
    crate::ffi_guard!({
        if text.is_null() {
            return SiderustStatus::NullPointer;
        }
        crate::check_out!(out_states, out_count);

        // SAFETY: caller guarantees NUL-terminated text.
        let s = match unsafe { CStr::from_ptr(text) }.to_str() {
            Ok(s) => s,
            Err(_) => return SiderustStatus::InvalidArgument,
        };

        let oem = match read_oem(s.as_bytes()) {
            Ok(f) => f,
            Err(_) => return SiderustStatus::InvalidArgument,
        };

        let mut states: Vec<SiderustOemState> = Vec::new();
        for seg in oem.segments {
            for st in seg.states {
                states.push(SiderustOemState {
                    epoch_jd: st.epoch_jd,
                    pos_km: st.position_km,
                    vel_kms: st.velocity_km_s,
                });
            }
        }

        let count = states.len() as c_ulong;
        // Transfer ownership to C caller via raw pointer.
        let ptr = if states.is_empty() {
            // Box<[]> has a well-defined non-null dangling pointer; use
            // a proper allocation so free is symmetric.
            std::ptr::null_mut()
        } else {
            let mut boxed: Box<[SiderustOemState]> = states.into_boxed_slice();
            let ptr = boxed.as_mut_ptr();
            std::mem::forget(boxed);
            ptr
        };

        // SAFETY: out_states / out_count checked non-null above.
        unsafe {
            *out_states = ptr;
            *out_count = count;
        }
        SiderustStatus::Ok
    })
}

/// Free an array of [`SiderustOemState`] returned by [`siderust_oem_parse_str`].
///
/// # Safety
///
/// `states` must have been returned by [`siderust_oem_parse_str`] and
/// `count` must match the value written into `out_count`.
/// Passing null or count 0 is a no-op.
#[no_mangle]
pub extern "C" fn siderust_oem_states_free(states: *mut SiderustOemState, count: c_ulong) {
    if states.is_null() || count == 0 {
        return;
    }
    // SAFETY: states/count match the values produced by siderust_oem_parse_str.
    unsafe {
        let slice = std::slice::from_raw_parts_mut(states, count as usize);
        drop(Box::from_raw(slice as *mut [SiderustOemState]));
    }
}
