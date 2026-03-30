// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Shared azimuth helpers for the unified subject FFI.
//!
//! The public azimuth query entry points live in [`crate::subject`]. This
//! module only contains C-array conversion and free helpers.

use crate::error::SiderustStatus;
use crate::ffi_utils::{free_boxed_slice, vec_to_c, FfiFrom};
use crate::types::*;
use siderust::calculus::azimuth::{AzimuthCrossingEvent, AzimuthExtremum};

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

/// Free an array of azimuth crossing events.
///
/// # Safety
///
/// * `ptr` must have been returned by a siderust FFI function that allocates
///   `SiderustAzimuthCrossingEvent` arrays.
/// * `count` must be the element count returned alongside `ptr`.
/// * The pointer must not have been freed before, and must not be used after
///   this call.
#[no_mangle]
pub unsafe extern "C" fn siderust_azimuth_crossings_free(
    ptr: *mut SiderustAzimuthCrossingEvent,
    count: usize,
) {
    unsafe { free_boxed_slice(ptr, count) };
}

/// Free an array of azimuth extrema.
///
/// # Safety
///
/// * `ptr` must have been returned by a siderust FFI function that allocates
///   `SiderustAzimuthExtremum` arrays.
/// * `count` must be the element count returned alongside `ptr`.
/// * The pointer must not have been freed before, and must not be used after
///   this call.
#[no_mangle]
pub unsafe extern "C" fn siderust_azimuth_extrema_free(
    ptr: *mut SiderustAzimuthExtremum,
    count: usize,
) {
    unsafe { free_boxed_slice(ptr, count) };
}
