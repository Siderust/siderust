// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Shared altitude helpers for the unified subject FFI.
//!
//! The public query entry points live in [`crate::subject`]. This module only
//! keeps the common window conversions and array-free helpers used by that
//! surface and by phase/event modules.

use crate::error::SiderustStatus;
use crate::ffi_utils::{free_boxed_slice, vec_to_c, FfiFrom};
use crate::types::*;
#[cfg(test)]
use qtty::*;
#[cfg(test)]
use siderust::coordinates::spherical;
use tempoch::{Interval, ModifiedJulianDate, Period, MJD};

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

#[cfg(test)]
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

/// Free an array of MJD periods.
///
/// # Safety
///
/// * `ptr` must have been returned by a siderust FFI function that allocates
///   `TempochPeriodMjd` arrays.
/// * `count` must be the element count returned alongside `ptr`.
/// * The pointer must not have been freed before, and must not be used after
///   this call.
#[no_mangle]
pub unsafe extern "C" fn siderust_periods_free(ptr: *mut TempochPeriodMjd, count: usize) {
    unsafe { free_boxed_slice(ptr, count) };
}

/// Free an array of crossing events.
///
/// # Safety
///
/// * `ptr` must have been returned by a siderust FFI function that allocates
///   `SiderustCrossingEvent` arrays.
/// * `count` must be the element count returned alongside `ptr`.
/// * The pointer must not have been freed before, and must not be used after
///   this call.
#[no_mangle]
pub unsafe extern "C" fn siderust_crossings_free(ptr: *mut SiderustCrossingEvent, count: usize) {
    unsafe { free_boxed_slice(ptr, count) };
}

/// Free an array of culmination events.
///
/// # Safety
///
/// * `ptr` must have been returned by a siderust FFI function that allocates
///   `SiderustCulminationEvent` arrays.
/// * `count` must be the element count returned alongside `ptr`.
/// * The pointer must not have been freed before, and must not be used after
///   this call.
#[no_mangle]
pub unsafe extern "C" fn siderust_culminations_free(
    ptr: *mut SiderustCulminationEvent,
    count: usize,
) {
    unsafe { free_boxed_slice(ptr, count) };
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_helpers::{icrs_vega, one_day_window};

    #[test]
    fn window_invalid_period() {
        let bad = TempochPeriodMjd {
            start_mjd: 60001.0,
            end_mjd: 60000.0,
        };
        assert_eq!(window_from_c(bad), Err(SiderustStatus::InvalidPeriod));
    }

    #[test]
    fn window_valid_period() {
        assert!(window_from_c(one_day_window()).is_ok());
    }

    #[test]
    fn icrs_from_c_wrong_frame() {
        let dir = SiderustSphericalDir {
            polar_deg: 0.0,
            azimuth_deg: 0.0,
            frame: SiderustFrame::EquatorialMeanJ2000,
        };
        assert!(matches!(
            icrs_from_c(dir),
            Err(SiderustStatus::InvalidFrame)
        ));
    }

    #[test]
    fn icrs_from_c_correct_frame() {
        assert!(icrs_from_c(icrs_vega()).is_ok());
    }
}
