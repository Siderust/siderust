// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! FFI bindings for **generic body dispatch**: altitude, azimuth, crossings,
//! and culminations for any solar-system body identified by [`SiderustBody`].
//!
//! This module provides a single set of `siderust_body_*` functions that accept
//! a [`SiderustBody`] discriminant and dispatch to the appropriate concrete
//! implementation (Sun, Moon, or VSOP87 planet).

use crate::altitude::{crossings_to_c, culminations_to_c, periods_to_c, window_from_c};
use crate::azimuth::{vec_az_crossings_to_c, vec_az_extrema_to_c};
use crate::error::SiderustStatus;
use crate::types::*;
use qtty::*;
use siderust::bodies::solar_system;
use siderust::bodies::{Moon, Sun};
use siderust::calculus::azimuth::{azimuth_crossings, azimuth_extrema, in_azimuth_range};
use siderust::AltitudePeriodsProvider;
use siderust::AzimuthProvider;
use tempoch::ModifiedJulianDate;

// ═══════════════════════════════════════════════════════════════════════════
// Dispatch helper
// ═══════════════════════════════════════════════════════════════════════════

/// Dispatches an operation over the concrete body type corresponding to
/// `body`.  The closure receives a `&dyn AltitudePeriodsProvider` and can
/// call any trait method on it.
///
/// This avoids duplicating a 10-arm match in every FFI function.
macro_rules! dispatch_body {
    ($body:expr, |$provider:ident| $action:expr) => {
        match $body {
            SiderustBody::Sun => {
                let $provider = Sun;
                $action
            }
            SiderustBody::Moon => {
                let $provider = Moon;
                $action
            }
            SiderustBody::Mercury => {
                let $provider = solar_system::Mercury;
                $action
            }
            SiderustBody::Venus => {
                let $provider = solar_system::Venus;
                $action
            }
            SiderustBody::Mars => {
                let $provider = solar_system::Mars;
                $action
            }
            SiderustBody::Jupiter => {
                let $provider = solar_system::Jupiter;
                $action
            }
            SiderustBody::Saturn => {
                let $provider = solar_system::Saturn;
                $action
            }
            SiderustBody::Uranus => {
                let $provider = solar_system::Uranus;
                $action
            }
            SiderustBody::Neptune => {
                let $provider = solar_system::Neptune;
                $action
            }
        }
    };
}

// ═══════════════════════════════════════════════════════════════════════════
// Altitude
// ═══════════════════════════════════════════════════════════════════════════

/// Altitude of a solar-system body at an instant (radians).
#[no_mangle]
pub extern "C" fn siderust_body_altitude_at(
    body: SiderustBody,
    observer: SiderustGeodetict,
    mjd: f64,
    out_rad: *mut f64,
) -> SiderustStatus {
    if out_rad.is_null() {
        return SiderustStatus::NullPointer;
    }
    dispatch_body!(body, |b| {
        unsafe {
            *out_rad = b
                .altitude_at(&observer.to_rust(), ModifiedJulianDate::new(mjd))
                .value();
        }
        SiderustStatus::Ok
    })
}

/// Periods when a solar-system body is above a threshold altitude.
#[no_mangle]
pub extern "C" fn siderust_body_above_threshold(
    body: SiderustBody,
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
    dispatch_body!(body, |b| {
        periods_to_c(
            siderust::above_threshold(
                &b,
                &observer.to_rust(),
                window,
                Degrees::new(threshold_deg),
                opts.to_rust(),
            ),
            out,
            count,
        )
    })
}

/// Periods when a solar-system body is below a threshold altitude.
#[no_mangle]
pub extern "C" fn siderust_body_below_threshold(
    body: SiderustBody,
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
    dispatch_body!(body, |b| {
        periods_to_c(
            siderust::below_threshold(
                &b,
                &observer.to_rust(),
                window,
                Degrees::new(threshold_deg),
                opts.to_rust(),
            ),
            out,
            count,
        )
    })
}

/// Threshold-crossing events for a solar-system body.
#[no_mangle]
pub extern "C" fn siderust_body_crossings(
    body: SiderustBody,
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
    dispatch_body!(body, |b| {
        crossings_to_c(
            siderust::crossings(
                &b,
                &observer.to_rust(),
                window,
                Degrees::new(threshold_deg),
                opts.to_rust(),
            ),
            out,
            count,
        )
    })
}

/// Culmination (local extrema) events for a solar-system body.
#[no_mangle]
pub extern "C" fn siderust_body_culminations(
    body: SiderustBody,
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
    dispatch_body!(body, |b| {
        culminations_to_c(
            siderust::culminations(&b, &observer.to_rust(), window, opts.to_rust()),
            out,
            count,
        )
    })
}

/// Periods when a solar-system body's altitude is within [min, max].
#[no_mangle]
pub extern "C" fn siderust_body_altitude_periods(
    body: SiderustBody,
    query: SiderustAltitudeQuery,
    out: *mut *mut TempochPeriodMjd,
    count: *mut usize,
) -> SiderustStatus {
    let q = query.to_rust();
    dispatch_body!(body, |b| {
        periods_to_c(b.altitude_periods(&q), out, count)
    })
}

// ═══════════════════════════════════════════════════════════════════════════
// Azimuth
// ═══════════════════════════════════════════════════════════════════════════

/// Azimuth of a solar-system body at an instant (radians).
#[no_mangle]
pub extern "C" fn siderust_body_azimuth_at(
    body: SiderustBody,
    observer: SiderustGeodetict,
    mjd: f64,
    out_rad: *mut f64,
) -> SiderustStatus {
    if out_rad.is_null() {
        return SiderustStatus::NullPointer;
    }
    dispatch_body!(body, |b| {
        unsafe {
            *out_rad = b
                .azimuth_at(&observer.to_rust(), ModifiedJulianDate::new(mjd))
                .value();
        }
        SiderustStatus::Ok
    })
}

/// Azimuth bearing-crossing events for a solar-system body.
#[no_mangle]
pub extern "C" fn siderust_body_azimuth_crossings(
    body: SiderustBody,
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
    dispatch_body!(body, |b| {
        vec_az_crossings_to_c(
            azimuth_crossings(
                &b,
                &observer.to_rust(),
                window,
                Degrees::new(bearing_deg),
                opts.to_rust(),
            ),
            out,
            count,
        )
    })
}

/// Azimuth extrema (northernmost and southernmost bearing) for a body.
#[no_mangle]
pub extern "C" fn siderust_body_azimuth_extrema(
    body: SiderustBody,
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
    dispatch_body!(body, |b| {
        vec_az_extrema_to_c(
            azimuth_extrema(&b, &observer.to_rust(), window, opts.to_rust()),
            out,
            count,
        )
    })
}

/// Periods when a body's azimuth is within [min_deg, max_deg].
#[no_mangle]
pub extern "C" fn siderust_body_in_azimuth_range(
    body: SiderustBody,
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
    dispatch_body!(body, |b| {
        periods_to_c(
            in_azimuth_range(
                &b,
                &observer.to_rust(),
                window,
                Degrees::new(min_deg),
                Degrees::new(max_deg),
                opts.to_rust(),
            ),
            out,
            count,
        )
    })
}

// ═══════════════════════════════════════════════════════════════════════════
// Tests
// ═══════════════════════════════════════════════════════════════════════════

#[cfg(test)]
mod tests {
    use super::*;
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

    #[test]
    fn body_altitude_at_sun_is_finite() {
        let mut out = 0.0f64;
        let st = siderust_body_altitude_at(SiderustBody::Sun, paris(), 60000.5, &mut out);
        assert_eq!(st, SiderustStatus::Ok);
        assert!(out.is_finite());
    }

    #[test]
    fn body_altitude_at_mars_is_finite() {
        let mut out = 0.0f64;
        let st = siderust_body_altitude_at(SiderustBody::Mars, paris(), 60000.5, &mut out);
        assert_eq!(st, SiderustStatus::Ok);
        assert!(out.is_finite());
    }

    #[test]
    fn body_altitude_at_null_ptr() {
        let st = siderust_body_altitude_at(SiderustBody::Sun, paris(), 60000.5, ptr::null_mut());
        assert_eq!(st, SiderustStatus::NullPointer);
    }

    #[test]
    fn body_above_threshold_sun() {
        let mut out: *mut TempochPeriodMjd = ptr::null_mut();
        let mut count: usize = 0;
        let st = siderust_body_above_threshold(
            SiderustBody::Sun,
            paris(),
            one_day_window(),
            0.0,
            default_opts(),
            &mut out,
            &mut count,
        );
        assert_eq!(st, SiderustStatus::Ok);
        assert!(count > 0, "Sun should be above horizon");
        unsafe { crate::altitude::siderust_periods_free(out, count) };
    }

    #[test]
    fn body_above_threshold_jupiter() {
        let mut out: *mut TempochPeriodMjd = ptr::null_mut();
        let mut count: usize = 0;
        let window = TempochPeriodMjd {
            start_mjd: 60000.0,
            end_mjd: 60007.0, // one week
        };
        let st = siderust_body_above_threshold(
            SiderustBody::Jupiter,
            paris(),
            window,
            0.0,
            default_opts(),
            &mut out,
            &mut count,
        );
        assert_eq!(st, SiderustStatus::Ok);
        assert!(count > 0, "Jupiter should be above horizon at some point");
        unsafe { crate::altitude::siderust_periods_free(out, count) };
    }

    #[test]
    fn body_crossings_mars() {
        let mut out: *mut SiderustCrossingEvent = ptr::null_mut();
        let mut count: usize = 0;
        let window = TempochPeriodMjd {
            start_mjd: 60000.0,
            end_mjd: 60003.0,
        };
        let st = siderust_body_crossings(
            SiderustBody::Mars,
            paris(),
            window,
            0.0,
            default_opts(),
            &mut out,
            &mut count,
        );
        assert_eq!(st, SiderustStatus::Ok);
        // Mars may or may not have crossings in 3 days, but the function should succeed
        if count > 0 {
            unsafe { crate::altitude::siderust_crossings_free(out, count) };
        }
    }

    #[test]
    fn body_culminations_venus() {
        let mut out: *mut SiderustCulminationEvent = ptr::null_mut();
        let mut count: usize = 0;
        let window = TempochPeriodMjd {
            start_mjd: 60000.0,
            end_mjd: 60003.0,
        };
        let st = siderust_body_culminations(
            SiderustBody::Venus,
            paris(),
            window,
            default_opts(),
            &mut out,
            &mut count,
        );
        assert_eq!(st, SiderustStatus::Ok);
        if count > 0 {
            unsafe { crate::altitude::siderust_culminations_free(out, count) };
        }
    }

    #[test]
    fn all_bodies_return_finite_altitude() {
        let bodies = [
            SiderustBody::Sun,
            SiderustBody::Moon,
            SiderustBody::Mercury,
            SiderustBody::Venus,
            SiderustBody::Mars,
            SiderustBody::Jupiter,
            SiderustBody::Saturn,
            SiderustBody::Uranus,
            SiderustBody::Neptune,
        ];
        for body in &bodies {
            let mut out = 0.0f64;
            let st = siderust_body_altitude_at(*body, paris(), 60000.5, &mut out);
            assert_eq!(st, SiderustStatus::Ok, "Failed for {:?}", body);
            assert!(
                out.is_finite(),
                "Non-finite altitude for {:?}: {}",
                body,
                out
            );
            assert!(
                out.abs() < std::f64::consts::FRAC_PI_2,
                "Altitude out of range for {:?}: {}",
                body,
                out
            );
        }
    }

    // --- Azimuth tests ---

    #[test]
    fn body_azimuth_at_sun_is_finite() {
        let mut out = 0.0f64;
        let st = siderust_body_azimuth_at(SiderustBody::Sun, paris(), 60000.5, &mut out);
        assert_eq!(st, SiderustStatus::Ok);
        assert!(out.is_finite());
    }

    #[test]
    fn body_azimuth_at_mars_is_finite() {
        let mut out = 0.0f64;
        let st = siderust_body_azimuth_at(SiderustBody::Mars, paris(), 60000.5, &mut out);
        assert_eq!(st, SiderustStatus::Ok);
        assert!(out.is_finite());
    }

    #[test]
    fn all_bodies_return_finite_azimuth() {
        let bodies = [
            SiderustBody::Sun,
            SiderustBody::Moon,
            SiderustBody::Mercury,
            SiderustBody::Venus,
            SiderustBody::Mars,
            SiderustBody::Jupiter,
            SiderustBody::Saturn,
            SiderustBody::Uranus,
            SiderustBody::Neptune,
        ];
        for body in &bodies {
            let mut out = 0.0f64;
            let st = siderust_body_azimuth_at(*body, paris(), 60000.5, &mut out);
            assert_eq!(st, SiderustStatus::Ok, "azimuth_at failed for {:?}", body);
            assert!(
                out.is_finite(),
                "Non-finite azimuth for {:?}: {}",
                body,
                out
            );
            assert!(
                out >= 0.0 && out < std::f64::consts::TAU,
                "Azimuth out of range for {:?}: {}",
                body,
                out
            );
        }
    }

    #[test]
    fn body_azimuth_crossings_sun() {
        let mut out: *mut SiderustAzimuthCrossingEvent = ptr::null_mut();
        let mut count: usize = 0;
        let st = siderust_body_azimuth_crossings(
            SiderustBody::Sun,
            paris(),
            one_day_window(),
            180.0,
            default_opts(),
            &mut out,
            &mut count,
        );
        assert_eq!(st, SiderustStatus::Ok);
        if count > 0 {
            unsafe { crate::azimuth::siderust_azimuth_crossings_free(out, count) };
        }
    }
}
