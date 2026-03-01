// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! FFI bindings for **generic body dispatch**: altitude, azimuth, crossings,
//! and culminations for any solar-system body identified by [`SiderustBody`].
//!
//! Every function here is a thin wrapper that constructs a
//! [`SiderustSubject::body`] and delegates to the corresponding unified
//! function in [`crate::subject`].

use crate::error::SiderustStatus;
use crate::types::*;

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
    crate::subject::siderust_altitude_at(SiderustSubject::body(body), observer, mjd, out_rad)
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
    crate::subject::siderust_above_threshold(
        SiderustSubject::body(body),
        observer,
        window,
        threshold_deg,
        opts,
        out,
        count,
    )
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
    crate::subject::siderust_below_threshold(
        SiderustSubject::body(body),
        observer,
        window,
        threshold_deg,
        opts,
        out,
        count,
    )
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
    crate::subject::siderust_crossings(
        SiderustSubject::body(body),
        observer,
        window,
        threshold_deg,
        opts,
        out,
        count,
    )
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
    crate::subject::siderust_culminations(
        SiderustSubject::body(body),
        observer,
        window,
        opts,
        out,
        count,
    )
}

/// Periods when a solar-system body's altitude is within [min, max].
#[no_mangle]
pub extern "C" fn siderust_body_altitude_periods(
    body: SiderustBody,
    query: SiderustAltitudeQuery,
    out: *mut *mut TempochPeriodMjd,
    count: *mut usize,
) -> SiderustStatus {
    crate::subject::siderust_altitude_periods(SiderustSubject::body(body), query, out, count)
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
    crate::subject::siderust_azimuth_at(SiderustSubject::body(body), observer, mjd, out_rad)
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
    crate::subject::siderust_azimuth_crossings(
        SiderustSubject::body(body),
        observer,
        window,
        bearing_deg,
        opts,
        out,
        count,
    )
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
    crate::subject::siderust_azimuth_extrema(
        SiderustSubject::body(body),
        observer,
        window,
        opts,
        out,
        count,
    )
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
    crate::subject::siderust_in_azimuth_range(
        SiderustSubject::body(body),
        observer,
        window,
        min_deg,
        max_deg,
        opts,
        out,
        count,
    )
}

// ═══════════════════════════════════════════════════════════════════════════
// Tests
// ═══════════════════════════════════════════════════════════════════════════

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_helpers::*;
    use std::ptr;

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

    // ── body_below_threshold ──────────────────────────────────────────────

    #[test]
    fn body_below_threshold_sun_succeeds() {
        let mut out: *mut TempochPeriodMjd = ptr::null_mut();
        let mut count: usize = 0;
        let st = siderust_body_below_threshold(
            SiderustBody::Sun,
            paris(),
            one_day_window(),
            0.0,
            default_opts(),
            &mut out,
            &mut count,
        );
        assert_eq!(st, SiderustStatus::Ok);
        if count > 0 {
            unsafe { crate::altitude::siderust_periods_free(out, count) };
        }
    }

    #[test]
    fn body_below_threshold_moon_returns_ok() {
        let mut out: *mut TempochPeriodMjd = ptr::null_mut();
        let mut count: usize = 0;
        let window = TempochPeriodMjd {
            start_mjd: 60000.0,
            end_mjd: 60003.0,
        };
        let st = siderust_body_below_threshold(
            SiderustBody::Moon,
            paris(),
            window,
            0.0,
            default_opts(),
            &mut out,
            &mut count,
        );
        assert_eq!(st, SiderustStatus::Ok);
        if count > 0 {
            unsafe { crate::altitude::siderust_periods_free(out, count) };
        }
    }

    // ── body_altitude_periods ─────────────────────────────────────────────

    #[test]
    fn body_altitude_periods_sun_succeeds() {
        let mut out: *mut TempochPeriodMjd = ptr::null_mut();
        let mut count: usize = 0;
        let query = SiderustAltitudeQuery {
            observer: paris(),
            start_mjd: 60000.0,
            end_mjd: 60001.0,
            min_altitude_deg: 0.0,
            max_altitude_deg: 90.0,
        };
        let st = siderust_body_altitude_periods(SiderustBody::Sun, query, &mut out, &mut count);
        assert_eq!(st, SiderustStatus::Ok);
        if count > 0 {
            unsafe { crate::altitude::siderust_periods_free(out, count) };
        }
    }

    #[test]
    fn body_altitude_periods_jupiter_succeeds() {
        let mut out: *mut TempochPeriodMjd = ptr::null_mut();
        let mut count: usize = 0;
        let query = SiderustAltitudeQuery {
            observer: paris(),
            start_mjd: 60000.0,
            end_mjd: 60007.0,
            min_altitude_deg: 10.0,
            max_altitude_deg: 90.0,
        };
        let st = siderust_body_altitude_periods(SiderustBody::Jupiter, query, &mut out, &mut count);
        assert_eq!(st, SiderustStatus::Ok);
        if count > 0 {
            unsafe { crate::altitude::siderust_periods_free(out, count) };
        }
    }

    // ── body_azimuth_extrema ──────────────────────────────────────────────

    #[test]
    fn body_azimuth_extrema_sun_succeeds() {
        let mut out: *mut SiderustAzimuthExtremum = ptr::null_mut();
        let mut count: usize = 0;
        let st = siderust_body_azimuth_extrema(
            SiderustBody::Sun,
            paris(),
            one_day_window(),
            default_opts(),
            &mut out,
            &mut count,
        );
        assert_eq!(st, SiderustStatus::Ok);
        if count > 0 {
            unsafe { crate::azimuth::siderust_azimuth_extrema_free(out, count) };
        }
    }

    #[test]
    fn body_azimuth_extrema_moon_succeeds() {
        let mut out: *mut SiderustAzimuthExtremum = ptr::null_mut();
        let mut count: usize = 0;
        let window = TempochPeriodMjd {
            start_mjd: 60000.0,
            end_mjd: 60003.0,
        };
        let st = siderust_body_azimuth_extrema(
            SiderustBody::Moon,
            paris(),
            window,
            default_opts(),
            &mut out,
            &mut count,
        );
        assert_eq!(st, SiderustStatus::Ok);
        if count > 0 {
            unsafe { crate::azimuth::siderust_azimuth_extrema_free(out, count) };
        }
    }

    // ── body_in_azimuth_range ─────────────────────────────────────────────

    #[test]
    fn body_in_azimuth_range_sun_succeeds() {
        let mut out: *mut TempochPeriodMjd = ptr::null_mut();
        let mut count: usize = 0;
        let st = siderust_body_in_azimuth_range(
            SiderustBody::Sun,
            paris(),
            one_day_window(),
            90.0,  // east
            270.0, // west
            default_opts(),
            &mut out,
            &mut count,
        );
        assert_eq!(st, SiderustStatus::Ok);
        if count > 0 {
            unsafe { crate::altitude::siderust_periods_free(out, count) };
        }
    }

    #[test]
    fn body_in_azimuth_range_mars_succeeds() {
        let mut out: *mut TempochPeriodMjd = ptr::null_mut();
        let mut count: usize = 0;
        let window = TempochPeriodMjd {
            start_mjd: 60000.0,
            end_mjd: 60007.0,
        };
        let st = siderust_body_in_azimuth_range(
            SiderustBody::Mars,
            paris(),
            window,
            0.0,
            360.0,
            default_opts(),
            &mut out,
            &mut count,
        );
        assert_eq!(st, SiderustStatus::Ok);
        if count > 0 {
            unsafe { crate::altitude::siderust_periods_free(out, count) };
        }
    }
}
