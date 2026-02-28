// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Lunar phase FFI — phase geometry, event finding, illumination periods.

use crate::altitude::{periods_to_c, window_from_c};
use crate::error::SiderustStatus;
use crate::ffi_utils::{free_boxed_slice, vec_to_c};
use crate::types::*;
use qtty::*;
use siderust::calculus::ephemeris::Vsop87Ephemeris;
use siderust::calculus::lunar::phase::{
    find_phase_events, illumination_above, illumination_below, illumination_range,
    moon_phase_geocentric, moon_phase_topocentric, MoonPhaseLabel, PhaseKind, PhaseSearchOpts,
};
use siderust::time::JulianDate;

// ═══════════════════════════════════════════════════════════════════════════
// Conversion helpers
// ═══════════════════════════════════════════════════════════════════════════

fn phase_geometry_from_rust(
    g: siderust::calculus::lunar::phase::MoonPhaseGeometry,
) -> SiderustMoonPhaseGeometry {
    SiderustMoonPhaseGeometry {
        phase_angle_rad: g.phase_angle.value(),
        illuminated_fraction: g.illuminated_fraction,
        elongation_rad: g.elongation.value(),
        waxing: if g.waxing { 1 } else { 0 },
        _pad: [0; 7],
    }
}

fn phase_kind_from_rust(k: PhaseKind) -> SiderustPhaseKind {
    match k {
        PhaseKind::NewMoon => SiderustPhaseKind::NewMoon,
        PhaseKind::FirstQuarter => SiderustPhaseKind::FirstQuarter,
        PhaseKind::FullMoon => SiderustPhaseKind::FullMoon,
        PhaseKind::LastQuarter => SiderustPhaseKind::LastQuarter,
    }
}

fn search_opts_to_phase(opts: SiderustSearchOpts) -> PhaseSearchOpts {
    PhaseSearchOpts {
        time_tolerance: Days::new(opts.time_tolerance_days),
        scan_step: if opts.has_scan_step && opts.scan_step_days > 0.0 {
            Days::new(opts.scan_step_days)
        } else {
            PhaseSearchOpts::default().scan_step
        },
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// Memory management
// ═══════════════════════════════════════════════════════════════════════════

/// Free an array of phase events.
#[no_mangle]
pub unsafe extern "C" fn siderust_phase_events_free(ptr: *mut SiderustPhaseEvent, count: usize) {
    free_boxed_slice(ptr, count);
}

// ═══════════════════════════════════════════════════════════════════════════
// Phase geometry
// ═══════════════════════════════════════════════════════════════════════════

/// Compute geocentric Moon phase geometry at `jd` (Julian Date).
#[no_mangle]
pub extern "C" fn siderust_moon_phase_geocentric(
    jd: f64,
    out: *mut SiderustMoonPhaseGeometry,
) -> SiderustStatus {
    ffi_guard! {{
        if out.is_null() {
            return SiderustStatus::NullPointer;
        }
        let geom = moon_phase_geocentric::<Vsop87Ephemeris>(JulianDate::new(jd));
        unsafe { *out = phase_geometry_from_rust(geom) };
        SiderustStatus::Ok

    }}
}

/// Compute topocentric Moon phase geometry at `jd` for `observer`.
#[no_mangle]
pub extern "C" fn siderust_moon_phase_topocentric(
    jd: f64,
    observer: SiderustGeodetict,
    out: *mut SiderustMoonPhaseGeometry,
) -> SiderustStatus {
    ffi_guard! {{
        if out.is_null() {
            return SiderustStatus::NullPointer;
        }
        let geom = moon_phase_topocentric::<Vsop87Ephemeris>(JulianDate::new(jd), observer.to_rust());
        unsafe { *out = phase_geometry_from_rust(geom) };
        SiderustStatus::Ok

    }}
}

/// Map phase geometry to a named phase label.
#[no_mangle]
pub extern "C" fn siderust_moon_phase_label(
    geom: SiderustMoonPhaseGeometry,
    out: *mut SiderustMoonPhaseLabel,
) -> SiderustStatus {
    ffi_guard! {{
        if out.is_null() {
            return SiderustStatus::NullPointer;
        }
        let elong_deg = Radians::new(geom.elongation_rad).to::<Degree>();
        let label = MoonPhaseLabel::from_elongation(elong_deg, &Default::default());
        let ffi_label = match label {
            MoonPhaseLabel::NewMoon => SiderustMoonPhaseLabel::NewMoon,
            MoonPhaseLabel::WaxingCrescent => SiderustMoonPhaseLabel::WaxingCrescent,
            MoonPhaseLabel::FirstQuarter => SiderustMoonPhaseLabel::FirstQuarter,
            MoonPhaseLabel::WaxingGibbous => SiderustMoonPhaseLabel::WaxingGibbous,
            MoonPhaseLabel::FullMoon => SiderustMoonPhaseLabel::FullMoon,
            MoonPhaseLabel::WaningGibbous => SiderustMoonPhaseLabel::WaningGibbous,
            MoonPhaseLabel::LastQuarter => SiderustMoonPhaseLabel::LastQuarter,
            MoonPhaseLabel::WaningCrescent => SiderustMoonPhaseLabel::WaningCrescent,
        };
        unsafe { *out = ffi_label };
        SiderustStatus::Ok

    }}
}

// ═══════════════════════════════════════════════════════════════════════════
// Phase event finding
// ═══════════════════════════════════════════════════════════════════════════

/// Find all principal lunar phase events in `window`.
///
/// Results are sorted chronologically. The caller must free the array
/// with [`siderust_phase_events_free`].
#[no_mangle]
pub extern "C" fn siderust_find_phase_events(
    window: TempochPeriodMjd,
    opts: SiderustSearchOpts,
    out: *mut *mut SiderustPhaseEvent,
    count: *mut usize,
) -> SiderustStatus {
    ffi_guard! {{
        let window = match window_from_c(window) {
            Ok(w) => w,
            Err(e) => return e,
        };
        let events = find_phase_events::<Vsop87Ephemeris>(window, search_opts_to_phase(opts));
        vec_to_c(
            events,
            |e| SiderustPhaseEvent {
                mjd: e.mjd.value(),
                kind: phase_kind_from_rust(e.kind),
                _pad: [0; 4],
            },
            out,
            count,
        )

    }}
}

// ═══════════════════════════════════════════════════════════════════════════
// Illumination period finding
// ═══════════════════════════════════════════════════════════════════════════

/// Find windows where geocentric Moon illumination is above `k_min` ∈ [0,1].
#[no_mangle]
pub extern "C" fn siderust_moon_illumination_above(
    window: TempochPeriodMjd,
    k_min: f64,
    opts: SiderustSearchOpts,
    out: *mut *mut TempochPeriodMjd,
    count: *mut usize,
) -> SiderustStatus {
    ffi_guard! {{
        let window = match window_from_c(window) {
            Ok(w) => w,
            Err(e) => return e,
        };
        periods_to_c(
            illumination_above::<Vsop87Ephemeris>(window, k_min, search_opts_to_phase(opts)),
            out,
            count,
        )

    }}
}

/// Find windows where geocentric Moon illumination is below `k_max` ∈ [0,1].
#[no_mangle]
pub extern "C" fn siderust_moon_illumination_below(
    window: TempochPeriodMjd,
    k_max: f64,
    opts: SiderustSearchOpts,
    out: *mut *mut TempochPeriodMjd,
    count: *mut usize,
) -> SiderustStatus {
    ffi_guard! {{
        let window = match window_from_c(window) {
            Ok(w) => w,
            Err(e) => return e,
        };
        periods_to_c(
            illumination_below::<Vsop87Ephemeris>(window, k_max, search_opts_to_phase(opts)),
            out,
            count,
        )

    }}
}

/// Find windows where Moon illumination is within [k_min, k_max].
#[no_mangle]
pub extern "C" fn siderust_moon_illumination_range(
    window: TempochPeriodMjd,
    k_min: f64,
    k_max: f64,
    opts: SiderustSearchOpts,
    out: *mut *mut TempochPeriodMjd,
    count: *mut usize,
) -> SiderustStatus {
    ffi_guard! {{
        let window = match window_from_c(window) {
            Ok(w) => w,
            Err(e) => return e,
        };
        periods_to_c(
            illumination_range::<Vsop87Ephemeris>(window, k_min, k_max, search_opts_to_phase(opts)),
            out,
            count,
        )

    }}
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::ptr;

    fn one_month() -> TempochPeriodMjd {
        TempochPeriodMjd {
            start_mjd: 60000.0,
            end_mjd: 60030.0,
        }
    }

    fn default_opts() -> SiderustSearchOpts {
        SiderustSearchOpts {
            time_tolerance_days: 1e-9,
            scan_step_days: 0.0,
            has_scan_step: false,
        }
    }

    fn paris() -> SiderustGeodetict {
        SiderustGeodetict {
            lon_deg: 2.35,
            lat_deg: 48.85,
            height_m: 35.0,
        }
    }

    // ── Geocentric phase ──────────────────────────────────────────────────

    #[test]
    fn geocentric_phase_null_out() {
        assert_eq!(
            siderust_moon_phase_geocentric(2451545.0, ptr::null_mut()),
            SiderustStatus::NullPointer
        );
    }

    #[test]
    fn geocentric_phase_valid() {
        let mut out = SiderustMoonPhaseGeometry {
            phase_angle_rad: 0.0,
            illuminated_fraction: 0.0,
            elongation_rad: 0.0,
            waxing: 0,
            _pad: [0; 7],
        };
        assert_eq!(
            siderust_moon_phase_geocentric(2451545.0, &mut out),
            SiderustStatus::Ok
        );
        assert!(out.phase_angle_rad.is_finite());
        assert!(out.illuminated_fraction >= 0.0 && out.illuminated_fraction <= 1.0);
        assert!(out.elongation_rad.is_finite());
    }

    // ── Topocentric phase ─────────────────────────────────────────────────

    #[test]
    fn topocentric_phase_null_out() {
        assert_eq!(
            siderust_moon_phase_topocentric(2451545.0, paris(), ptr::null_mut()),
            SiderustStatus::NullPointer
        );
    }

    #[test]
    fn topocentric_phase_valid() {
        let mut out = SiderustMoonPhaseGeometry {
            phase_angle_rad: 0.0,
            illuminated_fraction: 0.0,
            elongation_rad: 0.0,
            waxing: 0,
            _pad: [0; 7],
        };
        assert_eq!(
            siderust_moon_phase_topocentric(2451545.0, paris(), &mut out),
            SiderustStatus::Ok
        );
        assert!(out.illuminated_fraction >= 0.0 && out.illuminated_fraction <= 1.0);
    }

    // ── Phase label ───────────────────────────────────────────────────────

    #[test]
    fn phase_label_null_out() {
        let geom = SiderustMoonPhaseGeometry {
            phase_angle_rad: 0.0,
            illuminated_fraction: 0.0,
            elongation_rad: 0.0,
            waxing: 0,
            _pad: [0; 7],
        };
        assert_eq!(
            siderust_moon_phase_label(geom, ptr::null_mut()),
            SiderustStatus::NullPointer
        );
    }

    #[test]
    fn phase_label_new_moon() {
        let mut geom = SiderustMoonPhaseGeometry {
            phase_angle_rad: 0.0,
            illuminated_fraction: 0.0,
            elongation_rad: 0.0,
            waxing: 1,
            _pad: [0; 7],
        };
        let mut label = SiderustMoonPhaseLabel::NewMoon;
        assert_eq!(
            siderust_moon_phase_geocentric(2451545.0, &mut geom),
            SiderustStatus::Ok
        );
        assert_eq!(
            siderust_moon_phase_label(geom, &mut label),
            SiderustStatus::Ok
        );
    }

    // ── Phase events ──────────────────────────────────────────────────────

    #[test]
    fn find_phase_events_null_out() {
        let mut count = 0usize;
        assert_eq!(
            siderust_find_phase_events(one_month(), default_opts(), ptr::null_mut(), &mut count),
            SiderustStatus::NullPointer
        );
    }

    #[test]
    fn find_phase_events_null_count() {
        let mut out: *mut SiderustPhaseEvent = ptr::null_mut();
        assert_eq!(
            siderust_find_phase_events(one_month(), default_opts(), &mut out, ptr::null_mut()),
            SiderustStatus::NullPointer
        );
    }

    #[test]
    fn find_phase_events_returns_events() {
        let mut out: *mut SiderustPhaseEvent = ptr::null_mut();
        let mut count = 0usize;
        assert_eq!(
            siderust_find_phase_events(one_month(), default_opts(), &mut out, &mut count),
            SiderustStatus::Ok
        );
        // 30-day window must contain at least one principal phase
        assert!(count >= 1);
        unsafe { siderust_phase_events_free(out, count) };
    }

    #[test]
    fn find_phase_events_invalid_window() {
        let bad = TempochPeriodMjd {
            start_mjd: 60030.0,
            end_mjd: 60000.0,
        };
        let mut out: *mut SiderustPhaseEvent = ptr::null_mut();
        let mut count = 0usize;
        assert_eq!(
            siderust_find_phase_events(bad, default_opts(), &mut out, &mut count),
            SiderustStatus::InvalidPeriod
        );
    }

    // ── Illumination above ────────────────────────────────────────────────

    #[test]
    fn illumination_above_ok() {
        let mut out: *mut TempochPeriodMjd = ptr::null_mut();
        let mut count = 0usize;
        assert_eq!(
            siderust_moon_illumination_above(
                one_month(),
                0.5,
                default_opts(),
                &mut out,
                &mut count
            ),
            SiderustStatus::Ok
        );
        unsafe { crate::altitude::siderust_periods_free(out, count) };
    }

    #[test]
    fn illumination_above_invalid_window() {
        let bad = TempochPeriodMjd {
            start_mjd: 60030.0,
            end_mjd: 60000.0,
        };
        let mut out: *mut TempochPeriodMjd = ptr::null_mut();
        let mut count = 0usize;
        assert_eq!(
            siderust_moon_illumination_above(bad, 0.5, default_opts(), &mut out, &mut count),
            SiderustStatus::InvalidPeriod
        );
    }

    // ── Illumination below ────────────────────────────────────────────────

    #[test]
    fn illumination_below_ok() {
        let mut out: *mut TempochPeriodMjd = ptr::null_mut();
        let mut count = 0usize;
        assert_eq!(
            siderust_moon_illumination_below(
                one_month(),
                0.5,
                default_opts(),
                &mut out,
                &mut count
            ),
            SiderustStatus::Ok
        );
        unsafe { crate::altitude::siderust_periods_free(out, count) };
    }

    #[test]
    fn illumination_below_invalid_window() {
        let bad = TempochPeriodMjd {
            start_mjd: 60030.0,
            end_mjd: 60000.0,
        };
        let mut out: *mut TempochPeriodMjd = ptr::null_mut();
        let mut count = 0usize;
        assert_eq!(
            siderust_moon_illumination_below(bad, 0.5, default_opts(), &mut out, &mut count),
            SiderustStatus::InvalidPeriod
        );
    }

    // ── Illumination range ────────────────────────────────────────────────

    #[test]
    fn illumination_range_ok() {
        let mut out: *mut TempochPeriodMjd = ptr::null_mut();
        let mut count = 0usize;
        assert_eq!(
            siderust_moon_illumination_range(
                one_month(),
                0.25,
                0.75,
                default_opts(),
                &mut out,
                &mut count
            ),
            SiderustStatus::Ok
        );
        unsafe { crate::altitude::siderust_periods_free(out, count) };
    }

    #[test]
    fn illumination_range_invalid_window() {
        let bad = TempochPeriodMjd {
            start_mjd: 60030.0,
            end_mjd: 60000.0,
        };
        let mut out: *mut TempochPeriodMjd = ptr::null_mut();
        let mut count = 0usize;
        assert_eq!(
            siderust_moon_illumination_range(bad, 0.25, 0.75, default_opts(), &mut out, &mut count),
            SiderustStatus::InvalidPeriod
        );
    }

    // ── Free safety ───────────────────────────────────────────────────────

    #[test]
    fn phase_events_free_null_safe() {
        unsafe { siderust_phase_events_free(ptr::null_mut(), 0) };
    }
}
