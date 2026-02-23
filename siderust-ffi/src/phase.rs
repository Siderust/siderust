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
    if out.is_null() {
        return SiderustStatus::NullPointer;
    }
    let geom = moon_phase_geocentric::<Vsop87Ephemeris>(JulianDate::new(jd));
    unsafe { *out = phase_geometry_from_rust(geom) };
    SiderustStatus::Ok
}

/// Compute topocentric Moon phase geometry at `jd` for `observer`.
#[no_mangle]
pub extern "C" fn siderust_moon_phase_topocentric(
    jd: f64,
    observer: SiderustGeodetict,
    out: *mut SiderustMoonPhaseGeometry,
) -> SiderustStatus {
    if out.is_null() {
        return SiderustStatus::NullPointer;
    }
    let geom = moon_phase_topocentric::<Vsop87Ephemeris>(JulianDate::new(jd), observer.to_rust());
    unsafe { *out = phase_geometry_from_rust(geom) };
    SiderustStatus::Ok
}

/// Map phase geometry to a named phase label.
#[no_mangle]
pub extern "C" fn siderust_moon_phase_label(
    geom: SiderustMoonPhaseGeometry,
    out: *mut SiderustMoonPhaseLabel,
) -> SiderustStatus {
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
    let window = match window_from_c(window) {
        Ok(w) => w,
        Err(e) => return e,
    };
    periods_to_c(
        illumination_above::<Vsop87Ephemeris>(window, k_min, search_opts_to_phase(opts)),
        out,
        count,
    )
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
    let window = match window_from_c(window) {
        Ok(w) => w,
        Err(e) => return e,
    };
    periods_to_c(
        illumination_below::<Vsop87Ephemeris>(window, k_max, search_opts_to_phase(opts)),
        out,
        count,
    )
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
    let window = match window_from_c(window) {
        Ok(w) => w,
        Err(e) => return e,
    };
    periods_to_c(
        illumination_range::<Vsop87Ephemeris>(window, k_min, k_max, search_opts_to_phase(opts)),
        out,
        count,
    )
}
