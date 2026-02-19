// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! FFI bindings for ephemeris queries (VSOP87 / DE440 / DE441).

use crate::error::SiderustStatus;
use crate::types::*;
use siderust::calculus::ephemeris::{Ephemeris, Vsop87Ephemeris};
use siderust::time::JulianDate;

// ═══════════════════════════════════════════════════════════════════════════
// VSOP87 (always available)
// ═══════════════════════════════════════════════════════════════════════════

/// Get the Sun's barycentric position (EclipticMeanJ2000, AU) via VSOP87.
#[no_mangle]
pub extern "C" fn siderust_vsop87_sun_barycentric(
    jd: f64,
    out: *mut SiderustCartesianPos,
) -> SiderustStatus {
    if out.is_null() {
        return SiderustStatus::NullPointer;
    }
    let t = JulianDate::new(jd);
    let target = Vsop87Ephemeris::sun_barycentric(t);
    let pos = target.position;
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
}

/// Get the Earth's barycentric position (EclipticMeanJ2000, AU) via VSOP87.
#[no_mangle]
pub extern "C" fn siderust_vsop87_earth_barycentric(
    jd: f64,
    out: *mut SiderustCartesianPos,
) -> SiderustStatus {
    if out.is_null() {
        return SiderustStatus::NullPointer;
    }
    let t = JulianDate::new(jd);
    let target = Vsop87Ephemeris::earth_barycentric(t);
    let pos = target.position;
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
}

/// Get the Earth's heliocentric position (EclipticMeanJ2000, AU) via VSOP87.
#[no_mangle]
pub extern "C" fn siderust_vsop87_earth_heliocentric(
    jd: f64,
    out: *mut SiderustCartesianPos,
) -> SiderustStatus {
    if out.is_null() {
        return SiderustStatus::NullPointer;
    }
    let t = JulianDate::new(jd);
    let target = Vsop87Ephemeris::earth_heliocentric(t);
    let pos = target.position;
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
}

/// Get the Moon's geocentric position (EclipticMeanJ2000, km) via ELP2000.
#[no_mangle]
pub extern "C" fn siderust_vsop87_moon_geocentric(
    jd: f64,
    out: *mut SiderustCartesianPos,
) -> SiderustStatus {
    if out.is_null() {
        return SiderustStatus::NullPointer;
    }
    let t = JulianDate::new(jd);
    let pos = Vsop87Ephemeris::moon_geocentric(t);
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
}
