// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! FFI bindings for ephemeris queries (VSOP87 / DE440 / DE441).

use crate::error::SiderustStatus;
use crate::types::*;
use siderust::bodies::solar_system::{Jupiter, Mars, Mercury, Neptune, Saturn, Uranus, Venus};
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
    ffi_guard! {{
        if out.is_null() {
            return SiderustStatus::NullPointer;
        }
        let t = JulianDate::new(jd);
        let pos = Vsop87Ephemeris::sun_barycentric(t);
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

/// Get the Earth's barycentric position (EclipticMeanJ2000, AU) via VSOP87.
#[no_mangle]
pub extern "C" fn siderust_vsop87_earth_barycentric(
    jd: f64,
    out: *mut SiderustCartesianPos,
) -> SiderustStatus {
    ffi_guard! {{
        if out.is_null() {
            return SiderustStatus::NullPointer;
        }
        let t = JulianDate::new(jd);
        let pos = Vsop87Ephemeris::earth_barycentric(t);
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

/// Get the Earth's heliocentric position (EclipticMeanJ2000, AU) via VSOP87.
#[no_mangle]
pub extern "C" fn siderust_vsop87_earth_heliocentric(
    jd: f64,
    out: *mut SiderustCartesianPos,
) -> SiderustStatus {
    ffi_guard! {{
        if out.is_null() {
            return SiderustStatus::NullPointer;
        }
        let t = JulianDate::new(jd);
        let pos = Vsop87Ephemeris::earth_heliocentric(t);
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

/// Get Mars's heliocentric position (EclipticMeanJ2000, AU) via VSOP87.
#[no_mangle]
pub extern "C" fn siderust_vsop87_mars_heliocentric(
    jd: f64,
    out: *mut SiderustCartesianPos,
) -> SiderustStatus {
    ffi_guard! {{
        if out.is_null() {
            return SiderustStatus::NullPointer;
        }
        let t = JulianDate::new(jd);
        let pos = Mars::vsop87a(t);
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

/// Get Mars's barycentric position (EclipticMeanJ2000, AU) via VSOP87.
#[no_mangle]
pub extern "C" fn siderust_vsop87_mars_barycentric(
    jd: f64,
    out: *mut SiderustCartesianPos,
) -> SiderustStatus {
    ffi_guard! {{
        if out.is_null() {
            return SiderustStatus::NullPointer;
        }
        let t = JulianDate::new(jd);
        let pos = Mars::vsop87e(t);
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

/// Get Venus's heliocentric position (EclipticMeanJ2000, AU) via VSOP87.
#[no_mangle]
pub extern "C" fn siderust_vsop87_venus_heliocentric(
    jd: f64,
    out: *mut SiderustCartesianPos,
) -> SiderustStatus {
    ffi_guard! {{
        if out.is_null() {
            return SiderustStatus::NullPointer;
        }
        let t = JulianDate::new(jd);
        let pos = Venus::vsop87a(t);
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

// ── Remaining-planet VSOP87 helpers ─────────────────────────────────────────

/// Get Mercury's heliocentric position (EclipticMeanJ2000, AU) via VSOP87.
#[no_mangle]
pub extern "C" fn siderust_vsop87_mercury_heliocentric(
    jd: f64,
    out: *mut SiderustCartesianPos,
) -> SiderustStatus {
    ffi_guard! {{
        if out.is_null() { return SiderustStatus::NullPointer; }
        let t = JulianDate::new(jd);
        let pos = Mercury::vsop87a(t);
        unsafe { *out = SiderustCartesianPos { x: pos.x().value(), y: pos.y().value(), z: pos.z().value(), frame: SiderustFrame::EclipticMeanJ2000, center: SiderustCenter::Heliocentric }; }
        SiderustStatus::Ok

    }}
}

/// Get Mercury's barycentric position (EclipticMeanJ2000, AU) via VSOP87.
#[no_mangle]
pub extern "C" fn siderust_vsop87_mercury_barycentric(
    jd: f64,
    out: *mut SiderustCartesianPos,
) -> SiderustStatus {
    ffi_guard! {{
        if out.is_null() { return SiderustStatus::NullPointer; }
        let t = JulianDate::new(jd);
        let pos = Mercury::vsop87e(t);
        unsafe { *out = SiderustCartesianPos { x: pos.x().value(), y: pos.y().value(), z: pos.z().value(), frame: SiderustFrame::EclipticMeanJ2000, center: SiderustCenter::Barycentric }; }
        SiderustStatus::Ok

    }}
}

/// Venus barycentric VSOP87.
#[no_mangle]
pub extern "C" fn siderust_vsop87_venus_barycentric(
    jd: f64,
    out: *mut SiderustCartesianPos,
) -> SiderustStatus {
    ffi_guard! {{
        if out.is_null() { return SiderustStatus::NullPointer; }
        let t = JulianDate::new(jd);
        let pos = Venus::vsop87e(t);
        unsafe { *out = SiderustCartesianPos { x: pos.x().value(), y: pos.y().value(), z: pos.z().value(), frame: SiderustFrame::EclipticMeanJ2000, center: SiderustCenter::Barycentric }; }
        SiderustStatus::Ok

    }}
}

/// Get Jupiter's heliocentric position (EclipticMeanJ2000, AU) via VSOP87.
#[no_mangle]
pub extern "C" fn siderust_vsop87_jupiter_heliocentric(
    jd: f64,
    out: *mut SiderustCartesianPos,
) -> SiderustStatus {
    ffi_guard! {{
        if out.is_null() { return SiderustStatus::NullPointer; }
        let t = JulianDate::new(jd);
        let pos = Jupiter::vsop87a(t);
        unsafe { *out = SiderustCartesianPos { x: pos.x().value(), y: pos.y().value(), z: pos.z().value(), frame: SiderustFrame::EclipticMeanJ2000, center: SiderustCenter::Heliocentric }; }
        SiderustStatus::Ok

    }}
}

/// Get Jupiter's barycentric position (EclipticMeanJ2000, AU) via VSOP87.
#[no_mangle]
pub extern "C" fn siderust_vsop87_jupiter_barycentric(
    jd: f64,
    out: *mut SiderustCartesianPos,
) -> SiderustStatus {
    ffi_guard! {{
        if out.is_null() { return SiderustStatus::NullPointer; }
        let t = JulianDate::new(jd);
        let pos = Jupiter::vsop87e(t);
        unsafe { *out = SiderustCartesianPos { x: pos.x().value(), y: pos.y().value(), z: pos.z().value(), frame: SiderustFrame::EclipticMeanJ2000, center: SiderustCenter::Barycentric }; }
        SiderustStatus::Ok

    }}
}

/// Get Saturn's heliocentric position (EclipticMeanJ2000, AU) via VSOP87.
#[no_mangle]
pub extern "C" fn siderust_vsop87_saturn_heliocentric(
    jd: f64,
    out: *mut SiderustCartesianPos,
) -> SiderustStatus {
    ffi_guard! {{
        if out.is_null() { return SiderustStatus::NullPointer; }
        let t = JulianDate::new(jd);
        let pos = Saturn::vsop87a(t);
        unsafe { *out = SiderustCartesianPos { x: pos.x().value(), y: pos.y().value(), z: pos.z().value(), frame: SiderustFrame::EclipticMeanJ2000, center: SiderustCenter::Heliocentric }; }
        SiderustStatus::Ok

    }}
}

/// Get Saturn's barycentric position (EclipticMeanJ2000, AU) via VSOP87.
#[no_mangle]
pub extern "C" fn siderust_vsop87_saturn_barycentric(
    jd: f64,
    out: *mut SiderustCartesianPos,
) -> SiderustStatus {
    ffi_guard! {{
        if out.is_null() { return SiderustStatus::NullPointer; }
        let t = JulianDate::new(jd);
        let pos = Saturn::vsop87e(t);
        unsafe { *out = SiderustCartesianPos { x: pos.x().value(), y: pos.y().value(), z: pos.z().value(), frame: SiderustFrame::EclipticMeanJ2000, center: SiderustCenter::Barycentric }; }
        SiderustStatus::Ok

    }}
}

/// Get Uranus's heliocentric position (EclipticMeanJ2000, AU) via VSOP87.
#[no_mangle]
pub extern "C" fn siderust_vsop87_uranus_heliocentric(
    jd: f64,
    out: *mut SiderustCartesianPos,
) -> SiderustStatus {
    ffi_guard! {{
        if out.is_null() { return SiderustStatus::NullPointer; }
        let t = JulianDate::new(jd);
        let pos = Uranus::vsop87a(t);
        unsafe { *out = SiderustCartesianPos { x: pos.x().value(), y: pos.y().value(), z: pos.z().value(), frame: SiderustFrame::EclipticMeanJ2000, center: SiderustCenter::Heliocentric }; }
        SiderustStatus::Ok

    }}
}

/// Get Uranus's barycentric position (EclipticMeanJ2000, AU) via VSOP87.
#[no_mangle]
pub extern "C" fn siderust_vsop87_uranus_barycentric(
    jd: f64,
    out: *mut SiderustCartesianPos,
) -> SiderustStatus {
    ffi_guard! {{
        if out.is_null() { return SiderustStatus::NullPointer; }
        let t = JulianDate::new(jd);
        let pos = Uranus::vsop87e(t);
        unsafe { *out = SiderustCartesianPos { x: pos.x().value(), y: pos.y().value(), z: pos.z().value(), frame: SiderustFrame::EclipticMeanJ2000, center: SiderustCenter::Barycentric }; }
        SiderustStatus::Ok

    }}
}

/// Get Neptune's heliocentric position (EclipticMeanJ2000, AU) via VSOP87.
#[no_mangle]
pub extern "C" fn siderust_vsop87_neptune_heliocentric(
    jd: f64,
    out: *mut SiderustCartesianPos,
) -> SiderustStatus {
    ffi_guard! {{
        if out.is_null() { return SiderustStatus::NullPointer; }
        let t = JulianDate::new(jd);
        let pos = Neptune::vsop87a(t);
        unsafe { *out = SiderustCartesianPos { x: pos.x().value(), y: pos.y().value(), z: pos.z().value(), frame: SiderustFrame::EclipticMeanJ2000, center: SiderustCenter::Heliocentric }; }
        SiderustStatus::Ok

    }}
}

/// Get Neptune's barycentric position (EclipticMeanJ2000, AU) via VSOP87.
#[no_mangle]
pub extern "C" fn siderust_vsop87_neptune_barycentric(
    jd: f64,
    out: *mut SiderustCartesianPos,
) -> SiderustStatus {
    ffi_guard! {{
        if out.is_null() { return SiderustStatus::NullPointer; }
        let t = JulianDate::new(jd);
        let pos = Neptune::vsop87e(t);
        unsafe { *out = SiderustCartesianPos { x: pos.x().value(), y: pos.y().value(), z: pos.z().value(), frame: SiderustFrame::EclipticMeanJ2000, center: SiderustCenter::Barycentric }; }
        SiderustStatus::Ok

    }}
}

/// Get the Moon's geocentric position (EclipticMeanJ2000, km) via ELP2000.
#[no_mangle]
pub extern "C" fn siderust_vsop87_moon_geocentric(
    jd: f64,
    out: *mut SiderustCartesianPos,
) -> SiderustStatus {
    ffi_guard! {{
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

    }}
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::ptr;

    const J2000: f64 = 2_451_545.0;

    fn empty_pos() -> SiderustCartesianPos {
        SiderustCartesianPos {
            x: 0.0,
            y: 0.0,
            z: 0.0,
            frame: SiderustFrame::ICRS,
            center: SiderustCenter::Barycentric,
        }
    }

    #[test]
    fn sun_barycentric_at_j2000() {
        let mut out = empty_pos();
        let s = siderust_vsop87_sun_barycentric(J2000, &mut out);
        assert_eq!(s, SiderustStatus::Ok);
        // Sun barycentric at J2000 is close to origin (helio ~ bary for sun)
        let dist = (out.x * out.x + out.y * out.y + out.z * out.z).sqrt();
        assert!(dist < 0.01, "Sun barycentric too far from SSB: {dist} AU");
        assert_eq!(out.frame, SiderustFrame::EclipticMeanJ2000);
        assert_eq!(out.center, SiderustCenter::Barycentric);
    }

    #[test]
    fn sun_barycentric_null_ptr() {
        assert_eq!(
            siderust_vsop87_sun_barycentric(J2000, ptr::null_mut()),
            SiderustStatus::NullPointer
        );
    }

    #[test]
    fn earth_barycentric_at_j2000() {
        let mut out = empty_pos();
        let s = siderust_vsop87_earth_barycentric(J2000, &mut out);
        assert_eq!(s, SiderustStatus::Ok);
        // Earth is ~1 AU from SSB
        let dist = (out.x * out.x + out.y * out.y + out.z * out.z).sqrt();
        assert!(
            dist > 0.9 && dist < 1.1,
            "Earth should be ~1 AU from SSB, got {dist}"
        );
    }

    #[test]
    fn earth_barycentric_null_ptr() {
        assert_eq!(
            siderust_vsop87_earth_barycentric(J2000, ptr::null_mut()),
            SiderustStatus::NullPointer
        );
    }

    #[test]
    fn earth_heliocentric_at_j2000() {
        let mut out = empty_pos();
        let s = siderust_vsop87_earth_heliocentric(J2000, &mut out);
        assert_eq!(s, SiderustStatus::Ok);
        let dist = (out.x * out.x + out.y * out.y + out.z * out.z).sqrt();
        assert!(
            dist > 0.98 && dist < 1.02,
            "Earth heliocentric ~1 AU, got {dist}"
        );
        assert_eq!(out.center, SiderustCenter::Heliocentric);
    }

    #[test]
    fn earth_heliocentric_null_ptr() {
        assert_eq!(
            siderust_vsop87_earth_heliocentric(J2000, ptr::null_mut()),
            SiderustStatus::NullPointer
        );
    }

    #[test]
    fn moon_geocentric_at_j2000() {
        let mut out = empty_pos();
        let s = siderust_vsop87_moon_geocentric(J2000, &mut out);
        assert_eq!(s, SiderustStatus::Ok);
        // Moon is ~384,400 km from Earth
        let dist = (out.x * out.x + out.y * out.y + out.z * out.z).sqrt();
        assert!(
            dist > 350_000.0 && dist < 420_000.0,
            "Moon geocentric ~384400 km, got {dist}"
        );
        assert_eq!(out.center, SiderustCenter::Geocentric);
    }

    #[test]
    fn moon_geocentric_null_ptr() {
        assert_eq!(
            siderust_vsop87_moon_geocentric(J2000, ptr::null_mut()),
            SiderustStatus::NullPointer
        );
    }
}
