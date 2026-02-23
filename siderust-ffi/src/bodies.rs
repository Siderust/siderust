// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! FFI bindings for celestial bodies — Sun, Moon, Star, Planet, etc.

use crate::error::SiderustStatus;
use crate::types::{SiderustPlanet, SiderustProperMotion, SiderustRaConvention};
use qtty::length::nominal::SolarRadiuses;
use qtty::*;
use siderust::astro::proper_motion::{ProperMotion, RaProperMotionConvention};
use siderust::bodies::{self, Star};
use siderust::coordinates::centers::Geocentric;
use siderust::coordinates::frames::EquatorialMeanJ2000;
use siderust::coordinates::spherical;
use siderust::targets::Target;
use siderust::time::JulianDate;
use std::borrow::Cow;
use std::ffi::CStr;
use std::os::raw::c_char;

// ═══════════════════════════════════════════════════════════════════════════
// Opaque handle for Star
// ═══════════════════════════════════════════════════════════════════════════

/// Opaque handle to a Star. Created via `siderust_star_*` functions, freed
/// with `siderust_star_free`.
pub struct SiderustStar {
    pub(crate) inner: Star<'static>,
}

/// Look up a star from the built-in catalog by name.
///
/// Supported names: "VEGA", "SIRIUS", "POLARIS", "CANOPUS", "ARCTURUS",
/// "RIGEL", "BETELGEUSE", "PROCYON", "ALDEBARAN", "ALTAIR".
///
/// On success, `*out` receives a newly allocated handle.
/// The caller must free it with `siderust_star_free`.
#[no_mangle]
pub extern "C" fn siderust_star_catalog(
    name: *const c_char,
    out: *mut *mut SiderustStar,
) -> SiderustStatus {
    if name.is_null() || out.is_null() {
        return SiderustStatus::NullPointer;
    }
    let name_str = unsafe { CStr::from_ptr(name) };
    let name_str = match name_str.to_str() {
        Ok(s) => s,
        Err(_) => return SiderustStatus::InvalidArgument,
    };

    let star: Option<&Star<'static>> = match name_str.to_uppercase().as_str() {
        "VEGA" => Some(&bodies::VEGA),
        "SIRIUS" => Some(&bodies::SIRIUS),
        "POLARIS" => Some(&bodies::POLARIS),
        "CANOPUS" => Some(&bodies::CANOPUS),
        "ARCTURUS" => Some(&bodies::ARCTURUS),
        "RIGEL" => Some(&bodies::RIGEL),
        "BETELGEUSE" => Some(&bodies::BETELGEUSE),
        "PROCYON" => Some(&bodies::PROCYON),
        "ALDEBARAN" => Some(&bodies::ALDEBARAN),
        "ALTAIR" => Some(&bodies::ALTAIR),
        _ => None,
    };

    match star {
        Some(s) => {
            let handle = Box::new(SiderustStar { inner: s.clone() });
            unsafe { *out = Box::into_raw(handle) };
            SiderustStatus::Ok
        }
        None => SiderustStatus::UnknownStar,
    }
}

/// Create a custom star.
///
/// `name` — UTF-8 null-terminated string.
/// `ra_deg`, `dec_deg` — J2000 equatorial coordinates in degrees.
/// `distance_ly` — distance in light-years.
/// `epoch_jd` — epoch of the given coordinates (Julian Date).
/// `proper_motion` — optional: pass null for a fixed star.
#[no_mangle]
pub extern "C" fn siderust_star_create(
    name: *const c_char,
    distance_ly: f64,
    mass_solar: f64,
    radius_solar: f64,
    luminosity_solar: f64,
    ra_deg: f64,
    dec_deg: f64,
    epoch_jd: f64,
    proper_motion: *const SiderustProperMotion,
    out: *mut *mut SiderustStar,
) -> SiderustStatus {
    if name.is_null() || out.is_null() {
        return SiderustStatus::NullPointer;
    }
    let name_str = unsafe { CStr::from_ptr(name) };
    let name_str = match name_str.to_str() {
        Ok(s) => s.to_owned(),
        Err(_) => return SiderustStatus::InvalidArgument,
    };

    let pos = spherical::Position::<Geocentric, EquatorialMeanJ2000, LightYear>::new(
        Degrees::new(ra_deg),
        Degrees::new(dec_deg),
        LightYears::new(distance_ly),
    );

    let epoch = JulianDate::new(epoch_jd);

    let pm = if proper_motion.is_null() {
        None
    } else {
        let pm_c = unsafe { &*proper_motion };
        let convention = match pm_c.ra_convention {
            SiderustRaConvention::MuAlpha => RaProperMotionConvention::MuAlpha,
            SiderustRaConvention::MuAlphaStar => RaProperMotionConvention::MuAlphaStar,
        };
        Some(ProperMotion {
            pm_ra: Quantity::new(pm_c.pm_ra_deg_yr),
            pm_dec: Quantity::new(pm_c.pm_dec_deg_yr),
            ra_convention: convention,
        })
    };

    let target = match pm {
        Some(pm) => Target::new(pos, epoch, pm),
        None => Target::new_static(pos, epoch),
    };

    let star = Star::new(
        Cow::<'static, str>::Owned(name_str),
        LightYears::new(distance_ly),
        SolarMasses::new(mass_solar),
        SolarRadiuses::new(radius_solar),
        SolarLuminosities::new(luminosity_solar),
        target,
    );

    let handle = Box::new(SiderustStar { inner: star });
    unsafe { *out = Box::into_raw(handle) };
    SiderustStatus::Ok
}

/// Free a Star handle.
///
/// # Safety
/// The handle must have been allocated by `siderust_star_catalog` or
/// `siderust_star_create`, and must not be used after this call.
#[no_mangle]
pub unsafe extern "C" fn siderust_star_free(handle: *mut SiderustStar) {
    if !handle.is_null() {
        drop(Box::from_raw(handle));
    }
}

/// Get the star's name. Copies into `buf` up to `buf_len` bytes (including NUL).
/// Sets `*written` to the number of bytes written (excluding NUL).
/// Returns InvalidArgument if the buffer is too small.
#[no_mangle]
pub extern "C" fn siderust_star_name(
    handle: *const SiderustStar,
    buf: *mut c_char,
    buf_len: usize,
    written: *mut usize,
) -> SiderustStatus {
    if handle.is_null() || buf.is_null() {
        return SiderustStatus::NullPointer;
    }
    let star = unsafe { &*handle };
    let name = star.inner.name.as_bytes();
    if name.len() + 1 > buf_len {
        return SiderustStatus::InvalidArgument;
    }
    unsafe {
        std::ptr::copy_nonoverlapping(name.as_ptr(), buf as *mut u8, name.len());
        *buf.add(name.len()) = 0; // NUL terminator
        if !written.is_null() {
            *written = name.len();
        }
    }
    SiderustStatus::Ok
}

/// Get the star's distance in light-years.
#[no_mangle]
pub extern "C" fn siderust_star_distance_ly(handle: *const SiderustStar) -> f64 {
    if handle.is_null() {
        return f64::NAN;
    }
    unsafe { (*handle).inner.distance.value() }
}

/// Get the star's mass in solar masses.
#[no_mangle]
pub extern "C" fn siderust_star_mass_solar(handle: *const SiderustStar) -> f64 {
    if handle.is_null() {
        return f64::NAN;
    }
    unsafe { (*handle).inner.mass.value() }
}

/// Get the star's radius in solar radii.
#[no_mangle]
pub extern "C" fn siderust_star_radius_solar(handle: *const SiderustStar) -> f64 {
    if handle.is_null() {
        return f64::NAN;
    }
    unsafe { (*handle).inner.radius.value() }
}

/// Get the star's luminosity in solar luminosities.
#[no_mangle]
pub extern "C" fn siderust_star_luminosity_solar(handle: *const SiderustStar) -> f64 {
    if handle.is_null() {
        return f64::NAN;
    }
    unsafe { (*handle).inner.luminosity.value() }
}

// ═══════════════════════════════════════════════════════════════════════════
// Planet constants
// ═══════════════════════════════════════════════════════════════════════════

/// Get Mercury's orbital and physical parameters.
#[no_mangle]
pub extern "C" fn siderust_planet_mercury(out: *mut SiderustPlanet) -> SiderustStatus {
    if out.is_null() {
        return SiderustStatus::NullPointer;
    }
    unsafe { *out = SiderustPlanet::from_rust(&bodies::MERCURY) };
    SiderustStatus::Ok
}

/// Get Venus's orbital and physical parameters.
#[no_mangle]
pub extern "C" fn siderust_planet_venus(out: *mut SiderustPlanet) -> SiderustStatus {
    if out.is_null() {
        return SiderustStatus::NullPointer;
    }
    unsafe { *out = SiderustPlanet::from_rust(&bodies::VENUS) };
    SiderustStatus::Ok
}

/// Get Earth's orbital and physical parameters.
#[no_mangle]
pub extern "C" fn siderust_planet_earth(out: *mut SiderustPlanet) -> SiderustStatus {
    if out.is_null() {
        return SiderustStatus::NullPointer;
    }
    unsafe { *out = SiderustPlanet::from_rust(&bodies::EARTH) };
    SiderustStatus::Ok
}

/// Get Mars's orbital and physical parameters.
#[no_mangle]
pub extern "C" fn siderust_planet_mars(out: *mut SiderustPlanet) -> SiderustStatus {
    if out.is_null() {
        return SiderustStatus::NullPointer;
    }
    unsafe { *out = SiderustPlanet::from_rust(&bodies::MARS) };
    SiderustStatus::Ok
}

/// Get Jupiter's orbital and physical parameters.
#[no_mangle]
pub extern "C" fn siderust_planet_jupiter(out: *mut SiderustPlanet) -> SiderustStatus {
    if out.is_null() {
        return SiderustStatus::NullPointer;
    }
    unsafe { *out = SiderustPlanet::from_rust(&bodies::JUPITER) };
    SiderustStatus::Ok
}

/// Get Saturn's orbital and physical parameters.
#[no_mangle]
pub extern "C" fn siderust_planet_saturn(out: *mut SiderustPlanet) -> SiderustStatus {
    if out.is_null() {
        return SiderustStatus::NullPointer;
    }
    unsafe { *out = SiderustPlanet::from_rust(&bodies::SATURN) };
    SiderustStatus::Ok
}

/// Get Uranus's orbital and physical parameters.
#[no_mangle]
pub extern "C" fn siderust_planet_uranus(out: *mut SiderustPlanet) -> SiderustStatus {
    if out.is_null() {
        return SiderustStatus::NullPointer;
    }
    unsafe { *out = SiderustPlanet::from_rust(&bodies::URANUS) };
    SiderustStatus::Ok
}

/// Get Neptune's orbital and physical parameters.
#[no_mangle]
pub extern "C" fn siderust_planet_neptune(out: *mut SiderustPlanet) -> SiderustStatus {
    if out.is_null() {
        return SiderustStatus::NullPointer;
    }
    unsafe { *out = SiderustPlanet::from_rust(&bodies::NEPTUNE) };
    SiderustStatus::Ok
}
