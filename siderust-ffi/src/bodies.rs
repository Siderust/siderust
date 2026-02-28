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
    ffi_guard! {{
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

    }}
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
    ffi_guard! {{
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

    }}
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
    ffi_guard! {{
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

    }}
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
    ffi_guard! {{
        if out.is_null() {
            return SiderustStatus::NullPointer;
        }
        unsafe { *out = SiderustPlanet::from_rust(&bodies::MERCURY) };
        SiderustStatus::Ok

    }}
}

/// Get Venus's orbital and physical parameters.
#[no_mangle]
pub extern "C" fn siderust_planet_venus(out: *mut SiderustPlanet) -> SiderustStatus {
    ffi_guard! {{
        if out.is_null() {
            return SiderustStatus::NullPointer;
        }
        unsafe { *out = SiderustPlanet::from_rust(&bodies::VENUS) };
        SiderustStatus::Ok

    }}
}

/// Get Earth's orbital and physical parameters.
#[no_mangle]
pub extern "C" fn siderust_planet_earth(out: *mut SiderustPlanet) -> SiderustStatus {
    ffi_guard! {{
        if out.is_null() {
            return SiderustStatus::NullPointer;
        }
        unsafe { *out = SiderustPlanet::from_rust(&bodies::EARTH) };
        SiderustStatus::Ok

    }}
}

/// Get Mars's orbital and physical parameters.
#[no_mangle]
pub extern "C" fn siderust_planet_mars(out: *mut SiderustPlanet) -> SiderustStatus {
    ffi_guard! {{
        if out.is_null() {
            return SiderustStatus::NullPointer;
        }
        unsafe { *out = SiderustPlanet::from_rust(&bodies::MARS) };
        SiderustStatus::Ok

    }}
}

/// Get Jupiter's orbital and physical parameters.
#[no_mangle]
pub extern "C" fn siderust_planet_jupiter(out: *mut SiderustPlanet) -> SiderustStatus {
    ffi_guard! {{
        if out.is_null() {
            return SiderustStatus::NullPointer;
        }
        unsafe { *out = SiderustPlanet::from_rust(&bodies::JUPITER) };
        SiderustStatus::Ok

    }}
}

/// Get Saturn's orbital and physical parameters.
#[no_mangle]
pub extern "C" fn siderust_planet_saturn(out: *mut SiderustPlanet) -> SiderustStatus {
    ffi_guard! {{
        if out.is_null() {
            return SiderustStatus::NullPointer;
        }
        unsafe { *out = SiderustPlanet::from_rust(&bodies::SATURN) };
        SiderustStatus::Ok

    }}
}

/// Get Uranus's orbital and physical parameters.
#[no_mangle]
pub extern "C" fn siderust_planet_uranus(out: *mut SiderustPlanet) -> SiderustStatus {
    ffi_guard! {{
        if out.is_null() {
            return SiderustStatus::NullPointer;
        }
        unsafe { *out = SiderustPlanet::from_rust(&bodies::URANUS) };
        SiderustStatus::Ok

    }}
}

/// Get Neptune's orbital and physical parameters.
#[no_mangle]
pub extern "C" fn siderust_planet_neptune(out: *mut SiderustPlanet) -> SiderustStatus {
    ffi_guard! {{
        if out.is_null() {
            return SiderustStatus::NullPointer;
        }
        unsafe { *out = SiderustPlanet::from_rust(&bodies::NEPTUNE) };
        SiderustStatus::Ok

    }}
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::SiderustOrbit;
    use std::ffi::CString;
    use std::ptr;

    fn uninit_planet() -> SiderustPlanet {
        SiderustPlanet {
            mass_kg: 0.0,
            radius_km: 0.0,
            orbit: SiderustOrbit {
                semi_major_axis_au: 0.0,
                eccentricity: 0.0,
                inclination_deg: 0.0,
                lon_ascending_node_deg: 0.0,
                arg_perihelion_deg: 0.0,
                mean_anomaly_deg: 0.0,
                epoch_jd: 0.0,
            },
        }
    }

    // ── Planets ──────────────────────────────────────────────────────────

    macro_rules! planet_test {
        ($fn_name:ident, $ffi_fn:ident) => {
            #[test]
            fn $fn_name() {
                let mut out = uninit_planet();
                assert_eq!($ffi_fn(&mut out), SiderustStatus::Ok);
                assert!(out.mass_kg > 0.0);
                assert!(out.radius_km > 0.0);
            }
        };
    }

    planet_test!(mercury_ok, siderust_planet_mercury);
    planet_test!(venus_ok, siderust_planet_venus);
    planet_test!(earth_ok, siderust_planet_earth);
    planet_test!(mars_ok, siderust_planet_mars);
    planet_test!(jupiter_ok, siderust_planet_jupiter);
    planet_test!(saturn_ok, siderust_planet_saturn);
    planet_test!(uranus_ok, siderust_planet_uranus);
    planet_test!(neptune_ok, siderust_planet_neptune);

    #[test]
    fn planet_null_ptr() {
        assert_eq!(
            siderust_planet_earth(ptr::null_mut()),
            SiderustStatus::NullPointer
        );
        assert_eq!(
            siderust_planet_mars(ptr::null_mut()),
            SiderustStatus::NullPointer
        );
        assert_eq!(
            siderust_planet_jupiter(ptr::null_mut()),
            SiderustStatus::NullPointer
        );
        assert_eq!(
            siderust_planet_saturn(ptr::null_mut()),
            SiderustStatus::NullPointer
        );
        assert_eq!(
            siderust_planet_uranus(ptr::null_mut()),
            SiderustStatus::NullPointer
        );
        assert_eq!(
            siderust_planet_neptune(ptr::null_mut()),
            SiderustStatus::NullPointer
        );
        assert_eq!(
            siderust_planet_mercury(ptr::null_mut()),
            SiderustStatus::NullPointer
        );
        assert_eq!(
            siderust_planet_venus(ptr::null_mut()),
            SiderustStatus::NullPointer
        );
    }

    // ── Star catalog ─────────────────────────────────────────────────────

    fn catalog_star(name: &str) -> *mut SiderustStar {
        let cname = CString::new(name).unwrap();
        let mut handle: *mut SiderustStar = ptr::null_mut();
        let s = unsafe { siderust_star_catalog(cname.as_ptr(), &mut handle) };
        assert_eq!(s, SiderustStatus::Ok, "Catalog lookup for {name} failed");
        assert!(!handle.is_null());
        handle
    }

    #[test]
    fn catalog_all_known_stars() {
        for name in &[
            "VEGA",
            "SIRIUS",
            "POLARIS",
            "CANOPUS",
            "ARCTURUS",
            "RIGEL",
            "BETELGEUSE",
            "PROCYON",
            "ALDEBARAN",
            "ALTAIR",
        ] {
            let h = catalog_star(name);
            // Check distance is finite and positive
            let dist = unsafe { siderust_star_distance_ly(h) };
            assert!(dist > 0.0 && dist.is_finite(), "{name}: distance {dist}");
            unsafe { siderust_star_free(h) };
        }
    }

    #[test]
    fn catalog_lowercase_name() {
        let cname = CString::new("vega").unwrap();
        let mut handle: *mut SiderustStar = ptr::null_mut();
        let s = unsafe { siderust_star_catalog(cname.as_ptr(), &mut handle) };
        assert_eq!(s, SiderustStatus::Ok);
        unsafe { siderust_star_free(handle) };
    }

    #[test]
    fn catalog_unknown_star() {
        let cname = CString::new("NOTASTAR").unwrap();
        let mut handle: *mut SiderustStar = ptr::null_mut();
        assert_eq!(
            unsafe { siderust_star_catalog(cname.as_ptr(), &mut handle) },
            SiderustStatus::UnknownStar
        );
        assert!(handle.is_null());
    }

    #[test]
    fn catalog_null_name() {
        let mut handle: *mut SiderustStar = ptr::null_mut();
        assert_eq!(
            unsafe { siderust_star_catalog(ptr::null(), &mut handle) },
            SiderustStatus::NullPointer
        );
    }

    #[test]
    fn catalog_null_out() {
        let cname = CString::new("VEGA").unwrap();
        assert_eq!(
            unsafe { siderust_star_catalog(cname.as_ptr(), ptr::null_mut()) },
            SiderustStatus::NullPointer
        );
    }

    // ── Star properties ───────────────────────────────────────────────────

    #[test]
    fn star_null_handle_returns_nan() {
        let dist = siderust_star_distance_ly(ptr::null());
        assert!(dist.is_nan());
        let mass = siderust_star_mass_solar(ptr::null());
        assert!(mass.is_nan());
        let rad = siderust_star_radius_solar(ptr::null());
        assert!(rad.is_nan());
        let lum = siderust_star_luminosity_solar(ptr::null());
        assert!(lum.is_nan());
    }

    #[test]
    fn star_name_roundtrip() {
        let h = catalog_star("VEGA");
        let mut buf = vec![0i8; 32];
        let mut written = 0usize;
        let s = unsafe { siderust_star_name(h, buf.as_mut_ptr(), buf.len(), &mut written) };
        assert_eq!(s, SiderustStatus::Ok);
        assert!(written > 0);
        unsafe { siderust_star_free(h) };
    }

    #[test]
    fn star_name_null_ptr() {
        let h = catalog_star("SIRIUS");
        let s = unsafe { siderust_star_name(h, ptr::null_mut(), 32, ptr::null_mut()) };
        assert_eq!(s, SiderustStatus::NullPointer);
        // null handle
        let mut buf = vec![0i8; 32];
        let s2 = unsafe { siderust_star_name(ptr::null(), buf.as_mut_ptr(), 32, ptr::null_mut()) };
        assert_eq!(s2, SiderustStatus::NullPointer);
        unsafe { siderust_star_free(h) };
    }

    #[test]
    fn star_name_buffer_too_small() {
        let h = catalog_star("BETELGEUSE"); // 10 chars + NUL = 11 bytes
        let mut buf = vec![0i8; 4]; // too small
        let mut written = 0usize;
        let s = unsafe { siderust_star_name(h, buf.as_mut_ptr(), 4, &mut written) };
        assert_eq!(s, SiderustStatus::InvalidArgument);
        unsafe { siderust_star_free(h) };
    }

    // ── Custom star creation ──────────────────────────────────────────────

    #[test]
    fn custom_star_no_proper_motion() {
        let cname = CString::new("TestStar").unwrap();
        let mut handle: *mut SiderustStar = ptr::null_mut();
        let s = unsafe {
            siderust_star_create(
                cname.as_ptr(),
                10.0, // distance_ly
                1.0,  // mass_solar
                1.0,  // radius_solar
                1.0,  // luminosity_solar
                90.0, // ra_deg
                45.0, // dec_deg
                2_451_545.0,
                ptr::null(),
                &mut handle,
            )
        };
        assert_eq!(s, SiderustStatus::Ok);
        assert!(!handle.is_null());
        assert!((siderust_star_distance_ly(handle) - 10.0).abs() < 1e-9);
        assert!((siderust_star_mass_solar(handle) - 1.0).abs() < 1e-9);
        assert!((siderust_star_radius_solar(handle) - 1.0).abs() < 1e-9);
        assert!((siderust_star_luminosity_solar(handle) - 1.0).abs() < 1e-9);
        unsafe { siderust_star_free(handle) };
    }

    #[test]
    fn custom_star_with_proper_motion() {
        let cname = CString::new("PMstar").unwrap();
        let pm = SiderustProperMotion {
            pm_ra_deg_yr: 0.01,
            pm_dec_deg_yr: -0.005,
            ra_convention: SiderustRaConvention::MuAlpha,
        };
        let mut handle: *mut SiderustStar = ptr::null_mut();
        let s = unsafe {
            siderust_star_create(
                cname.as_ptr(),
                5.0,
                0.5,
                0.5,
                0.1,
                100.0,
                30.0,
                2_451_545.0,
                &pm as *const _,
                &mut handle,
            )
        };
        assert_eq!(s, SiderustStatus::Ok);
        unsafe { siderust_star_free(handle) };
    }

    #[test]
    fn custom_star_mu_alpha_star_convention() {
        let cname = CString::new("PMstar2").unwrap();
        let pm = SiderustProperMotion {
            pm_ra_deg_yr: 0.01,
            pm_dec_deg_yr: -0.005,
            ra_convention: SiderustRaConvention::MuAlphaStar,
        };
        let mut handle: *mut SiderustStar = ptr::null_mut();
        let s = unsafe {
            siderust_star_create(
                cname.as_ptr(),
                5.0,
                0.5,
                0.5,
                0.1,
                100.0,
                30.0,
                2_451_545.0,
                &pm as *const _,
                &mut handle,
            )
        };
        assert_eq!(s, SiderustStatus::Ok);
        unsafe { siderust_star_free(handle) };
    }

    #[test]
    fn custom_star_null_name_returns_null_pointer() {
        let mut handle: *mut SiderustStar = ptr::null_mut();
        let s = unsafe {
            siderust_star_create(
                ptr::null(),
                5.0,
                1.0,
                1.0,
                1.0,
                0.0,
                0.0,
                2_451_545.0,
                ptr::null(),
                &mut handle,
            )
        };
        assert_eq!(s, SiderustStatus::NullPointer);
    }

    #[test]
    fn custom_star_null_out_returns_null_pointer() {
        let cname = CString::new("X").unwrap();
        let s = unsafe {
            siderust_star_create(
                cname.as_ptr(),
                5.0,
                1.0,
                1.0,
                1.0,
                0.0,
                0.0,
                2_451_545.0,
                ptr::null(),
                ptr::null_mut(),
            )
        };
        assert_eq!(s, SiderustStatus::NullPointer);
    }

    #[test]
    fn star_free_null_is_safe() {
        // Should not crash
        unsafe { siderust_star_free(ptr::null_mut()) };
    }
}
