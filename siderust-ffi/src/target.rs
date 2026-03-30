// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Target FFI, opaque handles for celestial targets.
//!
//! This module exposes the forward-looking target handle for the C ABI:
//! [`SiderustGenericTarget`], which mirrors Rust's `CoordinateWithPM<T>`.

use crate::error::SiderustStatus;
use crate::types::*;
use qtty::*;
use siderust::astro::proper_motion::ProperMotion;
use siderust::coordinates::spherical;
use siderust::targets::CoordinateWithPM;
use siderust::time::JulianDate;
use tempoch_ffi::TempochJd;

/// Unit type alias for degrees per Julian year.
type DegreePerYear = qtty::Per<Degree, Year>;
/// Quantity alias for degrees per Julian year.
type DegreesPerYear = qtty::Quantity<DegreePerYear>;

// ═══════════════════════════════════════════════════════════════════════════
// Generic Target Handle (CoordinateWithPM<T>)
// ═══════════════════════════════════════════════════════════════════════════

/// Opaque handle representing a generic celestial target.
///
/// This mirrors Rust's `CoordinateWithPM<T>`, supporting:
/// - Direction-only coordinates (ICRS, equatorial, etc.)
/// - Full position coordinates with distance
/// - Optional proper motion
///
/// Use `siderust_generic_target_create_*` functions to construct,
/// and `siderust_generic_target_free` to release.
pub struct SiderustGenericTarget {
    /// The underlying coordinate-with-proper-motion.
    /// For FFI, we store the ICRS direction variant since that's what
    /// altitude/azimuth computations ultimately require.
    pub(crate) inner: CoordinateWithPM<spherical::direction::ICRS>,
    /// Original coordinate data for property access.
    pub(crate) data: SiderustGenericTargetData,
}

/// Create a generic target from an ICRS spherical direction.
///
/// # Parameters
/// - `ra_deg`: Right ascension in degrees.
/// - `dec_deg`: Declination in degrees.
/// - `epoch_jd`: Reference epoch as Julian Date.
/// - `out`: On success, receives the allocated handle pointer.
///
/// The caller must free the handle with [`siderust_generic_target_free`].
#[no_mangle]
pub extern "C" fn siderust_generic_target_create_icrs(
    ra_deg: f64,
    dec_deg: f64,
    epoch_jd: f64,
    out: *mut *mut SiderustGenericTarget,
) -> SiderustStatus {
    ffi_guard! {{
        if out.is_null() {
            return SiderustStatus::NullPointer;
        }
        let dir = spherical::direction::ICRS::new(Degrees::new(ra_deg), Degrees::new(dec_deg));
        let coord = CoordinateWithPM::new_static(dir, JulianDate::new(epoch_jd));
        let data = SiderustGenericTargetData {
            kind: SiderustTargetCoordKind::SphericalDir,
            _pad1: [0; 4],
            coord: SiderustTargetCoordUnion {
                spherical_dir: SiderustSphericalDir {
                    polar_deg: dec_deg,
                    azimuth_deg: ra_deg,
                    frame: SiderustFrame::ICRS,
                },
            },
            epoch_jd: TempochJd::new(epoch_jd),
            has_proper_motion: false,
            _pad2: [0; 7],
            proper_motion: SiderustProperMotion {
                pm_ra_deg_yr: 0.0,
                pm_dec_deg_yr: 0.0,
                ra_convention: SiderustRaConvention::MuAlpha,
            },
        };
        let handle = Box::new(SiderustGenericTarget { inner: coord, data });
        unsafe { *out = Box::into_raw(handle) };
        SiderustStatus::Ok
    }}
}

/// Create a generic target from an ICRS direction with proper motion.
///
/// # Parameters
/// - `ra_deg`: Right ascension in degrees.
/// - `dec_deg`: Declination in degrees.
/// - `epoch_jd`: Reference epoch as Julian Date.
/// - `pm_ra_deg_yr`: RA proper motion in degrees per Julian year.
/// - `pm_dec_deg_yr`: Dec proper motion in degrees per Julian year.
/// - `ra_convention`: Whether `pm_ra` is µα or µα★.
/// - `out`: On success, receives the allocated handle pointer.
///
/// The caller must free the handle with [`siderust_generic_target_free`].
#[no_mangle]
pub extern "C" fn siderust_generic_target_create_icrs_with_pm(
    ra_deg: f64,
    dec_deg: f64,
    epoch_jd: f64,
    pm_ra_deg_yr: f64,
    pm_dec_deg_yr: f64,
    ra_convention: SiderustRaConvention,
    out: *mut *mut SiderustGenericTarget,
) -> SiderustStatus {
    ffi_guard! {{
        if out.is_null() {
            return SiderustStatus::NullPointer;
        }
        let dir = spherical::direction::ICRS::new(Degrees::new(ra_deg), Degrees::new(dec_deg));

        // Convert FFI proper motion to Rust
        let pm = match ra_convention {
            SiderustRaConvention::MuAlpha => ProperMotion::from_mu_alpha::<DegreePerYear>(
                DegreesPerYear::new(pm_ra_deg_yr),
                DegreesPerYear::new(pm_dec_deg_yr),
            ),
            SiderustRaConvention::MuAlphaStar => ProperMotion::from_mu_alpha_star::<DegreePerYear>(
                DegreesPerYear::new(pm_ra_deg_yr),
                DegreesPerYear::new(pm_dec_deg_yr),
            ),
        };

        let coord = CoordinateWithPM::new(dir, JulianDate::new(epoch_jd), pm);
        let data = SiderustGenericTargetData {
            kind: SiderustTargetCoordKind::SphericalDir,
            _pad1: [0; 4],
            coord: SiderustTargetCoordUnion {
                spherical_dir: SiderustSphericalDir {
                    polar_deg: dec_deg,
                    azimuth_deg: ra_deg,
                    frame: SiderustFrame::ICRS,
                },
            },
            epoch_jd: TempochJd::new(epoch_jd),
            has_proper_motion: true,
            _pad2: [0; 7],
            proper_motion: SiderustProperMotion {
                pm_ra_deg_yr,
                pm_dec_deg_yr,
                ra_convention,
            },
        };
        let handle = Box::new(SiderustGenericTarget { inner: coord, data });
        unsafe { *out = Box::into_raw(handle) };
        SiderustStatus::Ok
    }}
}

/// Free a generic target handle.
///
/// # Safety
/// - `handle` must have been allocated by a `siderust_generic_target_create_*` function.
/// - The handle must not have been freed before, and must not be used after this call.
#[no_mangle]
pub unsafe extern "C" fn siderust_generic_target_free(handle: *mut SiderustGenericTarget) {
    if !handle.is_null() {
        drop(unsafe { Box::from_raw(handle) });
    }
}

/// Get the coordinate data from a generic target.
///
/// Writes the target's coordinate information to `*out`.
#[no_mangle]
pub extern "C" fn siderust_generic_target_get_data(
    handle: *const SiderustGenericTarget,
    out: *mut SiderustGenericTargetData,
) -> SiderustStatus {
    ffi_guard! {{
        if handle.is_null() || out.is_null() {
            return SiderustStatus::NullPointer;
        }
        unsafe { *out = (*handle).data };
        SiderustStatus::Ok
    }}
}

/// Get the epoch (Julian Date) of a generic target.
#[no_mangle]
pub extern "C" fn siderust_generic_target_epoch_jd(
    handle: *const SiderustGenericTarget,
    out: *mut f64,
) -> SiderustStatus {
    ffi_guard! {{
        if handle.is_null() || out.is_null() {
            return SiderustStatus::NullPointer;
        }
        unsafe { *out = (*handle).data.epoch_jd.value };
        SiderustStatus::Ok
    }}
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::ptr;

    fn create_icrs_target() -> *mut SiderustGenericTarget {
        let mut h: *mut SiderustGenericTarget = ptr::null_mut();
        assert_eq!(
            siderust_generic_target_create_icrs(279.23, 38.78, 2451545.0, &mut h),
            SiderustStatus::Ok
        );
        assert!(!h.is_null());
        h
    }

    #[test]
    fn create_null_out() {
        assert_eq!(
            siderust_generic_target_create_icrs(0.0, 0.0, 0.0, ptr::null_mut()),
            SiderustStatus::NullPointer
        );
    }

    #[test]
    fn create_and_free() {
        let h = create_icrs_target();
        unsafe { siderust_generic_target_free(h) };
    }

    #[test]
    fn free_null_safe() {
        unsafe { siderust_generic_target_free(ptr::null_mut()) };
    }

    #[test]
    fn get_data_roundtrip() {
        let h = create_icrs_target();
        let mut out = SiderustGenericTargetData {
            kind: SiderustTargetCoordKind::CartesianPos,
            _pad1: [0; 4],
            coord: SiderustTargetCoordUnion {
                cartesian_pos: SiderustCartesianPos {
                    x: 0.0,
                    y: 0.0,
                    z: 0.0,
                    frame: SiderustFrame::ICRS,
                    center: SiderustCenter::Geocentric,
                    length_unit: SiderustLengthUnit::Meter,
                },
            },
            epoch_jd: TempochJd::new(0.0),
            has_proper_motion: true,
            _pad2: [0; 7],
            proper_motion: SiderustProperMotion {
                pm_ra_deg_yr: 1.0,
                pm_dec_deg_yr: 1.0,
                ra_convention: SiderustRaConvention::MuAlphaStar,
            },
        };
        assert_eq!(
            siderust_generic_target_get_data(h, &mut out),
            SiderustStatus::Ok
        );
        assert_eq!(out.kind, SiderustTargetCoordKind::SphericalDir);
        let coord = unsafe { out.coord.spherical_dir };
        assert_eq!(coord.frame, SiderustFrame::ICRS);
        assert!((coord.azimuth_deg - 279.23).abs() < 1e-10);
        assert!((coord.polar_deg - 38.78).abs() < 1e-10);
        assert_eq!(out.epoch_jd.value, 2451545.0);
        assert!(!out.has_proper_motion);
        unsafe { siderust_generic_target_free(h) };
    }

    #[test]
    fn create_with_pm_sets_payload() {
        let mut h: *mut SiderustGenericTarget = ptr::null_mut();
        assert_eq!(
            siderust_generic_target_create_icrs_with_pm(
                279.23,
                38.78,
                2451545.0,
                0.0001,
                -0.0002,
                SiderustRaConvention::MuAlphaStar,
                &mut h,
            ),
            SiderustStatus::Ok
        );
        let mut out = SiderustGenericTargetData {
            kind: SiderustTargetCoordKind::CartesianPos,
            _pad1: [0; 4],
            coord: SiderustTargetCoordUnion {
                cartesian_pos: SiderustCartesianPos {
                    x: 0.0,
                    y: 0.0,
                    z: 0.0,
                    frame: SiderustFrame::ICRS,
                    center: SiderustCenter::Geocentric,
                    length_unit: SiderustLengthUnit::Meter,
                },
            },
            epoch_jd: TempochJd::new(0.0),
            has_proper_motion: false,
            _pad2: [0; 7],
            proper_motion: SiderustProperMotion {
                pm_ra_deg_yr: 0.0,
                pm_dec_deg_yr: 0.0,
                ra_convention: SiderustRaConvention::MuAlpha,
            },
        };
        assert_eq!(
            siderust_generic_target_get_data(h, &mut out),
            SiderustStatus::Ok
        );
        assert!(out.has_proper_motion);
        assert!((out.proper_motion.pm_ra_deg_yr - 0.0001).abs() < 1e-12);
        assert!((out.proper_motion.pm_dec_deg_yr + 0.0002).abs() < 1e-12);
        assert_eq!(
            out.proper_motion.ra_convention,
            SiderustRaConvention::MuAlphaStar
        );
        unsafe { siderust_generic_target_free(h) };
    }

    #[test]
    fn epoch_jd_roundtrip() {
        let h = create_icrs_target();
        let mut out = 0.0f64;
        assert_eq!(
            siderust_generic_target_epoch_jd(h, &mut out),
            SiderustStatus::Ok
        );
        assert!((out - 2451545.0).abs() < 1e-6);
        unsafe { siderust_generic_target_free(h) };
    }

    #[test]
    fn get_data_null_handle() {
        let mut out = SiderustGenericTargetData {
            kind: SiderustTargetCoordKind::CartesianPos,
            _pad1: [0; 4],
            coord: SiderustTargetCoordUnion {
                cartesian_pos: SiderustCartesianPos {
                    x: 0.0,
                    y: 0.0,
                    z: 0.0,
                    frame: SiderustFrame::ICRS,
                    center: SiderustCenter::Geocentric,
                    length_unit: SiderustLengthUnit::Meter,
                },
            },
            epoch_jd: TempochJd::new(0.0),
            has_proper_motion: false,
            _pad2: [0; 7],
            proper_motion: SiderustProperMotion {
                pm_ra_deg_yr: 0.0,
                pm_dec_deg_yr: 0.0,
                ra_convention: SiderustRaConvention::MuAlpha,
            },
        };
        assert_eq!(
            siderust_generic_target_get_data(ptr::null(), &mut out),
            SiderustStatus::NullPointer
        );
    }

    #[test]
    fn get_data_null_out() {
        let h = create_icrs_target();
        assert_eq!(
            siderust_generic_target_get_data(h, ptr::null_mut()),
            SiderustStatus::NullPointer
        );
        unsafe { siderust_generic_target_free(h) };
    }

    #[test]
    fn epoch_null_out() {
        let h = create_icrs_target();
        assert_eq!(
            siderust_generic_target_epoch_jd(h, ptr::null_mut()),
            SiderustStatus::NullPointer
        );
        unsafe { siderust_generic_target_free(h) };
    }
}
