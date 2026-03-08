// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Target FFI — opaque handles for celestial targets.
//!
//! This module provides two target abstractions:
//!
//! - [`SiderustTarget`]: Legacy ICRS direction-only target (RA/Dec).
//! - [`SiderustGenericTarget`]: Generic target supporting direction or position
//!   coordinates with optional proper motion, mirroring Rust's `CoordinateWithPM<T>`.
//!
//! Both expose altitude and azimuth queries by delegating to the unified
//! [`crate::subject`] module.

use crate::error::SiderustStatus;
use crate::types::*;
use qtty::*;
use siderust::astro::proper_motion::ProperMotion;
use siderust::coordinates::spherical;
use siderust::targets::CoordinateWithPM;
use siderust::time::JulianDate;

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
            epoch_jd,
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
            epoch_jd,
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
        unsafe { *out = (*handle).data.epoch_jd };
        SiderustStatus::Ok
    }}
}

// ═══════════════════════════════════════════════════════════════════════════
// Legacy ICRS Target Handle
// ═══════════════════════════════════════════════════════════════════════════

/// Opaque handle representing an ICRS celestial target (RA/Dec).
///
/// This is the legacy target type. New code should prefer [`SiderustGenericTarget`].
pub struct SiderustTarget {
    pub(crate) dir: spherical::direction::ICRS,
    pub(crate) ra_deg: f64,
    pub(crate) dec_deg: f64,
    pub(crate) epoch_jd: f64,
}

// ═══════════════════════════════════════════════════════════════════════════
// Lifecycle
// ═══════════════════════════════════════════════════════════════════════════

/// Create a new target from right ascension and declination (degrees) and epoch
/// (Julian Date).
///
/// On success, writes a newly-allocated pointer to `*out`.
/// The caller must free it with [`siderust_target_free`].
#[no_mangle]
pub extern "C" fn siderust_target_create(
    ra_deg: f64,
    dec_deg: f64,
    epoch_jd: f64,
    out: *mut *mut SiderustTarget,
) -> SiderustStatus {
    ffi_guard! {{
        if out.is_null() {
            return SiderustStatus::NullPointer;
        }
        let dir = spherical::direction::ICRS::new(Degrees::new(ra_deg), Degrees::new(dec_deg));
        let handle = Box::new(SiderustTarget {
            dir,
            ra_deg,
            dec_deg,
            epoch_jd,
        });
        unsafe { *out = Box::into_raw(handle) };
        SiderustStatus::Ok
    }}
}

/// Free a target handle created by [`siderust_target_create`].
///
/// # Safety
///
/// * `handle` must have been allocated by [`siderust_target_create`].
/// * The handle must not have been freed before, and must not be used after
///   this call.
#[no_mangle]
pub unsafe extern "C" fn siderust_target_free(handle: *mut SiderustTarget) {
    if !handle.is_null() {
        // SAFETY: caller guarantees the handle was allocated by
        // `siderust_target_create` and has not been freed before.
        drop(unsafe { Box::from_raw(handle) });
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// Property accessors
// ═══════════════════════════════════════════════════════════════════════════

/// Write the right ascension (degrees) to `*out`.
#[no_mangle]
pub extern "C" fn siderust_target_ra_deg(
    handle: *const SiderustTarget,
    out: *mut f64,
) -> SiderustStatus {
    ffi_guard! {{
        if handle.is_null() || out.is_null() {
            return SiderustStatus::NullPointer;
        }
        unsafe { *out = (*handle).ra_deg };
        SiderustStatus::Ok
    }}
}

/// Write the declination (degrees) to `*out`.
#[no_mangle]
pub extern "C" fn siderust_target_dec_deg(
    handle: *const SiderustTarget,
    out: *mut f64,
) -> SiderustStatus {
    ffi_guard! {{
        if handle.is_null() || out.is_null() {
            return SiderustStatus::NullPointer;
        }
        unsafe { *out = (*handle).dec_deg };
        SiderustStatus::Ok
    }}
}

/// Write the reference epoch (Julian Date) to `*out`.
#[no_mangle]
pub extern "C" fn siderust_target_epoch_jd(
    handle: *const SiderustTarget,
    out: *mut f64,
) -> SiderustStatus {
    ffi_guard! {{
        if handle.is_null() || out.is_null() {
            return SiderustStatus::NullPointer;
        }
        unsafe { *out = (*handle).epoch_jd };
        SiderustStatus::Ok
    }}
}

// ═══════════════════════════════════════════════════════════════════════════
// Altitude  — thin wrappers
// ═══════════════════════════════════════════════════════════════════════════

/// Altitude of the target (radians).
#[no_mangle]
pub extern "C" fn siderust_target_altitude_at(
    handle: *const SiderustTarget,
    observer: SiderustGeodetict,
    mjd: f64,
    out_rad: *mut f64,
) -> SiderustStatus {
    crate::subject::siderust_altitude_at(SiderustSubject::target(handle), observer, mjd, out_rad)
}

/// Periods when the target is above `threshold_deg`.
#[no_mangle]
pub extern "C" fn siderust_target_above_threshold(
    handle: *const SiderustTarget,
    observer: SiderustGeodetict,
    window: TempochPeriodMjd,
    threshold_deg: f64,
    opts: SiderustSearchOpts,
    out: *mut *mut TempochPeriodMjd,
    count: *mut usize,
) -> SiderustStatus {
    crate::subject::siderust_above_threshold(
        SiderustSubject::target(handle),
        observer,
        window,
        threshold_deg,
        opts,
        out,
        count,
    )
}

/// Periods when the target is below `threshold_deg`.
#[no_mangle]
pub extern "C" fn siderust_target_below_threshold(
    handle: *const SiderustTarget,
    observer: SiderustGeodetict,
    window: TempochPeriodMjd,
    threshold_deg: f64,
    opts: SiderustSearchOpts,
    out: *mut *mut TempochPeriodMjd,
    count: *mut usize,
) -> SiderustStatus {
    crate::subject::siderust_below_threshold(
        SiderustSubject::target(handle),
        observer,
        window,
        threshold_deg,
        opts,
        out,
        count,
    )
}

/// Altitude crossing events for the target.
#[no_mangle]
pub extern "C" fn siderust_target_crossings(
    handle: *const SiderustTarget,
    observer: SiderustGeodetict,
    window: TempochPeriodMjd,
    threshold_deg: f64,
    opts: SiderustSearchOpts,
    out: *mut *mut SiderustCrossingEvent,
    count: *mut usize,
) -> SiderustStatus {
    crate::subject::siderust_crossings(
        SiderustSubject::target(handle),
        observer,
        window,
        threshold_deg,
        opts,
        out,
        count,
    )
}

/// Culmination events for the target.
#[no_mangle]
pub extern "C" fn siderust_target_culminations(
    handle: *const SiderustTarget,
    observer: SiderustGeodetict,
    window: TempochPeriodMjd,
    opts: SiderustSearchOpts,
    out: *mut *mut SiderustCulminationEvent,
    count: *mut usize,
) -> SiderustStatus {
    crate::subject::siderust_culminations(
        SiderustSubject::target(handle),
        observer,
        window,
        opts,
        out,
        count,
    )
}

// ═══════════════════════════════════════════════════════════════════════════
// Azimuth  — thin wrappers
// ═══════════════════════════════════════════════════════════════════════════

/// Azimuth of the target (degrees, North-clockwise).
#[no_mangle]
pub extern "C" fn siderust_target_azimuth_at(
    handle: *const SiderustTarget,
    observer: SiderustGeodetict,
    mjd: f64,
    out_deg: *mut f64,
) -> SiderustStatus {
    crate::subject::siderust_azimuth_at(SiderustSubject::target(handle), observer, mjd, out_deg)
}

/// Azimuth bearing-crossing events for the target.
#[no_mangle]
pub extern "C" fn siderust_target_azimuth_crossings(
    handle: *const SiderustTarget,
    observer: SiderustGeodetict,
    window: TempochPeriodMjd,
    bearing_deg: f64,
    opts: SiderustSearchOpts,
    out: *mut *mut SiderustAzimuthCrossingEvent,
    count: *mut usize,
) -> SiderustStatus {
    crate::subject::siderust_azimuth_crossings(
        SiderustSubject::target(handle),
        observer,
        window,
        bearing_deg,
        opts,
        out,
        count,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_helpers::*;
    use std::ptr;

    // Vega: RA = 279.23°, Dec = 38.78° (approximate)
    fn create_vega() -> *mut SiderustTarget {
        let mut h: *mut SiderustTarget = ptr::null_mut();
        assert_eq!(
            siderust_target_create(279.23, 38.78, 2451545.0, &mut h),
            SiderustStatus::Ok
        );
        assert!(!h.is_null());
        h
    }

    // ── Lifecycle ─────────────────────────────────────────────────────────

    #[test]
    fn create_null_out() {
        assert_eq!(
            siderust_target_create(0.0, 0.0, 0.0, ptr::null_mut()),
            SiderustStatus::NullPointer
        );
    }

    #[test]
    fn create_and_free() {
        let h = create_vega();
        unsafe { siderust_target_free(h) };
    }

    #[test]
    fn free_null_safe() {
        unsafe { siderust_target_free(ptr::null_mut()) };
    }

    // ── Accessors ─────────────────────────────────────────────────────────

    #[test]
    fn ra_deg_roundtrip() {
        let h = create_vega();
        let mut out = 0.0f64;
        assert_eq!(siderust_target_ra_deg(h, &mut out), SiderustStatus::Ok);
        assert!((out - 279.23).abs() < 1e-10);
        unsafe { siderust_target_free(h) };
    }

    #[test]
    fn dec_deg_roundtrip() {
        let h = create_vega();
        let mut out = 0.0f64;
        assert_eq!(siderust_target_dec_deg(h, &mut out), SiderustStatus::Ok);
        assert!((out - 38.78).abs() < 1e-10);
        unsafe { siderust_target_free(h) };
    }

    #[test]
    fn epoch_jd_roundtrip() {
        let h = create_vega();
        let mut out = 0.0f64;
        assert_eq!(siderust_target_epoch_jd(h, &mut out), SiderustStatus::Ok);
        assert!((out - 2451545.0).abs() < 1e-6);
        unsafe { siderust_target_free(h) };
    }

    #[test]
    fn ra_null_handle() {
        let mut out = 0.0f64;
        assert_eq!(
            siderust_target_ra_deg(ptr::null(), &mut out),
            SiderustStatus::NullPointer
        );
    }

    #[test]
    fn ra_null_out() {
        let h = create_vega();
        assert_eq!(
            siderust_target_ra_deg(h, ptr::null_mut()),
            SiderustStatus::NullPointer
        );
        unsafe { siderust_target_free(h) };
    }

    #[test]
    fn dec_null_out() {
        let h = create_vega();
        assert_eq!(
            siderust_target_dec_deg(h, ptr::null_mut()),
            SiderustStatus::NullPointer
        );
        unsafe { siderust_target_free(h) };
    }

    #[test]
    fn epoch_null_out() {
        let h = create_vega();
        assert_eq!(
            siderust_target_epoch_jd(h, ptr::null_mut()),
            SiderustStatus::NullPointer
        );
        unsafe { siderust_target_free(h) };
    }

    // ── Altitude ──────────────────────────────────────────────────────────

    #[test]
    fn altitude_at_is_finite() {
        let h = create_vega();
        let mut out = 0.0f64;
        assert_eq!(
            siderust_target_altitude_at(h, paris(), 60000.0, &mut out),
            SiderustStatus::Ok
        );
        assert!(out.is_finite());
        unsafe { siderust_target_free(h) };
    }

    #[test]
    fn altitude_at_null_handle() {
        let mut out = 0.0f64;
        assert_eq!(
            siderust_target_altitude_at(ptr::null(), paris(), 60000.0, &mut out),
            SiderustStatus::NullPointer
        );
    }

    #[test]
    fn altitude_at_null_out() {
        let h = create_vega();
        assert_eq!(
            siderust_target_altitude_at(h, paris(), 60000.0, ptr::null_mut()),
            SiderustStatus::NullPointer
        );
        unsafe { siderust_target_free(h) };
    }

    #[test]
    fn above_threshold_ok() {
        let h = create_vega();
        let mut out: *mut TempochPeriodMjd = ptr::null_mut();
        let mut count = 0usize;
        let s = siderust_target_above_threshold(
            h,
            paris(),
            one_day_window(),
            0.0,
            default_opts(),
            &mut out,
            &mut count,
        );
        assert_eq!(s, SiderustStatus::Ok);
        unsafe {
            crate::altitude::siderust_periods_free(out, count);
            siderust_target_free(h);
        }
    }

    #[test]
    fn above_threshold_null_handle() {
        let mut out: *mut TempochPeriodMjd = ptr::null_mut();
        let mut count = 0usize;
        assert_eq!(
            siderust_target_above_threshold(
                ptr::null(),
                paris(),
                one_day_window(),
                0.0,
                default_opts(),
                &mut out,
                &mut count
            ),
            SiderustStatus::NullPointer
        );
    }

    #[test]
    fn above_threshold_invalid_window() {
        let h = create_vega();
        let bad = TempochPeriodMjd {
            start_mjd: 60001.0,
            end_mjd: 60000.0,
        };
        let mut out: *mut TempochPeriodMjd = ptr::null_mut();
        let mut count = 0usize;
        assert_eq!(
            siderust_target_above_threshold(
                h,
                paris(),
                bad,
                0.0,
                default_opts(),
                &mut out,
                &mut count
            ),
            SiderustStatus::InvalidPeriod
        );
        unsafe { siderust_target_free(h) };
    }

    #[test]
    fn crossings_ok() {
        let h = create_vega();
        let mut out: *mut SiderustCrossingEvent = ptr::null_mut();
        let mut count = 0usize;
        let s = siderust_target_crossings(
            h,
            paris(),
            one_day_window(),
            0.0,
            default_opts(),
            &mut out,
            &mut count,
        );
        assert_eq!(s, SiderustStatus::Ok);
        unsafe {
            crate::altitude::siderust_crossings_free(out, count);
            siderust_target_free(h);
        }
    }

    #[test]
    fn crossings_null_handle() {
        let mut out: *mut SiderustCrossingEvent = ptr::null_mut();
        let mut count = 0usize;
        assert_eq!(
            siderust_target_crossings(
                ptr::null(),
                paris(),
                one_day_window(),
                0.0,
                default_opts(),
                &mut out,
                &mut count
            ),
            SiderustStatus::NullPointer
        );
    }

    #[test]
    fn culminations_ok() {
        let h = create_vega();
        let mut out: *mut SiderustCulminationEvent = ptr::null_mut();
        let mut count = 0usize;
        let s = siderust_target_culminations(
            h,
            paris(),
            one_day_window(),
            default_opts(),
            &mut out,
            &mut count,
        );
        assert_eq!(s, SiderustStatus::Ok);
        unsafe {
            crate::altitude::siderust_culminations_free(out, count);
            siderust_target_free(h);
        }
    }

    #[test]
    fn culminations_null_handle() {
        let mut out: *mut SiderustCulminationEvent = ptr::null_mut();
        let mut count = 0usize;
        assert_eq!(
            siderust_target_culminations(
                ptr::null(),
                paris(),
                one_day_window(),
                default_opts(),
                &mut out,
                &mut count
            ),
            SiderustStatus::NullPointer
        );
    }

    // ── Azimuth ───────────────────────────────────────────────────────────

    #[test]
    fn azimuth_at_is_finite() {
        let h = create_vega();
        let mut out = 0.0f64;
        assert_eq!(
            siderust_target_azimuth_at(h, paris(), 60000.0, &mut out),
            SiderustStatus::Ok
        );
        assert!(out.is_finite() && out >= 0.0 && out < 360.0);
        unsafe { siderust_target_free(h) };
    }

    #[test]
    fn azimuth_at_null_handle() {
        let mut out = 0.0f64;
        assert_eq!(
            siderust_target_azimuth_at(ptr::null(), paris(), 60000.0, &mut out),
            SiderustStatus::NullPointer
        );
    }

    #[test]
    fn azimuth_at_null_out() {
        let h = create_vega();
        assert_eq!(
            siderust_target_azimuth_at(h, paris(), 60000.0, ptr::null_mut()),
            SiderustStatus::NullPointer
        );
        unsafe { siderust_target_free(h) };
    }

    #[test]
    fn azimuth_crossings_ok() {
        let h = create_vega();
        let mut out: *mut SiderustAzimuthCrossingEvent = ptr::null_mut();
        let mut count = 0usize;
        let s = siderust_target_azimuth_crossings(
            h,
            paris(),
            one_day_window(),
            180.0,
            default_opts(),
            &mut out,
            &mut count,
        );
        assert_eq!(s, SiderustStatus::Ok);
        unsafe {
            crate::azimuth::siderust_azimuth_crossings_free(out, count);
            siderust_target_free(h);
        }
    }

    #[test]
    fn azimuth_crossings_null_handle() {
        let mut out: *mut SiderustAzimuthCrossingEvent = ptr::null_mut();
        let mut count = 0usize;
        assert_eq!(
            siderust_target_azimuth_crossings(
                ptr::null(),
                paris(),
                one_day_window(),
                180.0,
                default_opts(),
                &mut out,
                &mut count
            ),
            SiderustStatus::NullPointer
        );
    }

    #[test]
    fn azimuth_crossings_invalid_window() {
        let h = create_vega();
        let bad = TempochPeriodMjd {
            start_mjd: 60001.0,
            end_mjd: 60000.0,
        };
        let mut out: *mut SiderustAzimuthCrossingEvent = ptr::null_mut();
        let mut count = 0usize;
        assert_eq!(
            siderust_target_azimuth_crossings(
                h,
                paris(),
                bad,
                180.0,
                default_opts(),
                &mut out,
                &mut count
            ),
            SiderustStatus::InvalidPeriod
        );
        unsafe { siderust_target_free(h) };
    }

    // ── Below threshold ───────────────────────────────────────────────────

    #[test]
    fn below_threshold_ok() {
        let h = create_vega();
        let mut out: *mut TempochPeriodMjd = ptr::null_mut();
        let mut count = 0usize;
        let s = siderust_target_below_threshold(
            h,
            paris(),
            one_day_window(),
            0.0,
            default_opts(),
            &mut out,
            &mut count,
        );
        assert_eq!(s, SiderustStatus::Ok);
        unsafe {
            crate::altitude::siderust_periods_free(out, count);
            siderust_target_free(h);
        }
    }

    #[test]
    fn below_threshold_null_handle() {
        let mut out: *mut TempochPeriodMjd = ptr::null_mut();
        let mut count = 0usize;
        assert_eq!(
            siderust_target_below_threshold(
                ptr::null(),
                paris(),
                one_day_window(),
                0.0,
                default_opts(),
                &mut out,
                &mut count
            ),
            SiderustStatus::NullPointer
        );
    }

    #[test]
    fn below_threshold_invalid_window() {
        let h = create_vega();
        let bad = TempochPeriodMjd {
            start_mjd: 60001.0,
            end_mjd: 60000.0,
        };
        let mut out: *mut TempochPeriodMjd = ptr::null_mut();
        let mut count = 0usize;
        assert_eq!(
            siderust_target_below_threshold(
                h,
                paris(),
                bad,
                0.0,
                default_opts(),
                &mut out,
                &mut count
            ),
            SiderustStatus::InvalidPeriod
        );
        unsafe { siderust_target_free(h) };
    }
}
