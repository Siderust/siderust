// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Target FFI, opaque handles for celestial targets.
//!
//! The C ABI stores the original tagged target payload for lossless roundtrips
//! while deriving one canonical ICRS direction internally so the unified
//! subject-based query APIs can consume every supported target form through a
//! single handle type.

use crate::error::SiderustStatus;
use crate::types::*;
use qtty::*;
use siderust::astro::proper_motion::ProperMotion;
use siderust::coordinates::frames::{
    EclipticMeanJ2000, EquatorialMeanJ2000, EquatorialMeanOfDate, EquatorialTrueOfDate, ICRF, ICRS,
};
use siderust::coordinates::spherical;
use siderust::coordinates::transform::SphericalDirectionAstroExt;
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
/// - direction-only coordinates,
/// - spherical positions with distance metadata,
/// - cartesian positions with frame/center/unit metadata,
/// - and optional proper motion payloads.
pub struct SiderustGenericTarget {
    /// Canonical line-of-sight representation used by the unified subject APIs.
    pub(crate) inner: CoordinateWithPM<spherical::direction::ICRS>,
    /// Original coordinate data for property access and adapter roundtrips.
    pub(crate) data: SiderustGenericTargetData,
}

#[inline]
fn default_proper_motion() -> SiderustProperMotion {
    SiderustProperMotion {
        pm_ra_deg_yr: 0.0,
        pm_dec_deg_yr: 0.0,
        ra_convention: SiderustRaConvention::MuAlpha,
    }
}

fn proper_motion_from_ffi(data: &SiderustGenericTargetData) -> Option<ProperMotion> {
    if !data.has_proper_motion {
        return None;
    }

    Some(match data.proper_motion.ra_convention {
        SiderustRaConvention::MuAlpha => ProperMotion::from_mu_alpha::<DegreePerYear>(
            DegreesPerYear::new(data.proper_motion.pm_ra_deg_yr),
            DegreesPerYear::new(data.proper_motion.pm_dec_deg_yr),
        ),
        SiderustRaConvention::MuAlphaStar => ProperMotion::from_mu_alpha_star::<DegreePerYear>(
            DegreesPerYear::new(data.proper_motion.pm_ra_deg_yr),
            DegreesPerYear::new(data.proper_motion.pm_dec_deg_yr),
        ),
    })
}

fn icrs_direction_from_spherical_dir(
    dir: SiderustSphericalDir,
    jd: JulianDate,
) -> Result<spherical::direction::ICRS, SiderustStatus> {
    match dir.frame {
        SiderustFrame::ICRS => Ok(spherical::Direction::<ICRS>::new(
            Degrees::new(dir.azimuth_deg),
            Degrees::new(dir.polar_deg),
        )),
        SiderustFrame::ICRF => Ok(spherical::Direction::<ICRF>::new_raw(
            Degrees::new(dir.azimuth_deg),
            Degrees::new(dir.polar_deg),
        )
        .to_frame::<ICRS>(&jd)),
        SiderustFrame::EclipticMeanJ2000 => Ok(spherical::Direction::<EclipticMeanJ2000>::new(
            Degrees::new(dir.azimuth_deg),
            Degrees::new(dir.polar_deg),
        )
        .to_frame::<ICRS>(&jd)),
        SiderustFrame::EquatorialMeanJ2000 => Ok(spherical::Direction::<EquatorialMeanJ2000>::new(
            Degrees::new(dir.azimuth_deg),
            Degrees::new(dir.polar_deg),
        )
        .to_frame::<ICRS>(&jd)),
        SiderustFrame::EquatorialMeanOfDate => {
            Ok(spherical::Direction::<EquatorialMeanOfDate>::new(
                Degrees::new(dir.azimuth_deg),
                Degrees::new(dir.polar_deg),
            )
            .to_frame::<ICRS>(&jd))
        }
        SiderustFrame::EquatorialTrueOfDate => {
            Ok(spherical::Direction::<EquatorialTrueOfDate>::new(
                Degrees::new(dir.azimuth_deg),
                Degrees::new(dir.polar_deg),
            )
            .to_frame::<ICRS>(&jd))
        }
        _ => Err(SiderustStatus::InvalidFrame),
    }
}

fn direction_payload_from_spherical_pos(
    pos: SiderustSphericalPos,
) -> Result<SiderustSphericalDir, SiderustStatus> {
    if !pos.lon_deg.is_finite() || !pos.lat_deg.is_finite() || !pos.distance.is_finite() {
        return Err(SiderustStatus::InvalidArgument);
    }
    if pos.distance <= 0.0 {
        return Err(SiderustStatus::InvalidArgument);
    }

    Ok(SiderustSphericalDir {
        polar_deg: pos.lat_deg,
        azimuth_deg: pos.lon_deg,
        frame: pos.frame,
    })
}

fn direction_payload_from_cartesian_pos(
    pos: SiderustCartesianPos,
) -> Result<SiderustSphericalDir, SiderustStatus> {
    if !pos.x.is_finite() || !pos.y.is_finite() || !pos.z.is_finite() {
        return Err(SiderustStatus::InvalidArgument);
    }

    let norm = (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z).sqrt();
    if !norm.is_finite() || norm <= f64::EPSILON {
        return Err(SiderustStatus::InvalidArgument);
    }

    let x = pos.x / norm;
    let y = pos.y / norm;
    let z = pos.z / norm;
    let azimuth_deg = y.atan2(x).to_degrees().rem_euclid(360.0);
    let polar_deg = z.atan2((x * x + y * y).sqrt()).to_degrees();

    Ok(SiderustSphericalDir {
        polar_deg,
        azimuth_deg,
        frame: pos.frame,
    })
}

fn build_inner_target(
    data: &SiderustGenericTargetData,
) -> Result<CoordinateWithPM<spherical::direction::ICRS>, SiderustStatus> {
    let epoch = JulianDate::new(data.epoch_jd);
    let source_dir = match data.kind {
        SiderustTargetCoordKind::SphericalDir => unsafe { data.coord.spherical_dir },
        SiderustTargetCoordKind::SphericalPos => {
            direction_payload_from_spherical_pos(unsafe { data.coord.spherical_pos })?
        }
        SiderustTargetCoordKind::CartesianPos => {
            direction_payload_from_cartesian_pos(unsafe { data.coord.cartesian_pos })?
        }
    };

    let icrs_dir = icrs_direction_from_spherical_dir(source_dir, epoch)?;
    Ok(CoordinateWithPM::new_raw(
        icrs_dir,
        epoch,
        proper_motion_from_ffi(data),
    ))
}

/// Create a generic target from the tagged target payload.
///
/// This is the canonical constructor for the compact FFI target ABI. Existing
/// convenience helpers remain available as compatibility shims and forward here.
#[no_mangle]
pub extern "C" fn siderust_generic_target_create(
    data: SiderustGenericTargetData,
    out: *mut *mut SiderustGenericTarget,
) -> SiderustStatus {
    ffi_guard! {{
        if out.is_null() {
            return SiderustStatus::NullPointer;
        }

        let inner = match build_inner_target(&data) {
            Ok(inner) => inner,
            Err(status) => return status,
        };

        let handle = Box::new(SiderustGenericTarget { inner, data });
        unsafe { *out = Box::into_raw(handle) };
        SiderustStatus::Ok
    }}
}

/// Create a generic target from an ICRS spherical direction.
#[no_mangle]
pub extern "C" fn siderust_generic_target_create_icrs(
    ra_deg: f64,
    dec_deg: f64,
    epoch_jd: f64,
    out: *mut *mut SiderustGenericTarget,
) -> SiderustStatus {
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
        proper_motion: default_proper_motion(),
    };
    siderust_generic_target_create(data, out)
}

/// Create a generic target from an ICRS direction with proper motion.
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
    siderust_generic_target_create(data, out)
}

/// Free a generic target handle.
///
/// # Safety
/// `handle` must have been allocated by a `siderust_generic_target_create*`
/// function in this crate and must not have been freed before.
#[no_mangle]
pub unsafe extern "C" fn siderust_generic_target_free(handle: *mut SiderustGenericTarget) {
    if !handle.is_null() {
        drop(unsafe { Box::from_raw(handle) });
    }
}

/// Get the coordinate payload from a generic target.
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::ptr;

    fn dummy_target_data() -> SiderustGenericTargetData {
        SiderustGenericTargetData {
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
            epoch_jd: 0.0,
            has_proper_motion: false,
            _pad2: [0; 7],
            proper_motion: default_proper_motion(),
        }
    }

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
        let mut out = dummy_target_data();
        assert_eq!(
            siderust_generic_target_get_data(h, &mut out),
            SiderustStatus::Ok
        );
        assert_eq!(out.kind, SiderustTargetCoordKind::SphericalDir);
        let coord = unsafe { out.coord.spherical_dir };
        assert_eq!(coord.frame, SiderustFrame::ICRS);
        assert!((coord.azimuth_deg - 279.23).abs() < 1e-10);
        assert!((coord.polar_deg - 38.78).abs() < 1e-10);
        assert_eq!(out.epoch_jd, 2451545.0);
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
        let mut out = dummy_target_data();
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
    fn create_from_spherical_position_preserves_payload() {
        let data = SiderustGenericTargetData {
            kind: SiderustTargetCoordKind::SphericalPos,
            _pad1: [0; 4],
            coord: SiderustTargetCoordUnion {
                spherical_pos: SiderustSphericalPos {
                    lon_deg: 88.792939,
                    lat_deg: 7.407064,
                    distance: 548.0,
                    frame: SiderustFrame::EquatorialMeanJ2000,
                    center: SiderustCenter::Geocentric,
                    length_unit: SiderustLengthUnit::LightYear,
                },
            },
            epoch_jd: 2451545.0,
            has_proper_motion: false,
            _pad2: [0; 7],
            proper_motion: default_proper_motion(),
        };

        let mut handle: *mut SiderustGenericTarget = ptr::null_mut();
        assert_eq!(
            siderust_generic_target_create(data, &mut handle),
            SiderustStatus::Ok
        );

        let mut out = dummy_target_data();
        assert_eq!(
            siderust_generic_target_get_data(handle, &mut out),
            SiderustStatus::Ok
        );
        assert_eq!(out.kind, SiderustTargetCoordKind::SphericalPos);
        let pos = unsafe { out.coord.spherical_pos };
        assert_eq!(pos.frame, SiderustFrame::EquatorialMeanJ2000);
        assert_eq!(pos.center, SiderustCenter::Geocentric);
        assert_eq!(pos.length_unit, SiderustLengthUnit::LightYear);
        assert!((pos.distance - 548.0).abs() < 1e-12);

        unsafe { siderust_generic_target_free(handle) };
    }

    #[test]
    fn create_from_cartesian_position_preserves_payload() {
        let data = SiderustGenericTargetData {
            kind: SiderustTargetCoordKind::CartesianPos,
            _pad1: [0; 4],
            coord: SiderustTargetCoordUnion {
                cartesian_pos: SiderustCartesianPos {
                    x: 0.1,
                    y: 0.2,
                    z: 0.3,
                    frame: SiderustFrame::EquatorialMeanJ2000,
                    center: SiderustCenter::Barycentric,
                    length_unit: SiderustLengthUnit::AU,
                },
            },
            epoch_jd: 2451545.0,
            has_proper_motion: false,
            _pad2: [0; 7],
            proper_motion: default_proper_motion(),
        };

        let mut handle: *mut SiderustGenericTarget = ptr::null_mut();
        assert_eq!(
            siderust_generic_target_create(data, &mut handle),
            SiderustStatus::Ok
        );

        let mut out = dummy_target_data();
        assert_eq!(
            siderust_generic_target_get_data(handle, &mut out),
            SiderustStatus::Ok
        );
        assert_eq!(out.kind, SiderustTargetCoordKind::CartesianPos);
        let pos = unsafe { out.coord.cartesian_pos };
        assert_eq!(pos.frame, SiderustFrame::EquatorialMeanJ2000);
        assert_eq!(pos.center, SiderustCenter::Barycentric);
        assert_eq!(pos.length_unit, SiderustLengthUnit::AU);

        unsafe { siderust_generic_target_free(handle) };
    }

    #[test]
    fn create_rejects_invalid_position_payloads() {
        let mut handle: *mut SiderustGenericTarget = ptr::null_mut();
        let bad_cart = SiderustGenericTargetData {
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
            epoch_jd: 2451545.0,
            has_proper_motion: false,
            _pad2: [0; 7],
            proper_motion: default_proper_motion(),
        };
        assert_eq!(
            siderust_generic_target_create(bad_cart, &mut handle),
            SiderustStatus::InvalidArgument
        );
        assert!(handle.is_null());

        let bad_frame = SiderustGenericTargetData {
            kind: SiderustTargetCoordKind::SphericalDir,
            _pad1: [0; 4],
            coord: SiderustTargetCoordUnion {
                spherical_dir: SiderustSphericalDir {
                    polar_deg: 10.0,
                    azimuth_deg: 20.0,
                    frame: SiderustFrame::Horizontal,
                },
            },
            epoch_jd: 2451545.0,
            has_proper_motion: false,
            _pad2: [0; 7],
            proper_motion: default_proper_motion(),
        };
        assert_eq!(
            siderust_generic_target_create(bad_frame, &mut handle),
            SiderustStatus::InvalidFrame
        );
        assert!(handle.is_null());
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
        let mut out = dummy_target_data();
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
