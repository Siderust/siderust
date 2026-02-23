// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! FFI bindings for coordinate types and transformations.
//!
//! Coordinates in siderust are heavily generic (parameterized by frame, center, unit).
//! This module monomorphizes them into concrete C-callable functions using runtime
//! frame/center enum dispatch.

use crate::error::SiderustStatus;
use crate::types::*;
use qtty::*;
use siderust::coordinates::centers::Geodetic;
use siderust::coordinates::frames::{
    EclipticMeanJ2000, EquatorialMeanJ2000, EquatorialMeanOfDate, EquatorialTrueOfDate,
    ReferenceFrame, ECEF, ICRS,
};
use siderust::coordinates::spherical;
use siderust::coordinates::transform::{DirectionAstroExt, SphericalDirectionAstroExt};
use siderust::time::JulianDate;

// ═══════════════════════════════════════════════════════════════════════════
// Spherical Direction — frame transforms
// ═══════════════════════════════════════════════════════════════════════════

/// Trait object proxy for spherical direction transforms.
/// This allows runtime dispatch over the frame type parameter.
///
/// Supported source frames: ICRS, EclipticMeanJ2000, EquatorialMeanJ2000,
/// EquatorialMeanOfDate, EquatorialTrueOfDate.
///
/// Galactic, Horizontal, ECEF and other frames are not supported as source
/// frames because no FrameRotationProvider exists for them.
trait SphericalDirProxy {
    fn to_icrs(&self, jd: &JulianDate) -> (f64, f64);
    fn to_ecliptic_j2000(&self, jd: &JulianDate) -> (f64, f64);
    fn to_equatorial_j2000(&self, jd: &JulianDate) -> (f64, f64);
    fn to_equatorial_mean_of_date(&self, jd: &JulianDate) -> (f64, f64);
    fn to_equatorial_true_of_date(&self, jd: &JulianDate) -> (f64, f64);
    fn to_horizontal(&self, jd: &JulianDate, site: &Geodetic<ECEF>) -> (f64, f64);
}

/// Helper: extract (polar_deg, azimuth_deg) from a spherical Direction.
#[inline]
fn sph_to_pair<F: ReferenceFrame>(d: &spherical::Direction<F>) -> (f64, f64) {
    (d.polar.value(), d.azimuth.value())
}

/// Implement the proxy by always routing through ICRS as a hub.
///
/// For any frame F → target T, we do: F → ICRS → T.
/// When F == ICRS or T == ICRS, only one step is needed (the compiler
/// will optimize out the identity transform for same-frame).
macro_rules! impl_sph_dir_proxy {
    ($frame_ty:ty) => {
        impl SphericalDirProxy for spherical::Direction<$frame_ty> {
            fn to_icrs(&self, jd: &JulianDate) -> (f64, f64) {
                let r = SphericalDirectionAstroExt::to_frame::<ICRS>(self, jd);
                sph_to_pair(&r)
            }
            fn to_ecliptic_j2000(&self, jd: &JulianDate) -> (f64, f64) {
                // Route through ICRS: F→ICRS→EclipticMeanJ2000
                let icrs = SphericalDirectionAstroExt::to_frame::<ICRS>(self, jd);
                let r = SphericalDirectionAstroExt::to_frame::<EclipticMeanJ2000>(&icrs, jd);
                sph_to_pair(&r)
            }
            fn to_equatorial_j2000(&self, jd: &JulianDate) -> (f64, f64) {
                let icrs = SphericalDirectionAstroExt::to_frame::<ICRS>(self, jd);
                let r = SphericalDirectionAstroExt::to_frame::<EquatorialMeanJ2000>(&icrs, jd);
                sph_to_pair(&r)
            }
            fn to_equatorial_mean_of_date(&self, jd: &JulianDate) -> (f64, f64) {
                let icrs = SphericalDirectionAstroExt::to_frame::<ICRS>(self, jd);
                let r = SphericalDirectionAstroExt::to_frame::<EquatorialMeanOfDate>(&icrs, jd);
                sph_to_pair(&r)
            }
            fn to_equatorial_true_of_date(&self, jd: &JulianDate) -> (f64, f64) {
                let icrs = SphericalDirectionAstroExt::to_frame::<ICRS>(self, jd);
                let r = SphericalDirectionAstroExt::to_frame::<EquatorialTrueOfDate>(&icrs, jd);
                sph_to_pair(&r)
            }
            fn to_horizontal(&self, jd: &JulianDate, site: &Geodetic<ECEF>) -> (f64, f64) {
                // Route: F → ICRS → EquatorialTrueOfDate (cartesian) → Horizontal
                let icrs = SphericalDirectionAstroExt::to_frame::<ICRS>(self, jd);
                let eq_tod =
                    SphericalDirectionAstroExt::to_frame::<EquatorialTrueOfDate>(&icrs, jd);
                let cart = eq_tod.to_cartesian();
                let hz = DirectionAstroExt::to_horizontal(&cart, jd, site);
                let hz_sph = spherical::Direction::from_cartesian(&hz);
                sph_to_pair(&hz_sph)
            }
        }
    };
}

impl_sph_dir_proxy!(ICRS);
impl_sph_dir_proxy!(EclipticMeanJ2000);
impl_sph_dir_proxy!(EquatorialMeanJ2000);
impl_sph_dir_proxy!(EquatorialMeanOfDate);
impl_sph_dir_proxy!(EquatorialTrueOfDate);

/// Helper: construct a direction in a given frame as a trait object.
fn make_sph_dir_in_frame(
    frame: SiderustFrame,
    polar_deg: f64,
    azimuth_deg: f64,
) -> Result<Box<dyn SphericalDirProxy>, SiderustStatus> {
    match frame {
        SiderustFrame::ICRS => Ok(Box::new(spherical::Direction::<ICRS>::new(
            Degrees::new(azimuth_deg),
            Degrees::new(polar_deg),
        ))),
        SiderustFrame::EclipticMeanJ2000 => {
            Ok(Box::new(spherical::Direction::<EclipticMeanJ2000>::new(
                Degrees::new(azimuth_deg),
                Degrees::new(polar_deg),
            )))
        }
        SiderustFrame::EquatorialMeanJ2000 => {
            Ok(Box::new(spherical::Direction::<EquatorialMeanJ2000>::new(
                Degrees::new(azimuth_deg),
                Degrees::new(polar_deg),
            )))
        }
        SiderustFrame::EquatorialMeanOfDate => {
            Ok(Box::new(spherical::Direction::<EquatorialMeanOfDate>::new(
                Degrees::new(azimuth_deg),
                Degrees::new(polar_deg),
            )))
        }
        SiderustFrame::EquatorialTrueOfDate => {
            Ok(Box::new(spherical::Direction::<EquatorialTrueOfDate>::new(
                Degrees::new(azimuth_deg),
                Degrees::new(polar_deg),
            )))
        }
        _ => Err(SiderustStatus::InvalidFrame),
    }
}

/// Transform a spherical direction from one frame to another.
///
/// `jd` is required for time-dependent frames (mean/true of date, horizontal).
/// For time-independent transforms between fixed-epoch frames, the JD value
/// is ignored but must still be provided.
///
/// Supported source/target frames: ICRS, EclipticMeanJ2000, EquatorialMeanJ2000,
/// EquatorialMeanOfDate, EquatorialTrueOfDate.
///
/// For conversion to Horizontal, use `siderust_spherical_dir_to_horizontal`.
#[no_mangle]
pub extern "C" fn siderust_spherical_dir_transform_frame(
    polar_deg: f64,
    azimuth_deg: f64,
    src_frame: SiderustFrame,
    dst_frame: SiderustFrame,
    jd: f64,
    out: *mut SiderustSphericalDir,
) -> SiderustStatus {
    if out.is_null() {
        return SiderustStatus::NullPointer;
    }

    let dir = match make_sph_dir_in_frame(src_frame, polar_deg, azimuth_deg) {
        Ok(d) => d,
        Err(e) => return e,
    };

    let jd_val = JulianDate::new(jd);

    let (out_polar, out_azimuth) = match dst_frame {
        SiderustFrame::ICRS => dir.to_icrs(&jd_val),
        SiderustFrame::EclipticMeanJ2000 => dir.to_ecliptic_j2000(&jd_val),
        SiderustFrame::EquatorialMeanJ2000 => dir.to_equatorial_j2000(&jd_val),
        SiderustFrame::EquatorialMeanOfDate => dir.to_equatorial_mean_of_date(&jd_val),
        SiderustFrame::EquatorialTrueOfDate => dir.to_equatorial_true_of_date(&jd_val),
        SiderustFrame::Horizontal => {
            // Horizontal needs an observer — use siderust_spherical_dir_to_horizontal instead
            return SiderustStatus::InvalidFrame;
        }
        _ => return SiderustStatus::InvalidFrame,
    };

    unsafe {
        *out = SiderustSphericalDir {
            polar_deg: out_polar,
            azimuth_deg: out_azimuth,
            frame: dst_frame,
        };
    }
    SiderustStatus::Ok
}

/// Transform a spherical direction to the horizontal (alt-az) frame.
///
/// Requires an observer location (geodetic WGS84) and a Julian Date.
/// Supported source frames: ICRS, EclipticMeanJ2000, EquatorialMeanJ2000,
/// EquatorialMeanOfDate, EquatorialTrueOfDate.
#[no_mangle]
pub extern "C" fn siderust_spherical_dir_to_horizontal(
    polar_deg: f64,
    azimuth_deg: f64,
    src_frame: SiderustFrame,
    jd: f64,
    observer: SiderustGeodetict,
    out: *mut SiderustSphericalDir,
) -> SiderustStatus {
    if out.is_null() {
        return SiderustStatus::NullPointer;
    }

    let dir = match make_sph_dir_in_frame(src_frame, polar_deg, azimuth_deg) {
        Ok(d) => d,
        Err(e) => return e,
    };

    let jd_val = JulianDate::new(jd);
    let site = observer.to_rust();
    let (az, alt) = dir.to_horizontal(&jd_val, &site);

    unsafe {
        *out = SiderustSphericalDir {
            polar_deg: alt,
            azimuth_deg: az,
            frame: SiderustFrame::Horizontal,
        };
    }
    SiderustStatus::Ok
}

/// Create a geodetic position and convert to ECEF Cartesian.
#[no_mangle]
pub extern "C" fn siderust_geodetic_to_cartesian_ecef(
    geodetic: SiderustGeodetict,
    out: *mut SiderustCartesianPos,
) -> SiderustStatus {
    if out.is_null() {
        return SiderustStatus::NullPointer;
    }
    let g = geodetic.to_rust();
    let cart = g.to_cartesian::<Meter>();

    unsafe {
        *out = SiderustCartesianPos {
            x: cart.x().value(),
            y: cart.y().value(),
            z: cart.z().value(),
            frame: SiderustFrame::ECEF,
            center: SiderustCenter::Geocentric,
        };
    }
    SiderustStatus::Ok
}
