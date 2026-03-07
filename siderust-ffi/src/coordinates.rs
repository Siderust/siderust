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
use siderust::calculus::ephemeris::{Ephemeris, Vsop87Ephemeris};
use siderust::coordinates::centers::{Barycentric, Geodetic};
use siderust::coordinates::frames::{
    EclipticMeanJ2000, EquatorialMeanJ2000, EquatorialMeanOfDate, EquatorialTrueOfDate,
    ReferenceFrame, ECEF, ICRF, ICRS,
};
use siderust::coordinates::transform::{
    DirectionAstroExt, PositionAstroExt, SphericalDirectionAstroExt,
};
use siderust::coordinates::{cartesian, spherical};
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
    ffi_guard! {{
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

    }}
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
    ffi_guard! {{
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

    }}
}

// ═══════════════════════════════════════════════════════════════════════════
// Cartesian Direction — frame transforms
// ═══════════════════════════════════════════════════════════════════════════

/// Trait object proxy for cartesian direction frame transforms.
///
/// Mirrors `SphericalDirProxy` but for cartesian unit-vector directions.
trait CartesianDirProxy {
    fn to_icrs(&self, jd: &JulianDate) -> (f64, f64, f64);
    fn to_ecliptic_j2000(&self, jd: &JulianDate) -> (f64, f64, f64);
    fn to_equatorial_j2000(&self, jd: &JulianDate) -> (f64, f64, f64);
    fn to_equatorial_mean_of_date(&self, jd: &JulianDate) -> (f64, f64, f64);
    fn to_equatorial_true_of_date(&self, jd: &JulianDate) -> (f64, f64, f64);
}

#[inline]
fn cart_dir_to_triple<F: ReferenceFrame>(d: &cartesian::Direction<F>) -> (f64, f64, f64) {
    (d.x(), d.y(), d.z())
}

macro_rules! impl_cart_dir_proxy {
    ($frame_ty:ty) => {
        impl CartesianDirProxy for cartesian::Direction<$frame_ty> {
            fn to_icrs(&self, jd: &JulianDate) -> (f64, f64, f64) {
                cart_dir_to_triple(&DirectionAstroExt::to_frame::<ICRS>(self, jd))
            }
            fn to_ecliptic_j2000(&self, jd: &JulianDate) -> (f64, f64, f64) {
                let icrs = DirectionAstroExt::to_frame::<ICRS>(self, jd);
                cart_dir_to_triple(&DirectionAstroExt::to_frame::<EclipticMeanJ2000>(&icrs, jd))
            }
            fn to_equatorial_j2000(&self, jd: &JulianDate) -> (f64, f64, f64) {
                let icrs = DirectionAstroExt::to_frame::<ICRS>(self, jd);
                cart_dir_to_triple(&DirectionAstroExt::to_frame::<EquatorialMeanJ2000>(
                    &icrs, jd,
                ))
            }
            fn to_equatorial_mean_of_date(&self, jd: &JulianDate) -> (f64, f64, f64) {
                let icrs = DirectionAstroExt::to_frame::<ICRS>(self, jd);
                cart_dir_to_triple(&DirectionAstroExt::to_frame::<EquatorialMeanOfDate>(
                    &icrs, jd,
                ))
            }
            fn to_equatorial_true_of_date(&self, jd: &JulianDate) -> (f64, f64, f64) {
                let icrs = DirectionAstroExt::to_frame::<ICRS>(self, jd);
                cart_dir_to_triple(&DirectionAstroExt::to_frame::<EquatorialTrueOfDate>(
                    &icrs, jd,
                ))
            }
        }
    };
}

impl_cart_dir_proxy!(ICRS);
impl_cart_dir_proxy!(EclipticMeanJ2000);
impl_cart_dir_proxy!(EquatorialMeanJ2000);
impl_cart_dir_proxy!(EquatorialMeanOfDate);
impl_cart_dir_proxy!(EquatorialTrueOfDate);

fn make_cart_dir_in_frame(
    frame: SiderustFrame,
    x: f64,
    y: f64,
    z: f64,
) -> Result<Box<dyn CartesianDirProxy>, SiderustStatus> {
    match frame {
        SiderustFrame::ICRS => Ok(Box::new(cartesian::Direction::<ICRS>::new_unchecked(
            x, y, z,
        ))),
        SiderustFrame::EclipticMeanJ2000 => Ok(Box::new(
            cartesian::Direction::<EclipticMeanJ2000>::new_unchecked(x, y, z),
        )),
        SiderustFrame::EquatorialMeanJ2000 => Ok(Box::new(cartesian::Direction::<
            EquatorialMeanJ2000,
        >::new_unchecked(x, y, z))),
        SiderustFrame::EquatorialMeanOfDate => Ok(Box::new(cartesian::Direction::<
            EquatorialMeanOfDate,
        >::new_unchecked(x, y, z))),
        SiderustFrame::EquatorialTrueOfDate => Ok(Box::new(cartesian::Direction::<
            EquatorialTrueOfDate,
        >::new_unchecked(x, y, z))),
        _ => Err(SiderustStatus::InvalidFrame),
    }
}

/// Transform a Cartesian unit-vector direction from one frame to another.
///
/// The input (x, y, z) need not be normalised; the rotation preserves the
/// vector length, so the output has the same magnitude as the input.
///
/// Supported source/target frames: ICRS, EclipticMeanJ2000, EquatorialMeanJ2000,
/// EquatorialMeanOfDate, EquatorialTrueOfDate.
#[no_mangle]
pub extern "C" fn siderust_cartesian_dir_transform_frame(
    x: f64,
    y: f64,
    z: f64,
    src_frame: SiderustFrame,
    dst_frame: SiderustFrame,
    jd: f64,
    out: *mut SiderustCartesianPos,
) -> SiderustStatus {
    ffi_guard! {{
        if out.is_null() {
            return SiderustStatus::NullPointer;
        }

        let dir = match make_cart_dir_in_frame(src_frame, x, y, z) {
            Ok(d) => d,
            Err(e) => return e,
        };

        let jd_val = JulianDate::new(jd);

        let (ox, oy, oz) = match dst_frame {
            SiderustFrame::ICRS => dir.to_icrs(&jd_val),
            SiderustFrame::EclipticMeanJ2000 => dir.to_ecliptic_j2000(&jd_val),
            SiderustFrame::EquatorialMeanJ2000 => dir.to_equatorial_j2000(&jd_val),
            SiderustFrame::EquatorialMeanOfDate => dir.to_equatorial_mean_of_date(&jd_val),
            SiderustFrame::EquatorialTrueOfDate => dir.to_equatorial_true_of_date(&jd_val),
            _ => return SiderustStatus::InvalidFrame,
        };

        unsafe {
            *out = SiderustCartesianPos {
                x: ox,
                y: oy,
                z: oz,
                frame: dst_frame,
                center: SiderustCenter::Barycentric, // directions have no center
            };
        }
        SiderustStatus::Ok

    }}
}

// ═══════════════════════════════════════════════════════════════════════════
// Cartesian Position — frame transforms
// ═══════════════════════════════════════════════════════════════════════════

/// Trait object proxy for Cartesian position frame transforms.
///
/// Works with a fixed internal unit (AU) for the rotation arithmetic;
/// the result x/y/z are in the same implicit unit as the input.
trait CartesianPosProxy {
    fn to_icrs(&self, jd: &JulianDate) -> (f64, f64, f64);
    fn to_icrf(&self, jd: &JulianDate) -> (f64, f64, f64);
    fn to_ecliptic_j2000(&self, jd: &JulianDate) -> (f64, f64, f64);
    fn to_equatorial_j2000(&self, jd: &JulianDate) -> (f64, f64, f64);
    fn to_equatorial_mean_of_date(&self, jd: &JulianDate) -> (f64, f64, f64);
    fn to_equatorial_true_of_date(&self, jd: &JulianDate) -> (f64, f64, f64);
}

/// Helper: extract raw (x, y, z) from a Cartesian position.
#[inline]
fn cart_pos_to_triple<F: siderust::coordinates::frames::ReferenceFrame>(
    p: &cartesian::Position<Barycentric, F, AstronomicalUnit>,
) -> (f64, f64, f64) {
    (p.x().value(), p.y().value(), p.z().value())
}

macro_rules! impl_cart_pos_proxy {
    ($frame_ty:ty) => {
        impl CartesianPosProxy for cartesian::Position<Barycentric, $frame_ty, AstronomicalUnit> {
            fn to_icrs(&self, jd: &JulianDate) -> (f64, f64, f64) {
                cart_pos_to_triple(&PositionAstroExt::to_frame::<ICRS>(self, jd))
            }
            fn to_icrf(&self, jd: &JulianDate) -> (f64, f64, f64) {
                let icrs: cartesian::Position<Barycentric, ICRS, AstronomicalUnit> =
                    PositionAstroExt::to_frame(self, jd);
                cart_pos_to_triple(&PositionAstroExt::to_frame::<ICRF>(&icrs, jd))
            }
            fn to_ecliptic_j2000(&self, jd: &JulianDate) -> (f64, f64, f64) {
                let icrs: cartesian::Position<Barycentric, ICRS, AstronomicalUnit> =
                    PositionAstroExt::to_frame(self, jd);
                cart_pos_to_triple(&PositionAstroExt::to_frame::<EclipticMeanJ2000>(&icrs, jd))
            }
            fn to_equatorial_j2000(&self, jd: &JulianDate) -> (f64, f64, f64) {
                let icrs: cartesian::Position<Barycentric, ICRS, AstronomicalUnit> =
                    PositionAstroExt::to_frame(self, jd);
                cart_pos_to_triple(&PositionAstroExt::to_frame::<EquatorialMeanJ2000>(
                    &icrs, jd,
                ))
            }
            fn to_equatorial_mean_of_date(&self, jd: &JulianDate) -> (f64, f64, f64) {
                let icrs: cartesian::Position<Barycentric, ICRS, AstronomicalUnit> =
                    PositionAstroExt::to_frame(self, jd);
                cart_pos_to_triple(&PositionAstroExt::to_frame::<EquatorialMeanOfDate>(
                    &icrs, jd,
                ))
            }
            fn to_equatorial_true_of_date(&self, jd: &JulianDate) -> (f64, f64, f64) {
                let icrs: cartesian::Position<Barycentric, ICRS, AstronomicalUnit> =
                    PositionAstroExt::to_frame(self, jd);
                cart_pos_to_triple(&PositionAstroExt::to_frame::<EquatorialTrueOfDate>(
                    &icrs, jd,
                ))
            }
        }
    };
}

impl_cart_pos_proxy!(ICRS);
impl_cart_pos_proxy!(ICRF);
impl_cart_pos_proxy!(EclipticMeanJ2000);
impl_cart_pos_proxy!(EquatorialMeanJ2000);
impl_cart_pos_proxy!(EquatorialMeanOfDate);
impl_cart_pos_proxy!(EquatorialTrueOfDate);

fn make_cart_pos_in_frame(
    frame: SiderustFrame,
    x: f64,
    y: f64,
    z: f64,
) -> Result<Box<dyn CartesianPosProxy>, SiderustStatus> {
    match frame {
        SiderustFrame::ICRS => Ok(Box::new(cartesian::Position::<
            Barycentric,
            ICRS,
            AstronomicalUnit,
        >::new(x, y, z))),
        SiderustFrame::ICRF => Ok(Box::new(cartesian::Position::<
            Barycentric,
            ICRF,
            AstronomicalUnit,
        >::new(x, y, z))),
        SiderustFrame::EclipticMeanJ2000 => Ok(Box::new(cartesian::Position::<
            Barycentric,
            EclipticMeanJ2000,
            AstronomicalUnit,
        >::new(x, y, z))),
        SiderustFrame::EquatorialMeanJ2000 => Ok(Box::new(cartesian::Position::<
            Barycentric,
            EquatorialMeanJ2000,
            AstronomicalUnit,
        >::new(x, y, z))),
        SiderustFrame::EquatorialMeanOfDate => Ok(Box::new(cartesian::Position::<
            Barycentric,
            EquatorialMeanOfDate,
            AstronomicalUnit,
        >::new(x, y, z))),
        SiderustFrame::EquatorialTrueOfDate => Ok(Box::new(cartesian::Position::<
            Barycentric,
            EquatorialTrueOfDate,
            AstronomicalUnit,
        >::new(x, y, z))),
        _ => Err(SiderustStatus::InvalidFrame),
    }
}

/// Transform a Cartesian position from one frame to another (frame-only, same center).
///
/// The rotation preserves the vector magnitude.  The `center` field of `pos`
/// is copied unchanged to `out`.
///
/// Supported source/target frames: ICRS, EclipticMeanJ2000, EquatorialMeanJ2000,
/// EquatorialMeanOfDate, EquatorialTrueOfDate.
#[no_mangle]
pub extern "C" fn siderust_cartesian_pos_transform_frame(
    pos: SiderustCartesianPos,
    dst_frame: SiderustFrame,
    jd: f64,
    out: *mut SiderustCartesianPos,
) -> SiderustStatus {
    ffi_guard! {{
        if out.is_null() {
            return SiderustStatus::NullPointer;
        }

        let proxy = match make_cart_pos_in_frame(pos.frame, pos.x, pos.y, pos.z) {
            Ok(p) => p,
            Err(e) => return e,
        };

        let jd_val = JulianDate::new(jd);

        let (ox, oy, oz) = match dst_frame {
            SiderustFrame::ICRS => proxy.to_icrs(&jd_val),
            SiderustFrame::ICRF => proxy.to_icrf(&jd_val),
            SiderustFrame::EclipticMeanJ2000 => proxy.to_ecliptic_j2000(&jd_val),
            SiderustFrame::EquatorialMeanJ2000 => proxy.to_equatorial_j2000(&jd_val),
            SiderustFrame::EquatorialMeanOfDate => proxy.to_equatorial_mean_of_date(&jd_val),
            SiderustFrame::EquatorialTrueOfDate => proxy.to_equatorial_true_of_date(&jd_val),
            _ => return SiderustStatus::InvalidFrame,
        };

        unsafe {
            *out = SiderustCartesianPos {
                x: ox,
                y: oy,
                z: oz,
                frame: dst_frame,
                center: pos.center, // center is unchanged for a frame-only transform
            };
        }
        SiderustStatus::Ok

    }}
}

/// Create a geodetic position and convert to ECEF Cartesian.
#[no_mangle]
pub extern "C" fn siderust_geodetic_to_cartesian_ecef(
    geodetic: SiderustGeodetict,
    out: *mut SiderustCartesianPos,
) -> SiderustStatus {
    ffi_guard! {{
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

    }}
}

// ═══════════════════════════════════════════════════════════════════════════
// Keplerian propagation & bodycentric transforms
// ═══════════════════════════════════════════════════════════════════════════

/// Shift a position (x, y, z) in AU/EclipticMeanJ2000 from one reference center
/// to another, using VSOP87 for Earth-helio and Sun-bary offsets.
///
/// Center codes: 0 = Barycentric, 1 = Heliocentric, 2 = Geocentric.
/// These match `OrbitReferenceCenter` (Barycentric=0, Heliocentric=1, Geocentric=2).
fn shift_center_xyz(x: f64, y: f64, z: f64, from: u8, to: u8, jd: f64) -> (f64, f64, f64) {
    if from == to {
        return (x, y, z);
    }
    let t = JulianDate::new(jd);

    // Sun's barycentric position — used for helio ↔ bary shifts
    let sun_b = Vsop87Ephemeris::sun_barycentric(t);
    let (sb_x, sb_y, sb_z) = (sun_b.x().value(), sun_b.y().value(), sun_b.z().value());

    // Earth's heliocentric position — used for helio ↔ geo shifts
    let earth_h = Vsop87Ephemeris::earth_heliocentric(t);
    let (eh_x, eh_y, eh_z) = (
        earth_h.x().value(),
        earth_h.y().value(),
        earth_h.z().value(),
    );

    // Convert to heliocentric first (hub = heliocentric)
    let (hx, hy, hz) = match from {
        0 => (x - sb_x, y - sb_y, z - sb_z), // bary → helio  (subtract Sun-bary pos)
        1 => (x, y, z),                      // already heliocentric
        2 => (x + eh_x, y + eh_y, z + eh_z), // geo  → helio  (add Earth-helio pos)
        _ => (x, y, z),
    };

    // Convert heliocentric to target center
    match to {
        0 => (hx + sb_x, hy + sb_y, hz + sb_z), // helio → bary
        1 => (hx, hy, hz),                      // helio → helio (no-op)
        2 => (hx - eh_x, hy - eh_y, hz - eh_z), // helio → geo
        _ => (hx, hy, hz),
    }
}

/// Map SiderustCenter enum values to the internal shift_center_xyz codes.
fn center_to_shift_code(c: SiderustCenter) -> Option<u8> {
    match c {
        SiderustCenter::Barycentric => Some(0),
        SiderustCenter::Heliocentric => Some(1),
        SiderustCenter::Geocentric => Some(2),
        _ => None,
    }
}

/// Transform a Cartesian position from one reference center to another
/// (frame-only, same frame).
///
/// Only works for EclipticMeanJ2000 frame positions.  The input frame must
/// already be EclipticMeanJ2000; this function performs only the center shift.
///
/// Supported centers: Barycentric, Heliocentric, Geocentric.
#[no_mangle]
pub extern "C" fn siderust_cartesian_pos_transform_center(
    pos: SiderustCartesianPos,
    dst_center: SiderustCenter,
    jd: f64,
    out: *mut SiderustCartesianPos,
) -> SiderustStatus {
    ffi_guard! {{
        if out.is_null() {
            return SiderustStatus::NullPointer;
        }

        let from_code = match center_to_shift_code(pos.center) {
            Some(c) => c,
            None => return SiderustStatus::InvalidCenter,
        };
        let to_code = match center_to_shift_code(dst_center) {
            Some(c) => c,
            None => return SiderustStatus::InvalidCenter,
        };

        let (ox, oy, oz) = shift_center_xyz(pos.x, pos.y, pos.z, from_code, to_code, jd);

        unsafe {
            *out = SiderustCartesianPos {
                x: ox,
                y: oy,
                z: oz,
                frame: pos.frame,
                center: dst_center,
            };
        }
        SiderustStatus::Ok

    }}
}

/// Compute the Keplerian orbital position at a given Julian date.
///
/// Returns position in EclipticMeanJ2000 frame (AU), where the reference
/// center equals the orbit's own reference center (e.g. heliocentric for a
/// planet's orbit).  The `center` field of `out` is set to `Heliocentric` as
/// a placeholder; callers should interpret it according to `orbit_center` from
/// the associated `siderust_bodycentric_params_t`.
#[no_mangle]
pub extern "C" fn siderust_kepler_position(
    orbit: SiderustOrbit,
    jd: f64,
    out: *mut SiderustCartesianPos,
) -> SiderustStatus {
    ffi_guard! {{
        if out.is_null() {
            return SiderustStatus::NullPointer;
        }
        let pos = orbit.to_rust().kepler_position(JulianDate::new(jd));
        unsafe {
            *out = SiderustCartesianPos {
                x: pos.x().value(),
                y: pos.y().value(),
                z: pos.z().value(),
                frame: SiderustFrame::EclipticMeanJ2000,
                center: SiderustCenter::Heliocentric, // placeholder; real center = orbit_center
            };
        }
        SiderustStatus::Ok

    }}
}

/// Transform a Cartesian position to body-centric coordinates.
///
/// Mirrors Rust's `ToBodycentricExt::to_bodycentric()`.
///
/// `pos`    – source position in EclipticMeanJ2000 / AU.  Center must be
///            Geocentric, Heliocentric, or Barycentric.
/// `params` – Keplerian orbit + reference center of the body.
/// `jd`     – Julian Date for Keplerian propagation and center-shift.
/// `out`    – relative position in EclipticMeanJ2000 / AU with center Bodycentric.
#[no_mangle]
pub extern "C" fn siderust_to_bodycentric(
    pos: SiderustCartesianPos,
    params: SiderustBodycentricParams,
    jd: f64,
    out: *mut SiderustCartesianPos,
) -> SiderustStatus {
    ffi_guard! {{
        if out.is_null() {
            return SiderustStatus::NullPointer;
        }

        // Map SiderustCenter → OrbitReferenceCenter coding (0=Bary,1=Helio,2=Geo)
        let input_center: u8 = match pos.center {
            SiderustCenter::Barycentric => 0,
            SiderustCenter::Heliocentric => 1,
            SiderustCenter::Geocentric => 2,
            _ => return SiderustStatus::InvalidCenter,
        };

        // Keplerian position of the body in its own orbit's reference center
        let body_kep = params.orbit.to_rust().kepler_position(JulianDate::new(jd));
        let (bkx, bky, bkz) = (
            body_kep.x().value(),
            body_kep.y().value(),
            body_kep.z().value(),
        );

        // Shift body position to match the input center
        let (body_x, body_y, body_z) =
            shift_center_xyz(bkx, bky, bkz, params.orbit_center, input_center, jd);

        // Relative position: input – body (vector from body to target)
        unsafe {
            *out = SiderustCartesianPos {
                x: pos.x - body_x,
                y: pos.y - body_y,
                z: pos.z - body_z,
                frame: pos.frame,
                center: SiderustCenter::Bodycentric,
            };
        }
        SiderustStatus::Ok

    }}
}

/// Transform a body-centric position back to geocentric coordinates.
///
/// Mirrors Rust's `FromBodycentricExt::to_geocentric()`.
///
/// `pos`    – body-centric position in EclipticMeanJ2000 / AU.
/// `params` – same orbital parameters used during `siderust_to_bodycentric`.
/// `jd`     – same Julian Date used during `siderust_to_bodycentric`.
/// `out`    – recovered geocentric position in EclipticMeanJ2000 / AU.
#[no_mangle]
pub extern "C" fn siderust_from_bodycentric(
    pos: SiderustCartesianPos,
    params: SiderustBodycentricParams,
    jd: f64,
    out: *mut SiderustCartesianPos,
) -> SiderustStatus {
    ffi_guard! {{
        if out.is_null() {
            return SiderustStatus::NullPointer;
        }

        // Keplerian position of the body in its own orbit's reference center
        let body_kep = params.orbit.to_rust().kepler_position(JulianDate::new(jd));
        let (bkx, bky, bkz) = (
            body_kep.x().value(),
            body_kep.y().value(),
            body_kep.z().value(),
        );

        // Convert body position to geocentric (target center code = 2)
        let (body_geo_x, body_geo_y, body_geo_z) =
            shift_center_xyz(bkx, bky, bkz, params.orbit_center, 2, jd);

        // Recover geocentric: bodycentric + body_geocentric
        unsafe {
            *out = SiderustCartesianPos {
                x: pos.x + body_geo_x,
                y: pos.y + body_geo_y,
                z: pos.z + body_geo_z,
                frame: pos.frame,
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

    fn empty_dir() -> SiderustSphericalDir {
        SiderustSphericalDir {
            polar_deg: 0.0,
            azimuth_deg: 0.0,
            frame: SiderustFrame::ICRS,
        }
    }

    fn empty_cart() -> SiderustCartesianPos {
        SiderustCartesianPos {
            x: 0.0,
            y: 0.0,
            z: 0.0,
            frame: SiderustFrame::ICRS,
            center: SiderustCenter::Geocentric,
        }
    }

    fn paris_observer() -> SiderustGeodetict {
        SiderustGeodetict {
            lon_deg: 2.35,
            lat_deg: 48.85,
            height_m: 35.0,
        }
    }

    // ── Frame transforms ─────────────────────────────────────────────────

    #[test]
    fn icrs_to_icrs_is_identity() {
        let mut out = empty_dir();
        let s = siderust_spherical_dir_transform_frame(
            10.0,
            45.0,
            SiderustFrame::ICRS,
            SiderustFrame::ICRS,
            J2000,
            &mut out,
        );
        assert_eq!(s, SiderustStatus::Ok);
        assert!((out.polar_deg - 10.0).abs() < 1e-8);
        assert!((out.azimuth_deg - 45.0).abs() < 1e-8);
    }

    #[test]
    fn icrs_to_ecliptic_j2000() {
        let mut out = empty_dir();
        let s = siderust_spherical_dir_transform_frame(
            30.0,
            90.0,
            SiderustFrame::ICRS,
            SiderustFrame::EclipticMeanJ2000,
            J2000,
            &mut out,
        );
        assert_eq!(s, SiderustStatus::Ok);
        assert_eq!(out.frame, SiderustFrame::EclipticMeanJ2000);
    }

    #[test]
    fn icrs_to_equatorial_j2000() {
        let mut out = empty_dir();
        let s = siderust_spherical_dir_transform_frame(
            0.0,
            180.0,
            SiderustFrame::ICRS,
            SiderustFrame::EquatorialMeanJ2000,
            J2000,
            &mut out,
        );
        assert_eq!(s, SiderustStatus::Ok);
        assert_eq!(out.frame, SiderustFrame::EquatorialMeanJ2000);
    }

    #[test]
    fn icrs_to_equatorial_mean_of_date() {
        let mut out = empty_dir();
        let s = siderust_spherical_dir_transform_frame(
            0.0,
            0.0,
            SiderustFrame::ICRS,
            SiderustFrame::EquatorialMeanOfDate,
            J2000,
            &mut out,
        );
        assert_eq!(s, SiderustStatus::Ok);
        assert_eq!(out.frame, SiderustFrame::EquatorialMeanOfDate);
    }

    #[test]
    fn icrs_to_equatorial_true_of_date() {
        let mut out = empty_dir();
        let s = siderust_spherical_dir_transform_frame(
            0.0,
            0.0,
            SiderustFrame::ICRS,
            SiderustFrame::EquatorialTrueOfDate,
            J2000,
            &mut out,
        );
        assert_eq!(s, SiderustStatus::Ok);
        assert_eq!(out.frame, SiderustFrame::EquatorialTrueOfDate);
    }

    #[test]
    fn ecliptic_j2000_to_icrs() {
        let mut out = empty_dir();
        let s = siderust_spherical_dir_transform_frame(
            10.0,
            90.0,
            SiderustFrame::EclipticMeanJ2000,
            SiderustFrame::ICRS,
            J2000,
            &mut out,
        );
        assert_eq!(s, SiderustStatus::Ok);
    }

    #[test]
    fn equatorial_j2000_to_icrs() {
        let mut out = empty_dir();
        let s = siderust_spherical_dir_transform_frame(
            10.0,
            90.0,
            SiderustFrame::EquatorialMeanJ2000,
            SiderustFrame::ICRS,
            J2000,
            &mut out,
        );
        assert_eq!(s, SiderustStatus::Ok);
    }

    #[test]
    fn equatorial_mean_of_date_to_icrs() {
        let mut out = empty_dir();
        let s = siderust_spherical_dir_transform_frame(
            10.0,
            90.0,
            SiderustFrame::EquatorialMeanOfDate,
            SiderustFrame::ICRS,
            J2000,
            &mut out,
        );
        assert_eq!(s, SiderustStatus::Ok);
    }

    #[test]
    fn equatorial_true_of_date_to_icrs() {
        let mut out = empty_dir();
        let s = siderust_spherical_dir_transform_frame(
            10.0,
            90.0,
            SiderustFrame::EquatorialTrueOfDate,
            SiderustFrame::ICRS,
            J2000,
            &mut out,
        );
        assert_eq!(s, SiderustStatus::Ok);
    }

    #[test]
    fn unsupported_src_frame_returns_error() {
        let mut out = empty_dir();
        let s = siderust_spherical_dir_transform_frame(
            0.0,
            0.0,
            SiderustFrame::Horizontal,
            SiderustFrame::ICRS,
            J2000,
            &mut out,
        );
        assert_eq!(s, SiderustStatus::InvalidFrame);
    }

    #[test]
    fn horizontal_dst_frame_returns_error() {
        let mut out = empty_dir();
        let s = siderust_spherical_dir_transform_frame(
            0.0,
            0.0,
            SiderustFrame::ICRS,
            SiderustFrame::Horizontal,
            J2000,
            &mut out,
        );
        assert_eq!(s, SiderustStatus::InvalidFrame);
    }

    #[test]
    fn galactic_dst_frame_returns_error() {
        let mut out = empty_dir();
        let s = siderust_spherical_dir_transform_frame(
            0.0,
            0.0,
            SiderustFrame::ICRS,
            SiderustFrame::Galactic,
            J2000,
            &mut out,
        );
        assert_eq!(s, SiderustStatus::InvalidFrame);
    }

    #[test]
    fn transform_frame_null_out() {
        let s = siderust_spherical_dir_transform_frame(
            0.0,
            0.0,
            SiderustFrame::ICRS,
            SiderustFrame::ICRS,
            J2000,
            ptr::null_mut(),
        );
        assert_eq!(s, SiderustStatus::NullPointer);
    }

    // ── To horizontal ─────────────────────────────────────────────────────

    #[test]
    fn icrs_to_horizontal_vega() {
        let mut out = empty_dir();
        // Vega: RA~279.2°, Dec~38.8°
        let s = siderust_spherical_dir_to_horizontal(
            38.8,
            279.2,
            SiderustFrame::ICRS,
            J2000,
            paris_observer(),
            &mut out,
        );
        assert_eq!(s, SiderustStatus::Ok);
        assert_eq!(out.frame, SiderustFrame::Horizontal);
        // alt/az should be finite
        assert!(out.polar_deg.is_finite());
        assert!(out.azimuth_deg.is_finite());
    }

    #[test]
    fn ecliptic_j2000_to_horizontal() {
        let mut out = empty_dir();
        let s = siderust_spherical_dir_to_horizontal(
            0.0,
            0.0,
            SiderustFrame::EclipticMeanJ2000,
            J2000,
            paris_observer(),
            &mut out,
        );
        assert_eq!(s, SiderustStatus::Ok);
    }

    #[test]
    fn to_horizontal_unsupported_src_frame() {
        let mut out = empty_dir();
        let s = siderust_spherical_dir_to_horizontal(
            0.0,
            0.0,
            SiderustFrame::ECEF,
            J2000,
            paris_observer(),
            &mut out,
        );
        assert_eq!(s, SiderustStatus::InvalidFrame);
    }

    #[test]
    fn to_horizontal_null_out() {
        let s = siderust_spherical_dir_to_horizontal(
            0.0,
            0.0,
            SiderustFrame::ICRS,
            J2000,
            paris_observer(),
            ptr::null_mut(),
        );
        assert_eq!(s, SiderustStatus::NullPointer);
    }

    // ── Geodetic to ECEF ─────────────────────────────────────────────────

    #[test]
    fn geodetic_to_ecef_greenwich() {
        let geodetic = SiderustGeodetict {
            lon_deg: 0.0,
            lat_deg: 51.5,
            height_m: 0.0,
        };
        let mut out = empty_cart();
        let s = siderust_geodetic_to_cartesian_ecef(geodetic, &mut out);
        assert_eq!(s, SiderustStatus::Ok);
        // ECEF position should be ~6350 km from origin
        let r = (out.x * out.x + out.y * out.y + out.z * out.z).sqrt();
        assert!(r > 6_300_000.0 && r < 6_400_000.0, "ECEF radius {r} m");
        assert_eq!(out.frame, SiderustFrame::ECEF);
        assert_eq!(out.center, SiderustCenter::Geocentric);
    }

    #[test]
    fn geodetic_to_ecef_null_out() {
        let geodetic = SiderustGeodetict {
            lon_deg: 0.0,
            lat_deg: 0.0,
            height_m: 0.0,
        };
        let s = siderust_geodetic_to_cartesian_ecef(geodetic, ptr::null_mut());
        assert_eq!(s, SiderustStatus::NullPointer);
    }

    // ── equatorial_mean_of_date and true_of_date as source ───────────────

    #[test]
    fn equatorial_mean_of_date_to_horizontal() {
        let mut out = empty_dir();
        let s = siderust_spherical_dir_to_horizontal(
            0.0,
            90.0,
            SiderustFrame::EquatorialMeanOfDate,
            J2000,
            paris_observer(),
            &mut out,
        );
        assert_eq!(s, SiderustStatus::Ok);
    }

    #[test]
    fn equatorial_true_of_date_to_horizontal() {
        let mut out = empty_dir();
        let s = siderust_spherical_dir_to_horizontal(
            0.0,
            90.0,
            SiderustFrame::EquatorialTrueOfDate,
            J2000,
            paris_observer(),
            &mut out,
        );
        assert_eq!(s, SiderustStatus::Ok);
    }

    // ── Cartesian direction frame transforms ─────────────────────────────

    fn unit_x_cart() -> (f64, f64, f64) {
        (1.0, 0.0, 0.0)
    }

    #[test]
    fn cartesian_dir_icrs_to_icrs_identity() {
        let mut out = empty_cart();
        let (x, y, z) = unit_x_cart();
        let s = siderust_cartesian_dir_transform_frame(
            x,
            y,
            z,
            SiderustFrame::ICRS,
            SiderustFrame::ICRS,
            J2000,
            &mut out,
        );
        assert_eq!(s, SiderustStatus::Ok);
        assert!((out.x - 1.0).abs() < 1e-8);
        assert!(out.y.abs() < 1e-8);
        assert!(out.z.abs() < 1e-8);
        assert_eq!(out.frame, SiderustFrame::ICRS);
    }

    #[test]
    fn cartesian_dir_icrs_to_ecliptic_j2000() {
        let mut out = empty_cart();
        let (x, y, z) = unit_x_cart();
        let s = siderust_cartesian_dir_transform_frame(
            x,
            y,
            z,
            SiderustFrame::ICRS,
            SiderustFrame::EclipticMeanJ2000,
            J2000,
            &mut out,
        );
        assert_eq!(s, SiderustStatus::Ok);
        assert_eq!(out.frame, SiderustFrame::EclipticMeanJ2000);
        // magnitude is preserved
        let mag = (out.x * out.x + out.y * out.y + out.z * out.z).sqrt();
        assert!((mag - 1.0).abs() < 1e-8);
    }

    #[test]
    fn cartesian_dir_icrs_to_equatorial_j2000() {
        let mut out = empty_cart();
        let s = siderust_cartesian_dir_transform_frame(
            1.0,
            0.0,
            0.0,
            SiderustFrame::ICRS,
            SiderustFrame::EquatorialMeanJ2000,
            J2000,
            &mut out,
        );
        assert_eq!(s, SiderustStatus::Ok);
        assert_eq!(out.frame, SiderustFrame::EquatorialMeanJ2000);
    }

    #[test]
    fn cartesian_dir_icrs_to_equatorial_mean_of_date() {
        let mut out = empty_cart();
        let s = siderust_cartesian_dir_transform_frame(
            1.0,
            0.0,
            0.0,
            SiderustFrame::ICRS,
            SiderustFrame::EquatorialMeanOfDate,
            J2000,
            &mut out,
        );
        assert_eq!(s, SiderustStatus::Ok);
        assert_eq!(out.frame, SiderustFrame::EquatorialMeanOfDate);
    }

    #[test]
    fn cartesian_dir_icrs_to_equatorial_true_of_date() {
        let mut out = empty_cart();
        let s = siderust_cartesian_dir_transform_frame(
            1.0,
            0.0,
            0.0,
            SiderustFrame::ICRS,
            SiderustFrame::EquatorialTrueOfDate,
            J2000,
            &mut out,
        );
        assert_eq!(s, SiderustStatus::Ok);
        assert_eq!(out.frame, SiderustFrame::EquatorialTrueOfDate);
    }

    #[test]
    fn cartesian_dir_ecliptic_to_icrs() {
        let mut out = empty_cart();
        let s = siderust_cartesian_dir_transform_frame(
            1.0,
            0.0,
            0.0,
            SiderustFrame::EclipticMeanJ2000,
            SiderustFrame::ICRS,
            J2000,
            &mut out,
        );
        assert_eq!(s, SiderustStatus::Ok);
    }

    #[test]
    fn cartesian_dir_invalid_src_frame() {
        let mut out = empty_cart();
        let s = siderust_cartesian_dir_transform_frame(
            1.0,
            0.0,
            0.0,
            SiderustFrame::Horizontal,
            SiderustFrame::ICRS,
            J2000,
            &mut out,
        );
        assert_eq!(s, SiderustStatus::InvalidFrame);
    }

    #[test]
    fn cartesian_dir_invalid_dst_frame() {
        let mut out = empty_cart();
        let s = siderust_cartesian_dir_transform_frame(
            1.0,
            0.0,
            0.0,
            SiderustFrame::ICRS,
            SiderustFrame::Galactic,
            J2000,
            &mut out,
        );
        assert_eq!(s, SiderustStatus::InvalidFrame);
    }

    #[test]
    fn cartesian_dir_null_out() {
        let s = siderust_cartesian_dir_transform_frame(
            1.0,
            0.0,
            0.0,
            SiderustFrame::ICRS,
            SiderustFrame::ICRS,
            J2000,
            ptr::null_mut(),
        );
        assert_eq!(s, SiderustStatus::NullPointer);
    }

    // ── Cartesian position frame transforms ──────────────────────────────

    fn earth_pos_icrs() -> SiderustCartesianPos {
        // Approximate Earth position at J2000 (AU) in ICRS
        SiderustCartesianPos {
            x: -0.177,
            y: 0.969,
            z: 0.0,
            frame: SiderustFrame::ICRS,
            center: SiderustCenter::Barycentric,
        }
    }

    #[test]
    fn cartesian_pos_icrs_to_icrs_identity() {
        let pos = earth_pos_icrs();
        let mut out = empty_cart();
        let s = siderust_cartesian_pos_transform_frame(pos, SiderustFrame::ICRS, J2000, &mut out);
        assert_eq!(s, SiderustStatus::Ok);
        assert!((out.x - pos.x).abs() < 1e-8);
        assert!((out.y - pos.y).abs() < 1e-8);
        assert_eq!(out.frame, SiderustFrame::ICRS);
        assert_eq!(out.center, SiderustCenter::Barycentric);
    }

    #[test]
    fn cartesian_pos_icrs_to_ecliptic_j2000() {
        let pos = earth_pos_icrs();
        let mut out = empty_cart();
        let s = siderust_cartesian_pos_transform_frame(
            pos,
            SiderustFrame::EclipticMeanJ2000,
            J2000,
            &mut out,
        );
        assert_eq!(s, SiderustStatus::Ok);
        assert_eq!(out.frame, SiderustFrame::EclipticMeanJ2000);
        // magnitude preserved
        let mag_in = (pos.x * pos.x + pos.y * pos.y + pos.z * pos.z).sqrt();
        let mag_out = (out.x * out.x + out.y * out.y + out.z * out.z).sqrt();
        assert!((mag_in - mag_out).abs() < 1e-6);
    }

    #[test]
    fn cartesian_pos_icrs_to_equatorial_j2000() {
        let pos = earth_pos_icrs();
        let mut out = empty_cart();
        let s = siderust_cartesian_pos_transform_frame(
            pos,
            SiderustFrame::EquatorialMeanJ2000,
            J2000,
            &mut out,
        );
        assert_eq!(s, SiderustStatus::Ok);
        assert_eq!(out.frame, SiderustFrame::EquatorialMeanJ2000);
    }

    #[test]
    fn cartesian_pos_icrs_to_equatorial_mean_of_date() {
        let pos = earth_pos_icrs();
        let mut out = empty_cart();
        let s = siderust_cartesian_pos_transform_frame(
            pos,
            SiderustFrame::EquatorialMeanOfDate,
            J2000,
            &mut out,
        );
        assert_eq!(s, SiderustStatus::Ok);
        assert_eq!(out.frame, SiderustFrame::EquatorialMeanOfDate);
    }

    #[test]
    fn cartesian_pos_icrs_to_equatorial_true_of_date() {
        let pos = earth_pos_icrs();
        let mut out = empty_cart();
        let s = siderust_cartesian_pos_transform_frame(
            pos,
            SiderustFrame::EquatorialTrueOfDate,
            J2000,
            &mut out,
        );
        assert_eq!(s, SiderustStatus::Ok);
        assert_eq!(out.frame, SiderustFrame::EquatorialTrueOfDate);
    }

    #[test]
    fn cartesian_pos_icrs_to_icrf() {
        let pos = earth_pos_icrs();
        let mut out = empty_cart();
        let s = siderust_cartesian_pos_transform_frame(pos, SiderustFrame::ICRF, J2000, &mut out);
        assert_eq!(s, SiderustStatus::Ok);
        assert_eq!(out.frame, SiderustFrame::ICRF);
    }

    #[test]
    fn cartesian_pos_ecliptic_to_icrs() {
        let pos = SiderustCartesianPos {
            x: 1.0,
            y: 0.0,
            z: 0.0,
            frame: SiderustFrame::EclipticMeanJ2000,
            center: SiderustCenter::Heliocentric,
        };
        let mut out = empty_cart();
        let s = siderust_cartesian_pos_transform_frame(pos, SiderustFrame::ICRS, J2000, &mut out);
        assert_eq!(s, SiderustStatus::Ok);
    }

    #[test]
    fn cartesian_pos_invalid_src_frame() {
        let pos = SiderustCartesianPos {
            x: 1.0,
            y: 0.0,
            z: 0.0,
            frame: SiderustFrame::Horizontal,
            center: SiderustCenter::Geocentric,
        };
        let mut out = empty_cart();
        let s = siderust_cartesian_pos_transform_frame(pos, SiderustFrame::ICRS, J2000, &mut out);
        assert_eq!(s, SiderustStatus::InvalidFrame);
    }

    #[test]
    fn cartesian_pos_invalid_dst_frame() {
        let pos = earth_pos_icrs();
        let mut out = empty_cart();
        let s =
            siderust_cartesian_pos_transform_frame(pos, SiderustFrame::Galactic, J2000, &mut out);
        assert_eq!(s, SiderustStatus::InvalidFrame);
    }

    #[test]
    fn cartesian_pos_null_out() {
        let pos = earth_pos_icrs();
        let s = siderust_cartesian_pos_transform_frame(
            pos,
            SiderustFrame::ICRS,
            J2000,
            ptr::null_mut(),
        );
        assert_eq!(s, SiderustStatus::NullPointer);
    }

    // ── Cartesian position center transforms ─────────────────────────────

    #[test]
    fn center_transform_same_center_is_identity() {
        let pos = SiderustCartesianPos {
            x: 1.0,
            y: 0.5,
            z: -0.3,
            frame: SiderustFrame::EclipticMeanJ2000,
            center: SiderustCenter::Heliocentric,
        };
        let mut out = empty_cart();
        let s = siderust_cartesian_pos_transform_center(
            pos,
            SiderustCenter::Heliocentric,
            J2000,
            &mut out,
        );
        assert_eq!(s, SiderustStatus::Ok);
        assert!((out.x - 1.0).abs() < 1e-8);
        assert!((out.y - 0.5).abs() < 1e-8);
        assert!((out.z + 0.3).abs() < 1e-8);
        assert_eq!(out.center, SiderustCenter::Heliocentric);
    }

    #[test]
    fn center_transform_helio_to_bary() {
        let pos = SiderustCartesianPos {
            x: 1.0,
            y: 0.0,
            z: 0.0,
            frame: SiderustFrame::EclipticMeanJ2000,
            center: SiderustCenter::Heliocentric,
        };
        let mut out = empty_cart();
        let s = siderust_cartesian_pos_transform_center(
            pos,
            SiderustCenter::Barycentric,
            J2000,
            &mut out,
        );
        assert_eq!(s, SiderustStatus::Ok);
        assert_eq!(out.center, SiderustCenter::Barycentric);
        assert!(out.x.is_finite() && out.y.is_finite() && out.z.is_finite());
    }

    #[test]
    fn center_transform_helio_to_geo() {
        let pos = SiderustCartesianPos {
            x: 1.0,
            y: 0.0,
            z: 0.0,
            frame: SiderustFrame::EclipticMeanJ2000,
            center: SiderustCenter::Heliocentric,
        };
        let mut out = empty_cart();
        let s = siderust_cartesian_pos_transform_center(
            pos,
            SiderustCenter::Geocentric,
            J2000,
            &mut out,
        );
        assert_eq!(s, SiderustStatus::Ok);
        assert_eq!(out.center, SiderustCenter::Geocentric);
        assert!(out.x.is_finite());
    }

    #[test]
    fn center_transform_invalid_src_center() {
        let pos = SiderustCartesianPos {
            x: 1.0,
            y: 0.0,
            z: 0.0,
            frame: SiderustFrame::EclipticMeanJ2000,
            center: SiderustCenter::Bodycentric,
        };
        let mut out = empty_cart();
        let s = siderust_cartesian_pos_transform_center(
            pos,
            SiderustCenter::Heliocentric,
            J2000,
            &mut out,
        );
        assert_eq!(s, SiderustStatus::InvalidCenter);
    }

    #[test]
    fn center_transform_invalid_dst_center() {
        let pos = SiderustCartesianPos {
            x: 1.0,
            y: 0.0,
            z: 0.0,
            frame: SiderustFrame::EclipticMeanJ2000,
            center: SiderustCenter::Heliocentric,
        };
        let mut out = empty_cart();
        let s = siderust_cartesian_pos_transform_center(
            pos,
            SiderustCenter::Bodycentric,
            J2000,
            &mut out,
        );
        assert_eq!(s, SiderustStatus::InvalidCenter);
    }

    #[test]
    fn center_transform_null_out() {
        let pos = SiderustCartesianPos {
            x: 1.0,
            y: 0.0,
            z: 0.0,
            frame: SiderustFrame::EclipticMeanJ2000,
            center: SiderustCenter::Heliocentric,
        };
        let s = siderust_cartesian_pos_transform_center(
            pos,
            SiderustCenter::Barycentric,
            J2000,
            ptr::null_mut(),
        );
        assert_eq!(s, SiderustStatus::NullPointer);
    }

    // ── Kepler position ───────────────────────────────────────────────────

    // Earth's approximate J2000 orbital elements
    fn earth_orbit() -> SiderustOrbit {
        SiderustOrbit {
            semi_major_axis_au: 1.0000,
            eccentricity: 0.0167,
            inclination_deg: 0.0,
            lon_ascending_node_deg: 0.0,
            arg_perihelion_deg: 102.94,
            mean_anomaly_deg: 357.53,
            epoch_jd: J2000,
        }
    }

    #[test]
    fn kepler_position_earth_at_j2000() {
        let mut out = empty_cart();
        let s = siderust_kepler_position(earth_orbit(), J2000, &mut out);
        assert_eq!(s, SiderustStatus::Ok);
        // Earth should be ~1 AU from the Sun
        let r = (out.x * out.x + out.y * out.y + out.z * out.z).sqrt();
        assert!(r > 0.98 && r < 1.02, "r = {r}");
        assert_eq!(out.frame, SiderustFrame::EclipticMeanJ2000);
    }

    #[test]
    fn kepler_position_null_out() {
        let s = siderust_kepler_position(earth_orbit(), J2000, ptr::null_mut());
        assert_eq!(s, SiderustStatus::NullPointer);
    }

    // ── Bodycentric transforms ─────────────────────────────────────────────

    fn earth_bodycentric_params() -> SiderustBodycentricParams {
        SiderustBodycentricParams {
            orbit: earth_orbit(),
            orbit_center: 1, // Heliocentric
            _pad: [0u8; 7],
        }
    }

    #[test]
    fn to_bodycentric_and_back_roundtrip() {
        // A geocentric position
        let pos = SiderustCartesianPos {
            x: 0.0,
            y: 0.0,
            z: 0.0,
            frame: SiderustFrame::EclipticMeanJ2000,
            center: SiderustCenter::Geocentric,
        };
        let params = earth_bodycentric_params();
        let mut body_pos = empty_cart();
        let s1 = siderust_to_bodycentric(pos, params, J2000, &mut body_pos);
        assert_eq!(s1, SiderustStatus::Ok);
        assert_eq!(body_pos.center, SiderustCenter::Bodycentric);

        let mut recovered = empty_cart();
        let s2 = siderust_from_bodycentric(body_pos, params, J2000, &mut recovered);
        assert_eq!(s2, SiderustStatus::Ok);
        assert_eq!(recovered.center, SiderustCenter::Geocentric);
        // Round-trip: recovered should be close to original pos (geocentric origin)
        let diff = ((recovered.x - pos.x).powi(2)
            + (recovered.y - pos.y).powi(2)
            + (recovered.z - pos.z).powi(2))
        .sqrt();
        assert!(diff < 1e-6, "round-trip error = {diff}");
    }

    #[test]
    fn to_bodycentric_invalid_center() {
        let pos = SiderustCartesianPos {
            x: 1.0,
            y: 0.0,
            z: 0.0,
            frame: SiderustFrame::EclipticMeanJ2000,
            center: SiderustCenter::Bodycentric, // invalid for input
        };
        let mut out = empty_cart();
        let s = siderust_to_bodycentric(pos, earth_bodycentric_params(), J2000, &mut out);
        assert_eq!(s, SiderustStatus::InvalidCenter);
    }

    #[test]
    fn to_bodycentric_null_out() {
        let pos = SiderustCartesianPos {
            x: 1.0,
            y: 0.0,
            z: 0.0,
            frame: SiderustFrame::EclipticMeanJ2000,
            center: SiderustCenter::Geocentric,
        };
        let s = siderust_to_bodycentric(pos, earth_bodycentric_params(), J2000, ptr::null_mut());
        assert_eq!(s, SiderustStatus::NullPointer);
    }

    #[test]
    fn from_bodycentric_null_out() {
        let pos = SiderustCartesianPos {
            x: 0.1,
            y: 0.0,
            z: 0.0,
            frame: SiderustFrame::EclipticMeanJ2000,
            center: SiderustCenter::Bodycentric,
        };
        let s = siderust_from_bodycentric(pos, earth_bodycentric_params(), J2000, ptr::null_mut());
        assert_eq!(s, SiderustStatus::NullPointer);
    }

    #[test]
    fn to_bodycentric_heliocentric_input() {
        let pos = SiderustCartesianPos {
            x: 1.0,
            y: 0.0,
            z: 0.0,
            frame: SiderustFrame::EclipticMeanJ2000,
            center: SiderustCenter::Heliocentric,
        };
        let mut out = empty_cart();
        let s = siderust_to_bodycentric(pos, earth_bodycentric_params(), J2000, &mut out);
        assert_eq!(s, SiderustStatus::Ok);
        assert_eq!(out.center, SiderustCenter::Bodycentric);
    }

    #[test]
    fn to_bodycentric_barycentric_input() {
        let pos = SiderustCartesianPos {
            x: 1.0,
            y: 0.0,
            z: 0.0,
            frame: SiderustFrame::EclipticMeanJ2000,
            center: SiderustCenter::Barycentric,
        };
        let mut out = empty_cart();
        let s = siderust_to_bodycentric(pos, earth_bodycentric_params(), J2000, &mut out);
        assert_eq!(s, SiderustStatus::Ok);
        assert_eq!(out.center, SiderustCenter::Bodycentric);
    }
}
