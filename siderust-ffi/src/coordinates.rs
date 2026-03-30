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
use siderust::coordinates::centers::Barycentric;
use siderust::coordinates::frames::{
    EclipticMeanJ2000, EquatorialMeanJ2000, EquatorialMeanOfDate, EquatorialTrueOfDate,
    ReferenceFrame, ICRF, ICRS,
};
use siderust::coordinates::transform::{
    DirectionAstroExt, PositionAstroExt, SphericalDirectionAstroExt,
};
use siderust::coordinates::{cartesian, spherical};
use siderust::time::JulianDate;

// ═══════════════════════════════════════════════════════════════════════════
// Spherical Direction, frame transforms
// ═══════════════════════════════════════════════════════════════════════════

/// Helper: extract (polar_deg, azimuth_deg) from a spherical Direction.
#[inline]
fn sph_to_pair<F: ReferenceFrame>(d: &spherical::Direction<F>) -> (f64, f64) {
    (d.polar.value(), d.azimuth.value())
}

macro_rules! transform_spherical_dir_from {
    ($dir:expr, $dst_frame:expr, $jd:expr) => {{
        match $dst_frame {
            SiderustFrame::ICRS => Ok(sph_to_pair(&SphericalDirectionAstroExt::to_frame::<ICRS>(
                &$dir, $jd,
            ))),
            SiderustFrame::EclipticMeanJ2000 => {
                let icrs = SphericalDirectionAstroExt::to_frame::<ICRS>(&$dir, $jd);
                Ok(sph_to_pair(&SphericalDirectionAstroExt::to_frame::<
                    EclipticMeanJ2000,
                >(&icrs, $jd)))
            }
            SiderustFrame::EquatorialMeanJ2000 => {
                let icrs = SphericalDirectionAstroExt::to_frame::<ICRS>(&$dir, $jd);
                Ok(sph_to_pair(&SphericalDirectionAstroExt::to_frame::<
                    EquatorialMeanJ2000,
                >(&icrs, $jd)))
            }
            SiderustFrame::EquatorialMeanOfDate => {
                let icrs = SphericalDirectionAstroExt::to_frame::<ICRS>(&$dir, $jd);
                Ok(sph_to_pair(&SphericalDirectionAstroExt::to_frame::<
                    EquatorialMeanOfDate,
                >(&icrs, $jd)))
            }
            SiderustFrame::EquatorialTrueOfDate => {
                let icrs = SphericalDirectionAstroExt::to_frame::<ICRS>(&$dir, $jd);
                Ok(sph_to_pair(&SphericalDirectionAstroExt::to_frame::<
                    EquatorialTrueOfDate,
                >(&icrs, $jd)))
            }
            _ => Err(SiderustStatus::InvalidFrame),
        }
    }};
}

macro_rules! spherical_dir_to_horizontal_from {
    ($dir:expr, $jd:expr, $site:expr) => {{
        let icrs = SphericalDirectionAstroExt::to_frame::<ICRS>(&$dir, $jd);
        let eq_tod = SphericalDirectionAstroExt::to_frame::<EquatorialTrueOfDate>(&icrs, $jd);
        let cart = eq_tod.to_cartesian();
        let hz = DirectionAstroExt::to_horizontal(&cart, $jd, $site);
        let hz_sph = spherical::Direction::from_cartesian(&hz);
        sph_to_pair(&hz_sph)
    }};
}

macro_rules! with_spherical_dir_in_frame {
    ($frame:expr, $polar_deg:expr, $azimuth_deg:expr, |$dir:ident| $body:expr) => {{
        match $frame {
            SiderustFrame::ICRS => {
                let $dir = spherical::Direction::<ICRS>::new(
                    Degrees::new($azimuth_deg),
                    Degrees::new($polar_deg),
                );
                Ok($body)
            }
            SiderustFrame::EclipticMeanJ2000 => {
                let $dir = spherical::Direction::<EclipticMeanJ2000>::new(
                    Degrees::new($azimuth_deg),
                    Degrees::new($polar_deg),
                );
                Ok($body)
            }
            SiderustFrame::EquatorialMeanJ2000 => {
                let $dir = spherical::Direction::<EquatorialMeanJ2000>::new(
                    Degrees::new($azimuth_deg),
                    Degrees::new($polar_deg),
                );
                Ok($body)
            }
            SiderustFrame::EquatorialMeanOfDate => {
                let $dir = spherical::Direction::<EquatorialMeanOfDate>::new(
                    Degrees::new($azimuth_deg),
                    Degrees::new($polar_deg),
                );
                Ok($body)
            }
            SiderustFrame::EquatorialTrueOfDate => {
                let $dir = spherical::Direction::<EquatorialTrueOfDate>::new(
                    Degrees::new($azimuth_deg),
                    Degrees::new($polar_deg),
                );
                Ok($body)
            }
            _ => Err(SiderustStatus::InvalidFrame),
        }
    }};
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

        let (out_polar, out_azimuth) = match with_spherical_dir_in_frame!(
            src_frame,
            polar_deg,
            azimuth_deg,
            |dir| transform_spherical_dir_from!(dir, dst_frame, &JulianDate::new(jd))
        ) {
            Ok(Ok(values)) => values,
            Err(e) => return e,
            Ok(Err(e)) => return e,
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

        let site = observer.to_rust();
        let (az, alt) = match with_spherical_dir_in_frame!(
            src_frame,
            polar_deg,
            azimuth_deg,
            |dir| spherical_dir_to_horizontal_from!(dir, &JulianDate::new(jd), &site)
        ) {
            Ok(values) => values,
            Err(e) => return e,
        };

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
// Cartesian Direction, frame transforms
// ═══════════════════════════════════════════════════════════════════════════

#[inline]
fn cart_dir_to_triple<F: ReferenceFrame>(d: &cartesian::Direction<F>) -> (f64, f64, f64) {
    (d.x(), d.y(), d.z())
}

macro_rules! transform_cartesian_dir_from {
    ($dir:expr, $dst_frame:expr, $jd:expr) => {{
        match $dst_frame {
            SiderustFrame::ICRS => Ok(cart_dir_to_triple(&DirectionAstroExt::to_frame::<ICRS>(
                &$dir, $jd,
            ))),
            SiderustFrame::EclipticMeanJ2000 => {
                let icrs = DirectionAstroExt::to_frame::<ICRS>(&$dir, $jd);
                Ok(cart_dir_to_triple(&DirectionAstroExt::to_frame::<
                    EclipticMeanJ2000,
                >(&icrs, $jd)))
            }
            SiderustFrame::EquatorialMeanJ2000 => {
                let icrs = DirectionAstroExt::to_frame::<ICRS>(&$dir, $jd);
                Ok(cart_dir_to_triple(&DirectionAstroExt::to_frame::<
                    EquatorialMeanJ2000,
                >(&icrs, $jd)))
            }
            SiderustFrame::EquatorialMeanOfDate => {
                let icrs = DirectionAstroExt::to_frame::<ICRS>(&$dir, $jd);
                Ok(cart_dir_to_triple(&DirectionAstroExt::to_frame::<
                    EquatorialMeanOfDate,
                >(&icrs, $jd)))
            }
            SiderustFrame::EquatorialTrueOfDate => {
                let icrs = DirectionAstroExt::to_frame::<ICRS>(&$dir, $jd);
                Ok(cart_dir_to_triple(&DirectionAstroExt::to_frame::<
                    EquatorialTrueOfDate,
                >(&icrs, $jd)))
            }
            _ => Err(SiderustStatus::InvalidFrame),
        }
    }};
}

macro_rules! with_cartesian_dir_in_frame {
    ($frame:expr, $x:expr, $y:expr, $z:expr, |$dir:ident| $body:expr) => {{
        match $frame {
            SiderustFrame::ICRS => {
                let $dir = cartesian::Direction::<ICRS>::new_unchecked($x, $y, $z);
                Ok($body)
            }
            SiderustFrame::EclipticMeanJ2000 => {
                let $dir = cartesian::Direction::<EclipticMeanJ2000>::new_unchecked($x, $y, $z);
                Ok($body)
            }
            SiderustFrame::EquatorialMeanJ2000 => {
                let $dir = cartesian::Direction::<EquatorialMeanJ2000>::new_unchecked($x, $y, $z);
                Ok($body)
            }
            SiderustFrame::EquatorialMeanOfDate => {
                let $dir = cartesian::Direction::<EquatorialMeanOfDate>::new_unchecked($x, $y, $z);
                Ok($body)
            }
            SiderustFrame::EquatorialTrueOfDate => {
                let $dir = cartesian::Direction::<EquatorialTrueOfDate>::new_unchecked($x, $y, $z);
                Ok($body)
            }
            _ => Err(SiderustStatus::InvalidFrame),
        }
    }};
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

        let (ox, oy, oz) = match with_cartesian_dir_in_frame!(
            src_frame,
            x,
            y,
            z,
            |dir| transform_cartesian_dir_from!(dir, dst_frame, &JulianDate::new(jd))
        ) {
            Ok(Ok(values)) => values,
            Err(e) => return e,
            Ok(Err(e)) => return e,
        };

        unsafe {
            *out = SiderustCartesianPos {
                x: ox,
                y: oy,
                z: oz,
                frame: dst_frame,
                center: SiderustCenter::Barycentric, // directions have no center
                length_unit: SiderustLengthUnit::AU,
            };
        }
        SiderustStatus::Ok

    }}
}

// ═══════════════════════════════════════════════════════════════════════════
// Cartesian Position, frame transforms
// ═══════════════════════════════════════════════════════════════════════════

/// Helper: extract raw (x, y, z) from a Cartesian position.
#[inline]
fn cart_pos_to_triple<F: siderust::coordinates::frames::ReferenceFrame>(
    p: &cartesian::Position<Barycentric, F, AstronomicalUnit>,
) -> (f64, f64, f64) {
    (p.x().value(), p.y().value(), p.z().value())
}

macro_rules! transform_cartesian_pos_from {
    ($pos:expr, $dst_frame:expr, $jd:expr) => {{
        match $dst_frame {
            SiderustFrame::ICRS => Ok(cart_pos_to_triple(&PositionAstroExt::to_frame::<ICRS>(
                &$pos, $jd,
            ))),
            SiderustFrame::ICRF => {
                let icrs: cartesian::Position<Barycentric, ICRS, AstronomicalUnit> =
                    PositionAstroExt::to_frame::<ICRS>(&$pos, $jd);
                Ok(cart_pos_to_triple(&PositionAstroExt::to_frame::<ICRF>(
                    &icrs, $jd,
                )))
            }
            SiderustFrame::EclipticMeanJ2000 => {
                let icrs: cartesian::Position<Barycentric, ICRS, AstronomicalUnit> =
                    PositionAstroExt::to_frame::<ICRS>(&$pos, $jd);
                Ok(cart_pos_to_triple(&PositionAstroExt::to_frame::<
                    EclipticMeanJ2000,
                >(&icrs, $jd)))
            }
            SiderustFrame::EquatorialMeanJ2000 => {
                let icrs: cartesian::Position<Barycentric, ICRS, AstronomicalUnit> =
                    PositionAstroExt::to_frame::<ICRS>(&$pos, $jd);
                Ok(cart_pos_to_triple(&PositionAstroExt::to_frame::<
                    EquatorialMeanJ2000,
                >(&icrs, $jd)))
            }
            SiderustFrame::EquatorialMeanOfDate => {
                let icrs: cartesian::Position<Barycentric, ICRS, AstronomicalUnit> =
                    PositionAstroExt::to_frame::<ICRS>(&$pos, $jd);
                Ok(cart_pos_to_triple(&PositionAstroExt::to_frame::<
                    EquatorialMeanOfDate,
                >(&icrs, $jd)))
            }
            SiderustFrame::EquatorialTrueOfDate => {
                let icrs: cartesian::Position<Barycentric, ICRS, AstronomicalUnit> =
                    PositionAstroExt::to_frame::<ICRS>(&$pos, $jd);
                Ok(cart_pos_to_triple(&PositionAstroExt::to_frame::<
                    EquatorialTrueOfDate,
                >(&icrs, $jd)))
            }
            _ => Err(SiderustStatus::InvalidFrame),
        }
    }};
}

macro_rules! with_cartesian_pos_in_frame {
    ($frame:expr, $x:expr, $y:expr, $z:expr, |$pos:ident| $body:expr) => {{
        match $frame {
            SiderustFrame::ICRS => {
                let $pos =
                    cartesian::Position::<Barycentric, ICRS, AstronomicalUnit>::new($x, $y, $z);
                Ok($body)
            }
            SiderustFrame::ICRF => {
                let $pos =
                    cartesian::Position::<Barycentric, ICRF, AstronomicalUnit>::new($x, $y, $z);
                Ok($body)
            }
            SiderustFrame::EclipticMeanJ2000 => {
                let $pos =
                    cartesian::Position::<Barycentric, EclipticMeanJ2000, AstronomicalUnit>::new(
                        $x, $y, $z,
                    );
                Ok($body)
            }
            SiderustFrame::EquatorialMeanJ2000 => {
                let $pos =
                    cartesian::Position::<Barycentric, EquatorialMeanJ2000, AstronomicalUnit>::new(
                        $x, $y, $z,
                    );
                Ok($body)
            }
            SiderustFrame::EquatorialMeanOfDate => {
                let $pos =
                    cartesian::Position::<Barycentric, EquatorialMeanOfDate, AstronomicalUnit>::new(
                        $x, $y, $z,
                    );
                Ok($body)
            }
            SiderustFrame::EquatorialTrueOfDate => {
                let $pos =
                    cartesian::Position::<Barycentric, EquatorialTrueOfDate, AstronomicalUnit>::new(
                        $x, $y, $z,
                    );
                Ok($body)
            }
            _ => Err(SiderustStatus::InvalidFrame),
        }
    }};
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

        let (ox, oy, oz) = match with_cartesian_pos_in_frame!(
            pos.frame,
            pos.x,
            pos.y,
            pos.z,
            |source| transform_cartesian_pos_from!(source, dst_frame, &JulianDate::new(jd))
        ) {
            Ok(Ok(values)) => values,
            Err(e) => return e,
            Ok(Err(e)) => return e,
        };

        unsafe {
            *out = SiderustCartesianPos {
                x: ox,
                y: oy,
                z: oz,
                frame: dst_frame,
                center: pos.center, // center is unchanged for a frame-only transform
                length_unit: pos.length_unit,
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
                length_unit: SiderustLengthUnit::Meter,
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

    // Sun's barycentric position, used for helio ↔ bary shifts
    let sun_b = Vsop87Ephemeris::sun_barycentric(t);
    let (sb_x, sb_y, sb_z) = (sun_b.x().value(), sun_b.y().value(), sun_b.z().value());

    // Earth's heliocentric position, used for helio ↔ geo shifts
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
                length_unit: pos.length_unit,
            };
        }
        SiderustStatus::Ok

    }}
}

// Shared serializer for orbit propagation functions that return
// EclipticMeanJ2000 positions in astronomical units.
fn write_ecliptic_au_position(
    out: *mut SiderustCartesianPos,
    pos: siderust::coordinates::cartesian::position::EclipticMeanJ2000<AstronomicalUnit>,
    center: SiderustCenter,
) {
    unsafe {
        *out = SiderustCartesianPos {
            x: pos.x().value(),
            y: pos.y().value(),
            z: pos.z().value(),
            frame: SiderustFrame::EclipticMeanJ2000,
            center,
            length_unit: SiderustLengthUnit::AU,
        };
    }
}

/// Map an `SiderustOrbitRefCenter` (0=Bary,1=Helio,2=Geo) to a `SiderustCenter`.
fn orbit_ref_center_to_siderust(orbit_center: SiderustOrbitRefCenter) -> Option<SiderustCenter> {
    match orbit_center {
        0 => Some(SiderustCenter::Barycentric),
        1 => Some(SiderustCenter::Heliocentric),
        2 => Some(SiderustCenter::Geocentric),
        _ => None,
    }
}

/// Compute the Keplerian orbital position at a given Julian date, with
/// explicit output center semantics.
///
/// `orbit_center` specifies which reference center the orbit elements are
/// relative to (0=Barycentric, 1=Heliocentric, 2=Geocentric). The `center`
/// field of `out` will be set accordingly.
#[no_mangle]
pub extern "C" fn siderust_kepler_position_ex(
    orbit: SiderustOrbit,
    orbit_center: SiderustOrbitRefCenter,
    jd: f64,
    out: *mut SiderustCartesianPos,
) -> SiderustStatus {
    ffi_guard! {{
        if out.is_null() {
            return SiderustStatus::NullPointer;
        }
        let center = match orbit_ref_center_to_siderust(orbit_center) {
            Some(c) => c,
            None => return SiderustStatus::InvalidCenter,
        };
        let pos = orbit.to_rust().kepler_position(JulianDate::new(jd));
        write_ecliptic_au_position(out, pos, center);
        SiderustStatus::Ok

    }}
}

/// Compute the position of an explicit mean-motion orbit at a given Julian date.
#[no_mangle]
pub extern "C" fn siderust_mean_motion_position(
    orbit: SiderustMeanMotionOrbit,
    jd: f64,
    out: *mut SiderustCartesianPos,
) -> SiderustStatus {
    ffi_guard! {{
        if out.is_null() {
            return SiderustStatus::NullPointer;
        }
        let rust_orbit = match orbit.try_to_rust() {
            Ok(o) => o,
            Err(_) => return SiderustStatus::InvalidArgument,
        };
        let pos = match rust_orbit.position_at(JulianDate::new(jd)) {
            Ok(pos) => pos,
            Err(_) => return SiderustStatus::InvalidArgument,
        };
        write_ecliptic_au_position(out, pos, SiderustCenter::Heliocentric);
        SiderustStatus::Ok

    }}
}

/// Compute the position of a unified conic orbit at a given Julian date.
#[no_mangle]
pub extern "C" fn siderust_conic_position(
    orbit: SiderustConicOrbit,
    jd: f64,
    out: *mut SiderustCartesianPos,
) -> SiderustStatus {
    ffi_guard! {{
        if out.is_null() {
            return SiderustStatus::NullPointer;
        }
        let rust_orbit = match orbit.try_to_rust() {
            Ok(o) => o,
            Err(_) => return SiderustStatus::InvalidArgument,
        };
        let pos = match rust_orbit.position_at(JulianDate::new(jd)) {
            Ok(pos) => pos,
            Err(_) => return SiderustStatus::InvalidArgument,
        };
        write_ecliptic_au_position(out, pos, SiderustCenter::Heliocentric);
        SiderustStatus::Ok

    }}
}

// =============================================================================
// Opaque PreparedOrbit handle for amortized propagation
// =============================================================================

/// Opaque handle to a validated, precomputed Keplerian orbit.
///
/// Created via [`siderust_prepared_orbit_create`], used for repeated
/// propagation via [`siderust_prepared_orbit_position`], and freed via
/// [`siderust_prepared_orbit_destroy`].
///
/// This avoids per-call reconstruction, validation, and trig precomputation
/// when propagating the same orbit at many epochs.
pub type SiderustPreparedOrbitHandle = *mut std::ffi::c_void;

/// Create a validated, precomputed orbit handle from Keplerian elements.
///
/// Returns `Ok` and writes a non-null handle on success. Returns
/// `InvalidArgument` if the orbit is not a valid ellipse.
#[no_mangle]
pub extern "C" fn siderust_prepared_orbit_create(
    orbit: SiderustOrbit,
    out: *mut SiderustPreparedOrbitHandle,
) -> SiderustStatus {
    ffi_guard! {{
        if out.is_null() {
            return SiderustStatus::NullPointer;
        }
        let prepared = match siderust::PreparedOrbit::try_from_elements(
            AstronomicalUnits::new(orbit.semi_major_axis_au),
            orbit.eccentricity,
            Degrees::new(orbit.inclination_deg),
            Degrees::new(orbit.lon_ascending_node_deg),
            Degrees::new(orbit.arg_periapsis_deg),
            Degrees::new(orbit.mean_anomaly_deg),
            JulianDate::new(orbit.epoch_jd),
        ) {
            Ok(p) => p,
            Err(_) => return SiderustStatus::InvalidArgument,
        };
        let boxed = Box::new(prepared);
        unsafe { *out = Box::into_raw(boxed) as *mut std::ffi::c_void };
        SiderustStatus::Ok
    }}
}

/// Propagate a prepared orbit to a given Julian date.
///
/// `handle` must have been created via [`siderust_prepared_orbit_create`].
/// This is the fastest FFI propagation path: no validation, no reconstruction.
#[no_mangle]
pub extern "C" fn siderust_prepared_orbit_position(
    handle: SiderustPreparedOrbitHandle,
    jd: f64,
    out: *mut SiderustCartesianPos,
) -> SiderustStatus {
    ffi_guard! {{
        if handle.is_null() || out.is_null() {
            return SiderustStatus::NullPointer;
        }
        let prepared = unsafe { &*(handle as *const siderust::PreparedOrbit) };
        let pos = prepared.position_at(JulianDate::new(jd));
        write_ecliptic_au_position(out, pos, SiderustCenter::Heliocentric);
        SiderustStatus::Ok
    }}
}

/// Destroy a prepared orbit handle, freeing its memory.
///
/// After this call, `handle` must not be used again.
#[no_mangle]
pub extern "C" fn siderust_prepared_orbit_destroy(
    handle: SiderustPreparedOrbitHandle,
) -> SiderustStatus {
    ffi_guard! {{
        if handle.is_null() {
            return SiderustStatus::NullPointer;
        }
        unsafe { drop(Box::from_raw(handle as *mut siderust::PreparedOrbit)) };
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
                length_unit: pos.length_unit,
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
            length_unit: SiderustLengthUnit::AU,
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
            length_unit: SiderustLengthUnit::AU,
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
            length_unit: SiderustLengthUnit::AU,
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
            length_unit: SiderustLengthUnit::AU,
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
            length_unit: SiderustLengthUnit::AU,
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
            length_unit: SiderustLengthUnit::AU,
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
            length_unit: SiderustLengthUnit::AU,
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
            length_unit: SiderustLengthUnit::AU,
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
            length_unit: SiderustLengthUnit::AU,
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
            length_unit: SiderustLengthUnit::AU,
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
            length_unit: SiderustLengthUnit::AU,
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
            arg_periapsis_deg: 102.94,
            mean_anomaly_deg: 357.53,
            epoch_jd: J2000,
        }
    }

    fn earth_mean_motion_orbit() -> SiderustMeanMotionOrbit {
        SiderustMeanMotionOrbit {
            semi_major_axis_au: 1.0000,
            eccentricity: 0.0167,
            inclination_deg: 0.0,
            lon_ascending_node_deg: 0.0,
            arg_periapsis_deg: 102.94,
            mean_motion_deg_per_day: 0.9856076686,
            epoch_jd: J2000,
        }
    }

    #[test]
    fn kepler_position_ex_earth_at_j2000() {
        let mut out = empty_cart();
        let s = siderust_kepler_position_ex(earth_orbit(), 1, J2000, &mut out);
        assert_eq!(s, SiderustStatus::Ok);
        // Earth should be ~1 AU from the Sun
        let r = (out.x * out.x + out.y * out.y + out.z * out.z).sqrt();
        assert!(r > 0.98 && r < 1.02, "r = {r}");
        assert_eq!(out.frame, SiderustFrame::EclipticMeanJ2000);
        assert_eq!(out.center, SiderustCenter::Heliocentric);
    }

    #[test]
    fn kepler_position_ex_null_out() {
        let s = siderust_kepler_position_ex(earth_orbit(), 1, J2000, ptr::null_mut());
        assert_eq!(s, SiderustStatus::NullPointer);
    }

    #[test]
    fn mean_motion_position_earth_at_j2000() {
        let mut out = empty_cart();
        let s = siderust_mean_motion_position(earth_mean_motion_orbit(), J2000, &mut out);
        assert_eq!(s, SiderustStatus::Ok);
        let r = (out.x * out.x + out.y * out.y + out.z * out.z).sqrt();
        assert!(r > 0.98 && r < 1.02, "r = {r}");
        assert_eq!(out.frame, SiderustFrame::EclipticMeanJ2000);
    }

    #[test]
    fn conic_position_rejects_parabolic() {
        let orbit = SiderustConicOrbit {
            periapsis_distance_au: 1.0,
            eccentricity: 1.0,
            inclination_deg: 0.0,
            lon_ascending_node_deg: 0.0,
            arg_periapsis_deg: 0.0,
            mean_anomaly_deg: 0.0,
            epoch_jd: J2000,
        };
        let mut out = empty_cart();
        let s = siderust_conic_position(orbit, J2000, &mut out);
        assert_eq!(s, SiderustStatus::InvalidArgument);
    }

    #[test]
    fn conic_position_rejects_invalid_periapsis() {
        // Zero periapsis distance must be caught at the FFI boundary, not silently
        // produce NaN coordinates.
        let orbit = SiderustConicOrbit {
            periapsis_distance_au: 0.0,
            eccentricity: 0.5,
            inclination_deg: 0.0,
            lon_ascending_node_deg: 0.0,
            arg_periapsis_deg: 0.0,
            mean_anomaly_deg: 0.0,
            epoch_jd: J2000,
        };
        let mut out = empty_cart();
        let s = siderust_conic_position(orbit, J2000, &mut out);
        assert_eq!(s, SiderustStatus::InvalidArgument);
    }

    #[test]
    fn mean_motion_position_rejects_invalid_semi_major_axis() {
        let orbit = SiderustMeanMotionOrbit {
            semi_major_axis_au: -1.0,
            eccentricity: 0.0167,
            inclination_deg: 0.0,
            lon_ascending_node_deg: 0.0,
            arg_periapsis_deg: 0.0,
            mean_motion_deg_per_day: 0.9856,
            epoch_jd: J2000,
        };
        let mut out = empty_cart();
        let s = siderust_mean_motion_position(orbit, J2000, &mut out);
        assert_eq!(s, SiderustStatus::InvalidArgument);
    }

    #[test]
    fn prepared_orbit_create_rejects_invalid_semi_major_axis() {
        // Negative SMA must return InvalidArgument, not silently create a broken handle.
        let orbit = SiderustOrbit {
            semi_major_axis_au: -1.0,
            eccentricity: 0.0167,
            inclination_deg: 0.0,
            lon_ascending_node_deg: 0.0,
            arg_periapsis_deg: 0.0,
            mean_anomaly_deg: 0.0,
            epoch_jd: J2000,
        };
        let mut handle = std::ptr::null_mut();
        let s = siderust_prepared_orbit_create(orbit, &mut handle);
        assert_eq!(s, SiderustStatus::InvalidArgument);
        assert!(handle.is_null());
    }

    #[test]
    fn prepared_orbit_create_rejects_hyperbolic() {
        let orbit = SiderustOrbit {
            semi_major_axis_au: 1.0,
            eccentricity: 1.5, // hyperbolic — KeplerianOrbit requires e < 1
            inclination_deg: 0.0,
            lon_ascending_node_deg: 0.0,
            arg_periapsis_deg: 0.0,
            mean_anomaly_deg: 0.0,
            epoch_jd: J2000,
        };
        let mut handle = std::ptr::null_mut();
        let s = siderust_prepared_orbit_create(orbit, &mut handle);
        assert_eq!(s, SiderustStatus::InvalidArgument);
        assert!(handle.is_null());
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
            length_unit: SiderustLengthUnit::AU,
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
            length_unit: SiderustLengthUnit::AU,
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
            length_unit: SiderustLengthUnit::AU,
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
            length_unit: SiderustLengthUnit::AU,
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
            length_unit: SiderustLengthUnit::AU,
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
            length_unit: SiderustLengthUnit::AU,
        };
        let mut out = empty_cart();
        let s = siderust_to_bodycentric(pos, earth_bodycentric_params(), J2000, &mut out);
        assert_eq!(s, SiderustStatus::Ok);
        assert_eq!(out.center, SiderustCenter::Bodycentric);
    }
}
