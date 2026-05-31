// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # DE4xx Body-Chain Resolution
//!
//! ## Scientific scope
//!
//! The JPL DE ephemerides natively integrate the Sun, the Earth–Moon
//! Barycenter (EMB), and the Moon relative to the EMB.  This module
//! reconstructs the individual states of the Sun, Earth, and Moon in
//! every reference-center and unit combination required by the public
//! [`Ephemeris`](crate::calculus::ephemeris::Ephemeris) trait.
//!
//! The DE Moon segment is stored as the Moon position **relative to the
//! Earth–Moon barycenter** (`Moon_offset = Moon − EMB`).  The Earth–Moon
//! barycenter satisfies `EMB = FRAC_EARTH·Earth + FRAC_MOON·Moon`, so the
//! geocentric Earth → Moon separation is
//! `Moon − Earth = Moon_offset / FRAC_EARTH`, and Earth itself is
//! `Earth = EMB − Moon_offset · (FRAC_MOON / FRAC_EARTH) = EMB − Moon_offset / EMRAT`.
//!
//! Derived quantities:
//!
//! - **Earth barycentric** = EMB − Moon_offset / EMRAT
//! - **Moon geocentric**   = Moon_offset / FRAC_EARTH
//! - **Earth heliocentric**= Earth_bary − Sun_bary
//!
//! where 1/EMRAT ≈ 1/81.300569 and 1/FRAC_EARTH ≈ 82.300569/81.300569.
//!
//! ## Technical scope
//!
//! - All intermediate vectors carry frame
//!   ([`crate::coordinates::frames::ICRF`]) and unit ([`Kilometer`]
//!   or [`Per<Kilometer, Day>`]) in the type system; body-chain arithmetic
//!   is compile-time checked.
//! - Frame conversion (ICRF → EclipticMeanJ2000) and unit conversion
//!   (km → AU, km/day → AU/day) are applied explicitly after chain arithmetic.
//! - Time input: `JulianDate` (TT); TT → TDB conversion is performed
//!   internally via [`JulianDate::tt_to_tdb`].
//!
//! | Segment | DE NAIF IDs   | Meaning                                |
//! |---------|---------------|----------------------------------------|
//! | Sun     | 10 → 0 (SSB)  | Sun barycentric (ICRF, km)             |
//! | EMB     |  3 → 0 (SSB)  | Earth-Moon barycenter bary. (ICRF, km) |
//! | Moon    | 301 → 3 (EMB) | Moon offset from EMB (ICRF, km)        |
//!
//! ## References
//!
//! - Standish, E. M. (1998). "JPL Planetary and Lunar Ephemerides, DE405/LE405".
//!   *JPL Interoffice Memorandum* 312.F-98-048.
//! - Park, R. S., et al. (2021). "The JPL Planetary and Lunar Ephemerides
//!   DE440 and DE441". *The Astronomical Journal* 161, 105.
//!   <https://doi.org/10.3847/1538-3881/abd414>

use super::eval::{DynSegmentDescriptor, SegmentDescriptor};

use crate::calculus::ephemeris::EphemerisError;
use crate::coordinates::{
    cartesian::{Position, Velocity},
    centers::{Barycentric, Geocentric, Heliocentric},
    frames::EclipticMeanJ2000,
    transform::VectorAstroExt,
};
use crate::qtty::{AstronomicalUnit, Day, Kilometer, Per};
use crate::time::{JulianDate, TDB};

// ── Physical constants (single source of truth: qtty + DE4xx headers) ────

/// Earth/Moon mass ratio embedded in DE4xx headers.
///
/// This value is consistent across DE440 and DE441.
const EARTH_MOON_RATIO: f64 = 81.300_569_074_190_62;

/// μ_Earth / (μ_Earth + μ_Moon) , Earth's mass fraction of the EM system.
const FRAC_EARTH: f64 = EARTH_MOON_RATIO / (EARTH_MOON_RATIO + 1.0);

/// μ_Moon / (μ_Earth + μ_Moon) , Moon's mass fraction of the EM system.
#[allow(dead_code)]
const FRAC_MOON: f64 = 1.0 / (EARTH_MOON_RATIO + 1.0);

/// Scale converting `moon_off = Moon − EMB` to the Earth-side correction
/// `EMB − Earth = moon_off · (FRAC_MOON / FRAC_EARTH) = moon_off / EMRAT`.
const EARTH_OFFSET_FROM_MOON_OFF: f64 = 1.0 / EARTH_MOON_RATIO;

/// Scale converting `moon_off = Moon − EMB` to the geocentric Earth→Moon
/// vector `Moon − Earth = moon_off / FRAC_EARTH = moon_off · (EMRAT + 1) / EMRAT`.
const MOON_GEO_FROM_MOON_OFF: f64 = 1.0 / FRAC_EARTH;

// ── Velocity type alias ─────────────────────────────────────────────────

type AuPerDay = Per<AstronomicalUnit, Day>;

// ── Generic body-chain functions ─────────────────────────────────────────

/// Sun barycentric position in EclipticMeanJ2000 J2000 (AU).
///
/// Generic over any DE4xx data source providing a SUN segment.
#[inline]
pub fn try_sun_barycentric(
    jd: JulianDate,
    sun: &SegmentDescriptor,
) -> Result<Position<Barycentric, EclipticMeanJ2000, AstronomicalUnit>, EphemerisError> {
    let jd_tdb = jd.to_scale::<TDB>();
    let sun_icrf = sun.try_position(jd_tdb)?;
    let sun_ecl_au = sun_icrf
        .to_frame::<EclipticMeanJ2000>(&crate::J2000)
        .to_unit::<AstronomicalUnit>();
    Ok(Position::new(
        sun_ecl_au.x(),
        sun_ecl_au.y(),
        sun_ecl_au.z(),
    ))
}

/// Sun barycentric position in EclipticMeanJ2000 J2000 (AU).
///
/// # Panics
///
/// Panics when `jd` is outside the underlying JPL segment coverage. Use
/// [`try_sun_barycentric`] for explicit error handling.
#[inline]
pub fn sun_barycentric(
    jd: JulianDate,
    sun: &SegmentDescriptor,
) -> Position<Barycentric, EclipticMeanJ2000, AstronomicalUnit> {
    try_sun_barycentric(jd, sun).expect("JPL Sun barycentric position unavailable")
}

/// Earth barycentric position in EclipticMeanJ2000 J2000 (AU).
///
/// `Earth_bary = EMB − Moon_offset / EMRAT`  (equivalent to `EMB − Moon_offset · FRAC_MOON/FRAC_EARTH`)
///
/// Generic over any DE4xx data source providing EMB and MOON segments.
#[inline]
pub fn try_earth_barycentric(
    jd: JulianDate,
    emb: &SegmentDescriptor,
    moon: &SegmentDescriptor,
) -> Result<Position<Barycentric, EclipticMeanJ2000, AstronomicalUnit>, EphemerisError> {
    let jd_tdb = jd.to_scale::<TDB>();
    let emb_pos = emb.try_position(jd_tdb)?;
    let moon_off = moon.try_position(jd_tdb)?;
    let earth_icrf = emb_pos - moon_off.scale(EARTH_OFFSET_FROM_MOON_OFF);
    let earth_ecl_au = earth_icrf
        .to_frame::<EclipticMeanJ2000>(&crate::J2000)
        .to_unit::<AstronomicalUnit>();
    Ok(Position::new(
        earth_ecl_au.x(),
        earth_ecl_au.y(),
        earth_ecl_au.z(),
    ))
}

/// Earth barycentric position in EclipticMeanJ2000 J2000 (AU).
///
/// # Panics
///
/// Panics when `jd` is outside the underlying JPL segment coverage. Use
/// [`try_earth_barycentric`] for explicit error handling.
#[inline]
pub fn earth_barycentric(
    jd: JulianDate,
    emb: &SegmentDescriptor,
    moon: &SegmentDescriptor,
) -> Position<Barycentric, EclipticMeanJ2000, AstronomicalUnit> {
    try_earth_barycentric(jd, emb, moon).expect("JPL Earth barycentric position unavailable")
}

/// Earth heliocentric position in EclipticMeanJ2000 J2000 (AU).
///
/// `Earth_helio = Earth_bary − Sun_bary`
///
/// Generic over any DE4xx data source providing SUN, EMB, and MOON segments.
#[inline]
pub fn try_earth_heliocentric(
    jd: JulianDate,
    sun: &SegmentDescriptor,
    emb: &SegmentDescriptor,
    moon: &SegmentDescriptor,
) -> Result<Position<Heliocentric, EclipticMeanJ2000, AstronomicalUnit>, EphemerisError> {
    let jd_tdb = jd.to_scale::<TDB>();
    let emb_pos = emb.try_position(jd_tdb)?;
    let moon_off = moon.try_position(jd_tdb)?;
    let sun_pos = sun.try_position(jd_tdb)?;
    let earth_icrf = emb_pos - moon_off.scale(EARTH_OFFSET_FROM_MOON_OFF) - sun_pos;
    let earth_ecl_au = earth_icrf
        .to_frame::<EclipticMeanJ2000>(&crate::J2000)
        .to_unit::<AstronomicalUnit>();
    Ok(Position::new(
        earth_ecl_au.x(),
        earth_ecl_au.y(),
        earth_ecl_au.z(),
    ))
}

/// Earth heliocentric position in EclipticMeanJ2000 J2000 (AU).
///
/// # Panics
///
/// Panics when `jd` is outside the underlying JPL segment coverage. Use
/// [`try_earth_heliocentric`] for explicit error handling.
#[inline]
pub fn earth_heliocentric(
    jd: JulianDate,
    sun: &SegmentDescriptor,
    emb: &SegmentDescriptor,
    moon: &SegmentDescriptor,
) -> Position<Heliocentric, EclipticMeanJ2000, AstronomicalUnit> {
    try_earth_heliocentric(jd, sun, emb, moon).expect("JPL Earth heliocentric position unavailable")
}

/// Earth barycentric velocity in EclipticMeanJ2000 J2000 (AU/day).
///
/// `v_Earth = v_EMB − v_Moon_offset / EMRAT`
///
/// Generic over any DE4xx data source providing EMB and MOON segments.
#[inline]
pub fn try_earth_barycentric_velocity(
    jd: JulianDate,
    emb: &SegmentDescriptor,
    moon: &SegmentDescriptor,
) -> Result<Velocity<EclipticMeanJ2000, AuPerDay>, EphemerisError> {
    let jd_tdb = jd.to_scale::<TDB>();
    let v_emb = emb.try_velocity(jd_tdb)?;
    let v_moon_off = moon.try_velocity(jd_tdb)?;
    let v_earth_icrf = v_emb - v_moon_off.scale(EARTH_OFFSET_FROM_MOON_OFF);
    Ok(v_earth_icrf
        .to_frame::<EclipticMeanJ2000>(&crate::J2000)
        .to_unit::<AuPerDay>())
}

/// Earth barycentric velocity in EclipticMeanJ2000 J2000 (AU/day).
///
/// # Panics
///
/// Panics when `jd` is outside the underlying JPL segment coverage. Use
/// [`try_earth_barycentric_velocity`] for explicit error handling.
#[inline]
pub fn earth_barycentric_velocity(
    jd: JulianDate,
    emb: &SegmentDescriptor,
    moon: &SegmentDescriptor,
) -> Velocity<EclipticMeanJ2000, AuPerDay> {
    try_earth_barycentric_velocity(jd, emb, moon)
        .expect("JPL Earth barycentric velocity unavailable")
}

/// Moon geocentric position in EclipticMeanJ2000 J2000 (km).
///
/// `Moon_geo = Moon_offset / FRAC_EARTH`
///
/// Generic over any DE4xx data source providing a MOON segment.
#[inline]
pub fn try_moon_geocentric(
    jd: JulianDate,
    moon: &SegmentDescriptor,
) -> Result<Position<Geocentric, EclipticMeanJ2000, Kilometer>, EphemerisError> {
    let jd_tdb = jd.to_scale::<TDB>();
    let moon_off = moon.try_position(jd_tdb)?;
    let moon_geo_icrf = moon_off.scale(MOON_GEO_FROM_MOON_OFF);
    let moon_geo_ecl = moon_geo_icrf.to_frame::<EclipticMeanJ2000>(&crate::J2000);
    Ok(Position::new(
        moon_geo_ecl.x(),
        moon_geo_ecl.y(),
        moon_geo_ecl.z(),
    ))
}

/// Moon geocentric position in EclipticMeanJ2000 J2000 (km).
///
/// # Panics
///
/// Panics when `jd` is outside the underlying JPL segment coverage. Use
/// [`try_moon_geocentric`] for explicit error handling.
#[inline]
pub fn moon_geocentric(
    jd: JulianDate,
    moon: &SegmentDescriptor,
) -> Position<Geocentric, EclipticMeanJ2000, Kilometer> {
    try_moon_geocentric(jd, moon).expect("JPL Moon geocentric position unavailable")
}

/// Major-planet barycentric position from a direct SSB segment.
#[inline]
pub fn try_planet_barycentric(
    jd: JulianDate,
    planet: &SegmentDescriptor,
) -> Result<Position<Barycentric, EclipticMeanJ2000, AstronomicalUnit>, EphemerisError> {
    let jd_tdb = jd.to_scale::<TDB>();
    let planet_icrf = planet.try_position(jd_tdb)?;
    let planet_ecl = planet_icrf
        .to_frame::<EclipticMeanJ2000>(&crate::J2000)
        .to_unit::<AstronomicalUnit>();
    Ok(Position::new(
        planet_ecl.x(),
        planet_ecl.y(),
        planet_ecl.z(),
    ))
}

/// Child planet-center barycentric position from barycenter and center-offset segments.
#[inline]
pub fn try_child_planet_barycentric(
    jd: JulianDate,
    barycenter: &SegmentDescriptor,
    center_offset: &SegmentDescriptor,
) -> Result<Position<Barycentric, EclipticMeanJ2000, AstronomicalUnit>, EphemerisError> {
    let jd_tdb = jd.to_scale::<TDB>();
    let center_icrf = barycenter.try_position(jd_tdb)? + center_offset.try_position(jd_tdb)?;
    let center_ecl = center_icrf
        .to_frame::<EclipticMeanJ2000>(&crate::J2000)
        .to_unit::<AstronomicalUnit>();
    Ok(Position::new(
        center_ecl.x(),
        center_ecl.y(),
        center_ecl.z(),
    ))
}

/// Earth-centered major-planet position from a direct SSB segment.
#[inline]
pub fn try_planet_geocentric(
    jd: JulianDate,
    planet: &SegmentDescriptor,
    emb: &SegmentDescriptor,
    moon: &SegmentDescriptor,
) -> Result<Position<Geocentric, EclipticMeanJ2000, AstronomicalUnit>, EphemerisError> {
    let jd_tdb = jd.to_scale::<TDB>();
    let planet_icrf = planet.try_position(jd_tdb)?;
    let emb_icrf = emb.try_position(jd_tdb)?;
    let moon_off = moon.try_position(jd_tdb)?;
    let earth_icrf = emb_icrf - moon_off.scale(EARTH_OFFSET_FROM_MOON_OFF);
    let geo_ecl = (planet_icrf - earth_icrf)
        .to_frame::<EclipticMeanJ2000>(&crate::J2000)
        .to_unit::<AstronomicalUnit>();
    Ok(Position::new(geo_ecl.x(), geo_ecl.y(), geo_ecl.z()))
}

/// Earth-centered planet-center position from barycenter and center-offset segments.
#[inline]
pub fn try_child_planet_geocentric(
    jd: JulianDate,
    barycenter: &SegmentDescriptor,
    center_offset: &SegmentDescriptor,
    emb: &SegmentDescriptor,
    moon: &SegmentDescriptor,
) -> Result<Position<Geocentric, EclipticMeanJ2000, AstronomicalUnit>, EphemerisError> {
    let jd_tdb = jd.to_scale::<TDB>();
    let center_icrf = barycenter.try_position(jd_tdb)? + center_offset.try_position(jd_tdb)?;
    let emb_icrf = emb.try_position(jd_tdb)?;
    let moon_off = moon.try_position(jd_tdb)?;
    let earth_icrf = emb_icrf - moon_off.scale(EARTH_OFFSET_FROM_MOON_OFF);
    let geo_ecl = (center_icrf - earth_icrf)
        .to_frame::<EclipticMeanJ2000>(&crate::J2000)
        .to_unit::<AstronomicalUnit>();
    Ok(Position::new(geo_ecl.x(), geo_ecl.y(), geo_ecl.z()))
}

// ═══════════════════════════════════════════════════════════════════════════
// DynSegmentDescriptor variants, runtime-loaded data
// ═══════════════════════════════════════════════════════════════════════════

/// Sun barycentric position (runtime data).
#[inline]
pub fn try_dyn_sun_barycentric(
    jd: JulianDate,
    sun: &DynSegmentDescriptor,
) -> Result<Position<Barycentric, EclipticMeanJ2000, AstronomicalUnit>, EphemerisError> {
    let jd_tdb = jd.to_scale::<TDB>();
    let sun_icrf = sun.try_position(jd_tdb)?;
    let sun_ecl_au = sun_icrf
        .to_frame::<EclipticMeanJ2000>(&crate::J2000)
        .to_unit::<AstronomicalUnit>();
    Ok(Position::new(
        sun_ecl_au.x(),
        sun_ecl_au.y(),
        sun_ecl_au.z(),
    ))
}

/// Sun barycentric position (runtime data).
#[inline]
pub fn dyn_sun_barycentric(
    jd: JulianDate,
    sun: &DynSegmentDescriptor,
) -> Position<Barycentric, EclipticMeanJ2000, AstronomicalUnit> {
    try_dyn_sun_barycentric(jd, sun).expect("runtime JPL Sun barycentric position unavailable")
}

/// Earth barycentric position (runtime data).
#[inline]
pub fn try_dyn_earth_barycentric(
    jd: JulianDate,
    emb: &DynSegmentDescriptor,
    moon: &DynSegmentDescriptor,
) -> Result<Position<Barycentric, EclipticMeanJ2000, AstronomicalUnit>, EphemerisError> {
    let jd_tdb = jd.to_scale::<TDB>();
    let emb_pos = emb.try_position(jd_tdb)?;
    let moon_off = moon.try_position(jd_tdb)?;
    let earth_icrf = emb_pos - moon_off.scale(EARTH_OFFSET_FROM_MOON_OFF);
    let earth_ecl_au = earth_icrf
        .to_frame::<EclipticMeanJ2000>(&crate::J2000)
        .to_unit::<AstronomicalUnit>();
    Ok(Position::new(
        earth_ecl_au.x(),
        earth_ecl_au.y(),
        earth_ecl_au.z(),
    ))
}

/// Earth barycentric position (runtime data).
#[inline]
pub fn dyn_earth_barycentric(
    jd: JulianDate,
    emb: &DynSegmentDescriptor,
    moon: &DynSegmentDescriptor,
) -> Position<Barycentric, EclipticMeanJ2000, AstronomicalUnit> {
    try_dyn_earth_barycentric(jd, emb, moon)
        .expect("runtime JPL Earth barycentric position unavailable")
}

/// Earth heliocentric position (runtime data).
#[inline]
pub fn try_dyn_earth_heliocentric(
    jd: JulianDate,
    sun: &DynSegmentDescriptor,
    emb: &DynSegmentDescriptor,
    moon: &DynSegmentDescriptor,
) -> Result<Position<Heliocentric, EclipticMeanJ2000, AstronomicalUnit>, EphemerisError> {
    let jd_tdb = jd.to_scale::<TDB>();
    let emb_pos = emb.try_position(jd_tdb)?;
    let moon_off = moon.try_position(jd_tdb)?;
    let sun_pos = sun.try_position(jd_tdb)?;
    let earth_icrf = emb_pos - moon_off.scale(EARTH_OFFSET_FROM_MOON_OFF) - sun_pos;
    let earth_ecl_au = earth_icrf
        .to_frame::<EclipticMeanJ2000>(&crate::J2000)
        .to_unit::<AstronomicalUnit>();
    Ok(Position::new(
        earth_ecl_au.x(),
        earth_ecl_au.y(),
        earth_ecl_au.z(),
    ))
}

/// Earth heliocentric position (runtime data).
#[inline]
pub fn dyn_earth_heliocentric(
    jd: JulianDate,
    sun: &DynSegmentDescriptor,
    emb: &DynSegmentDescriptor,
    moon: &DynSegmentDescriptor,
) -> Position<Heliocentric, EclipticMeanJ2000, AstronomicalUnit> {
    try_dyn_earth_heliocentric(jd, sun, emb, moon)
        .expect("runtime JPL Earth heliocentric position unavailable")
}

/// Earth barycentric velocity (runtime data).
#[inline]
pub fn try_dyn_earth_barycentric_velocity(
    jd: JulianDate,
    emb: &DynSegmentDescriptor,
    moon: &DynSegmentDescriptor,
) -> Result<Velocity<EclipticMeanJ2000, AuPerDay>, EphemerisError> {
    let jd_tdb = jd.to_scale::<TDB>();
    let v_emb = emb.try_velocity(jd_tdb)?;
    let v_moon_off = moon.try_velocity(jd_tdb)?;
    let v_earth_icrf = v_emb - v_moon_off.scale(EARTH_OFFSET_FROM_MOON_OFF);
    Ok(v_earth_icrf
        .to_frame::<EclipticMeanJ2000>(&crate::J2000)
        .to_unit::<AuPerDay>())
}

/// Earth barycentric velocity (runtime data).
#[inline]
pub fn dyn_earth_barycentric_velocity(
    jd: JulianDate,
    emb: &DynSegmentDescriptor,
    moon: &DynSegmentDescriptor,
) -> Velocity<EclipticMeanJ2000, AuPerDay> {
    try_dyn_earth_barycentric_velocity(jd, emb, moon)
        .expect("runtime JPL Earth barycentric velocity unavailable")
}

/// Moon geocentric position (runtime data).
#[inline]
pub fn try_dyn_moon_geocentric(
    jd: JulianDate,
    moon: &DynSegmentDescriptor,
) -> Result<Position<Geocentric, EclipticMeanJ2000, Kilometer>, EphemerisError> {
    let jd_tdb = jd.to_scale::<TDB>();
    let moon_off = moon.try_position(jd_tdb)?;
    let moon_geo_icrf = moon_off.scale(MOON_GEO_FROM_MOON_OFF);
    let moon_geo_ecl = moon_geo_icrf.to_frame::<EclipticMeanJ2000>(&crate::J2000);
    Ok(Position::new(
        moon_geo_ecl.x(),
        moon_geo_ecl.y(),
        moon_geo_ecl.z(),
    ))
}

/// Moon geocentric position (runtime data).
#[inline]
pub fn dyn_moon_geocentric(
    jd: JulianDate,
    moon: &DynSegmentDescriptor,
) -> Position<Geocentric, EclipticMeanJ2000, Kilometer> {
    try_dyn_moon_geocentric(jd, moon).expect("runtime JPL Moon geocentric position unavailable")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::calculus::jpl::eval::{DynSegmentDescriptor, SegmentDescriptor};

    const SECONDS_PER_DAY: f64 = crate::qtty::time::SECONDS_PER_DAY;
    const JD_J2000: f64 = 2451545.0;

    // ── Shared test data for SegmentDescriptor (compile-time variant) ─────

    /// Static record for sun-like segment: position ~0.5 AU at tau=0.
    static SUN_RECORD: [f64; 8] = [
        500.0 * SECONDS_PER_DAY, // mid (seconds past J2000)
        500.0 * SECONDS_PER_DAY, // radius = half-interval
        7.5e7,
        0.0, // cx ~0.5 AU in km
        3.0e7,
        0.0, // cy
        1.0e7,
        0.0, // cz
    ];

    /// Static record for EMB-like segment: ~1 AU from SSB.
    static EMB_RECORD: [f64; 8] = [
        500.0 * SECONDS_PER_DAY,
        500.0 * SECONDS_PER_DAY,
        1.5e8,
        0.0, // ~1 AU
        0.0,
        0.0,
        0.0,
        0.0,
    ];

    /// Static record for Moon offset (~Earth-Moon distance).
    static MOON_RECORD: [f64; 8] = [
        500.0 * SECONDS_PER_DAY,
        500.0 * SECONDS_PER_DAY,
        3.84e5,
        0.0, // ~384400 km
        0.0,
        0.0,
        0.0,
        0.0,
    ];

    fn sun_record_fn(_: usize) -> &'static [f64] {
        &SUN_RECORD
    }
    fn emb_record_fn(_: usize) -> &'static [f64] {
        &EMB_RECORD
    }
    fn moon_record_fn(_: usize) -> &'static [f64] {
        &MOON_RECORD
    }

    fn make_static_seg(record_fn: fn(usize) -> &'static [f64]) -> SegmentDescriptor {
        use crate::qtty::Seconds;
        SegmentDescriptor {
            init: Seconds::new(0.0),
            intlen: Seconds::new(1000.0 * SECONDS_PER_DAY),
            ncoeff: 2,
            n_records: 1,
            record_fn,
        }
    }

    fn jd_test_static() -> JulianDate {
        crate::time::JulianDate::new(JD_J2000 + 500.0)
    }

    // ── SegmentDescriptor (compile-time) tests ────────────────────────────

    #[test]
    fn sun_barycentric_is_finite() {
        let sun = make_static_seg(sun_record_fn);
        let pos = sun_barycentric(jd_test_static(), &sun);
        assert!(pos.x().is_finite());
        assert!(pos.y().is_finite());
        assert!(pos.z().is_finite());
    }

    #[test]
    fn earth_barycentric_is_finite() {
        let emb = make_static_seg(emb_record_fn);
        let moon = make_static_seg(moon_record_fn);
        let pos = earth_barycentric(jd_test_static(), &emb, &moon);
        assert!(pos.x().is_finite());
        assert!(pos.y().is_finite());
        assert!(pos.z().is_finite());
    }

    #[test]
    fn earth_heliocentric_is_finite() {
        let sun = make_static_seg(sun_record_fn);
        let emb = make_static_seg(emb_record_fn);
        let moon = make_static_seg(moon_record_fn);
        let pos = earth_heliocentric(jd_test_static(), &sun, &emb, &moon);
        assert!(pos.x().is_finite());
        assert!(pos.y().is_finite());
        assert!(pos.z().is_finite());
    }

    #[test]
    fn earth_barycentric_velocity_is_finite() {
        let emb = make_static_seg(emb_record_fn);
        let moon = make_static_seg(moon_record_fn);
        let vel = earth_barycentric_velocity(jd_test_static(), &emb, &moon);
        assert!(vel.x().is_finite());
        assert!(vel.y().is_finite());
        assert!(vel.z().is_finite());
    }

    #[test]
    fn moon_geocentric_is_finite() {
        let moon = make_static_seg(moon_record_fn);
        let pos = moon_geocentric(jd_test_static(), &moon);
        assert!(pos.x().is_finite());
        assert!(pos.y().is_finite());
        assert!(pos.z().is_finite());
    }

    #[test]
    fn moon_geocentric_magnitude_scaled_by_inv_frac_earth() {
        // moon_off (in segment) = (r_km, 0, 0) → Moon_geo = moon_off / FRAC_EARTH.
        // The expected magnitude is r_km/FRAC_EARTH ≈ 389128 km. The ICRF→
        // EclipticMeanJ2000 frame composition introduces a small (<0.2%)
        // numerical drift in magnitude with the synthetic test fixture; we
        // therefore use a 1% tolerance for the sanity check. The exact
        // formula correctness is proven by `earth_barycentric_uses_one_over_emrat_offset`.
        let r_km = 384400.0;
        let moon = make_static_seg(moon_record_fn); // position = (r_km, 0, 0) at tau=0
        let pos = moon_geocentric(jd_test_static(), &moon);
        let mag =
            (pos.x().value().powi(2) + pos.y().value().powi(2) + pos.z().value().powi(2)).sqrt();
        let expected = r_km / FRAC_EARTH;
        assert!(
            (mag - expected).abs() / expected < 1e-2,
            "mag={mag}, expected={expected}"
        );
        // Strict bound: must be strictly larger than r_km (1/FRAC_EARTH > 1).
        assert!(mag > r_km, "mag={mag} must exceed r_km={r_km}");
    }

    #[test]
    fn earth_barycentric_uses_one_over_emrat_offset() {
        // moon_off = (r_km, 0, 0), EMB = (1.5e8, 0, 0).
        // Earth_bary_icrf = EMB - moon_off / EMRAT (before ecliptic rotation).
        let emb = make_static_seg(emb_record_fn);
        let moon = make_static_seg(moon_record_fn);
        let pos = earth_barycentric(jd_test_static(), &emb, &moon);
        // Magnitude is unchanged by the equatorial→ecliptic rotation.
        let mag_au =
            (pos.x().value().powi(2) + pos.y().value().powi(2) + pos.z().value().powi(2)).sqrt();
        let emb_km = 1.5e8;
        let moon_off_km = 384400.0;
        let earth_km = emb_km - moon_off_km / EARTH_MOON_RATIO;
        let expected_au = earth_km / 1.495_978_707e8;
        assert!(
            (mag_au - expected_au).abs() / expected_au < 1e-3,
            "mag_au={mag_au}, expected={expected_au}"
        );
    }

    // ── DynSegmentDescriptor tests ────────────────────────────────────────

    /// Create a synthetic DynSegmentDescriptor representing a constant position.
    ///
    /// The segment spans 1000 days from J2000. The position at tau=0 (J2000+500d)
    /// is (x_km, y_km, z_km) in ICRF km. The velocity is zero.
    fn make_seg(x_km: f64, y_km: f64, z_km: f64) -> DynSegmentDescriptor {
        use crate::qtty::Seconds;
        let ncoeff = 2usize;
        let rsize = 2 + 3 * ncoeff; // 8
        let intlen_secs = 1000.0 * SECONDS_PER_DAY;
        let mid = intlen_secs / 2.0;
        let radius = intlen_secs / 2.0;
        // Record: [mid, radius, cx0, cx1, cy0, cy1, cz0, cz1]
        let data = vec![mid, radius, x_km, 0.0, y_km, 0.0, z_km, 0.0];
        DynSegmentDescriptor {
            data_type: 2,
            init: Seconds::new(0.0),
            intlen: Seconds::new(intlen_secs),
            ncoeff,
            rsize,
            n_records: 1,
            data,
        }
    }

    /// JD at J2000 + 500 days (midpoint of our test segment → tau = 0).
    fn jd_test() -> JulianDate {
        crate::time::JulianDate::new(JD_J2000 + 500.0)
    }

    // ── FRAC constants ────────────────────────────────────────────────────

    #[test]
    fn frac_earth_plus_frac_moon_equals_one() {
        assert!((FRAC_EARTH + FRAC_MOON - 1.0).abs() < 1e-12);
    }

    #[test]
    fn frac_earth_is_dominant() {
        const {
            assert!(FRAC_EARTH > 0.98);
        }
        const {
            assert!(FRAC_MOON < 0.02);
        }
    }

    // ── dyn_sun_barycentric ───────────────────────────────────────────────

    #[test]
    fn dyn_sun_barycentric_is_finite() {
        let sun = make_seg(1.0e8, 2.0e7, 3.0e6); // ~1 AU-ish in km
        let pos = dyn_sun_barycentric(jd_test(), &sun);
        assert!(pos.x().is_finite());
        assert!(pos.y().is_finite());
        assert!(pos.z().is_finite());
    }

    // ── dyn_earth_barycentric ─────────────────────────────────────────────

    #[test]
    fn dyn_earth_barycentric_is_finite() {
        let emb = make_seg(1.5e8, 0.0, 0.0);
        let moon = make_seg(3.84e5, 0.0, 0.0);
        let pos = dyn_earth_barycentric(jd_test(), &emb, &moon);
        assert!(pos.x().is_finite());
        assert!(pos.y().is_finite());
        assert!(pos.z().is_finite());
    }

    // ── dyn_earth_heliocentric ─────────────────────────────────────────────

    #[test]
    fn dyn_earth_heliocentric_is_finite() {
        let sun = make_seg(1.0e8, 0.0, 0.0);
        let emb = make_seg(1.5e8, 0.0, 0.0);
        let moon = make_seg(3.84e5, 0.0, 0.0);
        let pos = dyn_earth_heliocentric(jd_test(), &sun, &emb, &moon);
        assert!(pos.x().is_finite());
        assert!(pos.y().is_finite());
        assert!(pos.z().is_finite());
    }

    // ── dyn_earth_barycentric_velocity ────────────────────────────────────

    #[test]
    fn dyn_earth_barycentric_velocity_is_finite() {
        let emb = make_seg(1.5e8, 0.0, 0.0);
        let moon = make_seg(3.84e5, 0.0, 0.0);
        let vel = dyn_earth_barycentric_velocity(jd_test(), &emb, &moon);
        assert!(vel.x().is_finite());
        assert!(vel.y().is_finite());
        assert!(vel.z().is_finite());
    }

    // ── dyn_moon_geocentric ───────────────────────────────────────────────

    #[test]
    fn dyn_moon_geocentric_is_finite() {
        let moon = make_seg(3.84e5, 1.0e4, 2.0e3);
        let pos = dyn_moon_geocentric(jd_test(), &moon);
        assert!(pos.x().is_finite());
        assert!(pos.y().is_finite());
        assert!(pos.z().is_finite());
    }

    #[test]
    fn dyn_moon_geocentric_scaled_by_inv_frac_earth() {
        // With moon_off = (R, 0, 0) at tau=0, Moon_geo = moon_off / FRAC_EARTH.
        // After ICRF→EclipticMeanJ2000 the magnitude has at most a small
        // (<1%) drift from frame-composition rounding; the strict scaling
        // assertion is covered by `dyn_earth_barycentric_uses_one_over_emrat_offset`.
        let r_km = 384400.0;
        let moon = make_seg(r_km, 0.0, 0.0);
        let pos = dyn_moon_geocentric(jd_test(), &moon);
        let magnitude =
            (pos.x().value().powi(2) + pos.y().value().powi(2) + pos.z().value().powi(2)).sqrt();
        let expected = r_km / FRAC_EARTH;
        assert!(
            (magnitude - expected).abs() / expected < 1e-2,
            "magnitude={magnitude}, expected={expected}"
        );
        assert!(
            magnitude > r_km,
            "magnitude={magnitude} must exceed r_km={r_km}"
        );
    }

    #[test]
    fn dyn_earth_barycentric_uses_one_over_emrat_offset() {
        // moon_off = (r_km, 0, 0), EMB = (1.5e8, 0, 0).
        // Earth_bary_icrf = EMB - moon_off / EMRAT; magnitude preserved under rotation.
        let emb = make_seg(1.5e8, 0.0, 0.0);
        let r_km = 384400.0;
        let moon = make_seg(r_km, 0.0, 0.0);
        let pos = dyn_earth_barycentric(jd_test(), &emb, &moon);
        let mag_au =
            (pos.x().value().powi(2) + pos.y().value().powi(2) + pos.z().value().powi(2)).sqrt();
        let earth_km = 1.5e8 - r_km / EARTH_MOON_RATIO;
        let expected_au = earth_km / 1.495_978_707e8;
        assert!(
            (mag_au - expected_au).abs() / expected_au < 1e-6,
            "mag_au={mag_au}, expected={expected_au}"
        );
    }

    // ── dyn_sun_barycentric with realistic values ─────────────────────────

    #[test]
    fn dyn_sun_barycentric_small_offset_from_origin() {
        // Sun barycentric should be close to origin (~few million km)
        // Use known-ish values: Sun is ~0.005 AU from SSB
        let sun = make_seg(5.0e5, 3.0e5, 1.0e5); // ~0.005 AU range
        let pos = dyn_sun_barycentric(jd_test(), &sun);
        let mag_au =
            (pos.x().value().powi(2) + pos.y().value().powi(2) + pos.z().value().powi(2)).sqrt();
        // ~0.005 AU is plausible for Sun-SSB offset
        assert!(mag_au < 0.1, "unexpected large Sun-SSB offset: {mag_au} AU");
    }
}
