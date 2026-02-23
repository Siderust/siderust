// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Generic body-chain resolution for DE4xx ephemerides.
//!
//! ## Responsibility
//!
//! This module provides **generic** implementations of the body-chain arithmetic
//! required by the [`Ephemeris`](crate::calculus::ephemeris::Ephemeris) trait,
//! parameterized over the specific data source (DE440, DE441, etc.).
//!
//! It does **not** own:
//! - Physical constants — those come from `qtty` (single source of truth).
//! - Time-scale conversion — delegated to [`JulianDate::tt_to_tdb`].
//! - Frame-model definitions — delegated to the coordinate transform providers.
//!
//! ## Type safety
//!
//! All intermediate vectors carry their reference frame ([`ICRF`]) and unit
//! ([`Kilometer`] or [`Per<Kilometer, Day>`]) in the type system, so
//! body-chain arithmetic (subtraction, scaling) is checked at compile time.
//! Frame conversion (ICRF → EclipticMeanJ2000) and unit conversion (km → AU,
//! km/day → AU/day) are composed explicitly via transform/unit adapters.
//!
//! ## DE4xx segment semantics
//!
//! | Segment | NAIF IDs       | Meaning                                |
//! |---------|----------------|----------------------------------------|
//! | Sun     | 10 → 0 (SSB)  | Sun barycentric (ICRF, km)             |
//! | EMB     | 3  → 0 (SSB)  | Earth-Moon barycenter bary. (ICRF, km) |
//! | Moon    | 301 → 3 (EMB) | Moon offset from EMB (ICRF, km)        |
//!
//! Derived quantities:
//! - **Earth bary.**   = EMB − Moon_offset × μ_Moon / (μ_Earth + μ_Moon)
//! - **Moon geocentric**= Moon_offset × μ_Earth / (μ_Earth + μ_Moon)
//! - **Earth helio.**  = Earth_bary − Sun_bary

use super::eval::SegmentDescriptor;

use crate::coordinates::{
    cartesian::{Position, Velocity},
    centers::{Barycentric, Geocentric, Heliocentric},
    frames::EclipticMeanJ2000,
    transform::VectorAstroExt,
};
use crate::time::JulianDate;
use qtty::{AstronomicalUnit, Day, Kilometer, Per};

// ── Physical constants (single source of truth: qtty + DE4xx headers) ────

/// Earth/Moon mass ratio embedded in DE4xx headers.
///
/// This value is consistent across DE440 and DE441.
const EARTH_MOON_RATIO: f64 = 81.300_569_074_190_62;

/// μ_Earth / (μ_Earth + μ_Moon)  — Earth's mass fraction of the EM system.
const FRAC_EARTH: f64 = EARTH_MOON_RATIO / (EARTH_MOON_RATIO + 1.0);

/// μ_Moon / (μ_Earth + μ_Moon)  — Moon's mass fraction of the EM system.
const FRAC_MOON: f64 = 1.0 / (EARTH_MOON_RATIO + 1.0);

// ── Velocity type alias ─────────────────────────────────────────────────

type AuPerDay = Per<AstronomicalUnit, Day>;

// ── Generic body-chain functions ─────────────────────────────────────────

/// Sun barycentric position in EclipticMeanJ2000 J2000 (AU).
///
/// Generic over any DE4xx data source providing a SUN segment.
#[inline]
pub fn sun_barycentric(
    jd: JulianDate,
    sun: &SegmentDescriptor,
) -> Position<Barycentric, EclipticMeanJ2000, AstronomicalUnit> {
    let jd_tdb = JulianDate::tt_to_tdb(jd);
    let sun_icrf = sun.position(jd_tdb);
    let sun_ecl_au = sun_icrf
        .to_frame::<EclipticMeanJ2000>(&JulianDate::J2000)
        .to_unit::<AstronomicalUnit>();
    Position::new(sun_ecl_au.x(), sun_ecl_au.y(), sun_ecl_au.z())
}

/// Earth barycentric position in EclipticMeanJ2000 J2000 (AU).
///
/// `Earth_bary = EMB − Moon_offset × FRAC_MOON`
///
/// Generic over any DE4xx data source providing EMB and MOON segments.
#[inline]
pub fn earth_barycentric(
    jd: JulianDate,
    emb: &SegmentDescriptor,
    moon: &SegmentDescriptor,
) -> Position<Barycentric, EclipticMeanJ2000, AstronomicalUnit> {
    let jd_tdb = JulianDate::tt_to_tdb(jd);
    let emb_pos = emb.position(jd_tdb);
    let moon_off = moon.position(jd_tdb);
    let earth_icrf = emb_pos - moon_off.scale(FRAC_MOON);
    let earth_ecl_au = earth_icrf
        .to_frame::<EclipticMeanJ2000>(&JulianDate::J2000)
        .to_unit::<AstronomicalUnit>();
    Position::new(earth_ecl_au.x(), earth_ecl_au.y(), earth_ecl_au.z())
}

/// Earth heliocentric position in EclipticMeanJ2000 J2000 (AU).
///
/// `Earth_helio = Earth_bary − Sun_bary`
///
/// Generic over any DE4xx data source providing SUN, EMB, and MOON segments.
#[inline]
pub fn earth_heliocentric(
    jd: JulianDate,
    sun: &SegmentDescriptor,
    emb: &SegmentDescriptor,
    moon: &SegmentDescriptor,
) -> Position<Heliocentric, EclipticMeanJ2000, AstronomicalUnit> {
    let jd_tdb = JulianDate::tt_to_tdb(jd);
    let emb_pos = emb.position(jd_tdb);
    let moon_off = moon.position(jd_tdb);
    let sun_pos = sun.position(jd_tdb);
    let earth_icrf = emb_pos - moon_off.scale(FRAC_MOON) - sun_pos;
    let earth_ecl_au = earth_icrf
        .to_frame::<EclipticMeanJ2000>(&JulianDate::J2000)
        .to_unit::<AstronomicalUnit>();
    Position::new(earth_ecl_au.x(), earth_ecl_au.y(), earth_ecl_au.z())
}

/// Earth barycentric velocity in EclipticMeanJ2000 J2000 (AU/day).
///
/// `v_Earth = v_EMB − v_Moon_offset × FRAC_MOON`
///
/// Generic over any DE4xx data source providing EMB and MOON segments.
#[inline]
pub fn earth_barycentric_velocity(
    jd: JulianDate,
    emb: &SegmentDescriptor,
    moon: &SegmentDescriptor,
) -> Velocity<EclipticMeanJ2000, AuPerDay> {
    let jd_tdb = JulianDate::tt_to_tdb(jd);
    let v_emb = emb.velocity(jd_tdb);
    let v_moon_off = moon.velocity(jd_tdb);
    let v_earth_icrf = v_emb - v_moon_off.scale(FRAC_MOON);
    v_earth_icrf
        .to_frame::<EclipticMeanJ2000>(&JulianDate::J2000)
        .to_unit::<AuPerDay>()
}

/// Moon geocentric position in EclipticMeanJ2000 J2000 (km).
///
/// `Moon_geo = Moon_offset × FRAC_EARTH`
///
/// Generic over any DE4xx data source providing a MOON segment.
#[inline]
pub fn moon_geocentric(
    jd: JulianDate,
    moon: &SegmentDescriptor,
) -> Position<Geocentric, EclipticMeanJ2000, Kilometer> {
    let jd_tdb = JulianDate::tt_to_tdb(jd);
    let moon_off = moon.position(jd_tdb);
    let moon_geo_icrf = moon_off.scale(FRAC_EARTH);
    let moon_geo_ecl = moon_geo_icrf.to_frame::<EclipticMeanJ2000>(&JulianDate::J2000);
    Position::new(moon_geo_ecl.x(), moon_geo_ecl.y(), moon_geo_ecl.z())
}
