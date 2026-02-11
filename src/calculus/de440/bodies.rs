// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Body-chain resolution for DE440 ephemeris.
//!
//! ## Responsibility
//!
//! This module is **only** responsible for combining the raw ICRF segment
//! outputs (km) into the composite body states required by the
//! [`Ephemeris`](crate::calculus::ephemeris::Ephemeris) trait.
//!
//! It does **not** own:
//! - Physical constants — those come from `qtty` (single source of truth).
//! - Time-scale conversion — delegated to [`JulianDate::tt_to_tdb`].
//! - Frame rotation — uses [`Rotation3::from_x_rotation`] with the
//!   J2000 mean obliquity to rotate ICRF → Ecliptic.
//!
//! ## Type safety
//!
//! All intermediate vectors carry their reference frame ([`ICRF`]) and unit
//! ([`Kilometer`] or [`Per<Kilometer, Day>`]) in the type system, so
//! body-chain arithmetic (subtraction, scaling) is checked at compile time.
//! Frame conversion (ICRF → Ecliptic) happens through an explicit rotation
//! step that produces vectors / positions in the target frame.
//!
//! ## DE440 segment semantics
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

use super::data;

use crate::coordinates::{
    cartesian::{Position, Velocity},
    centers::{Barycentric, Geocentric, Heliocentric},
    frames::{Ecliptic, ICRF},
};
use crate::targets::Target;
use crate::time::JulianDate;
use affn::{Displacement, Rotation3};
use qtty::{AstronomicalUnit, Day, Kilometer, Per};

// ── Physical constants (single source of truth: qtty + DE440 header) ─────

/// Earth/Moon mass ratio embedded in the DE440 header.
///
/// This value is specific to DE440 and is not a general physical constant.
const EARTH_MOON_RATIO: f64 = 81.300_569_074_190_62;

/// μ_Earth / (μ_Earth + μ_Moon)  — Earth's mass fraction of the EM system.
const FRAC_EARTH: f64 = EARTH_MOON_RATIO / (EARTH_MOON_RATIO + 1.0);

/// μ_Moon / (μ_Earth + μ_Moon)  — Moon's mass fraction of the EM system.
const FRAC_MOON: f64 = 1.0 / (EARTH_MOON_RATIO + 1.0);

/// J2000 mean obliquity (IAU 2006): 84 381.406″ expressed as radians via qtty.
///
/// Using `Arcseconds → Radian` conversion through unit ratios ensures the
/// value stays consistent with the rest of the codebase.
const OBLIQUITY_RAD: qtty::Radians = qtty::Arcseconds::new(84_381.406).to_const::<qtty::Radian>();

// ── ICRF → Ecliptic rotation ─────────────────────────────────────────────

/// Pre-computed rotation matrix from ICRF to Ecliptic J2000.
///
/// Pure rotation about the +X axis by **−ε** (negative obliquity), because
/// going from equatorial to ecliptic tilts the reference plane downward
/// relative to the right-hand sense around X.
fn icrf_to_ecliptic_rotation() -> Rotation3 {
    Rotation3::from_x_rotation(-OBLIQUITY_RAD.value())
}

// ── Velocity type alias ─────────────────────────────────────────────────

type AuPerDay = Per<AstronomicalUnit, Day>;

// ── Typed conversion helpers ────────────────────────────────────────────

/// Rotate a displacement from ICRF to Ecliptic and convert km → AU,
/// wrapping the result in a typed `Position`.
#[inline]
fn icrf_km_to_ecliptic_au_position<
    C: crate::coordinates::centers::ReferenceCenter<Params = ()>,
>(
    v: Displacement<ICRF, Kilometer>,
) -> Position<C, Ecliptic, AstronomicalUnit> {
    let rot = icrf_to_ecliptic_rotation();
    let [ex, ey, ez] = rot * [v.x(), v.y(), v.z()];
    Position::new(
        ex.to::<AstronomicalUnit>(),
        ey.to::<AstronomicalUnit>(),
        ez.to::<AstronomicalUnit>(),
    )
}

/// Rotate a velocity from ICRF to Ecliptic and convert km/day → AU/day.
#[inline]
fn icrf_to_ecliptic_velocity(
    v: affn::Velocity<ICRF, Per<Kilometer, Day>>,
) -> Velocity<Ecliptic, AuPerDay> {
    let rot = icrf_to_ecliptic_rotation();
    let [ex, ey, ez] = rot * [v.x(), v.y(), v.z()];
    Velocity::new(
        ex.to::<AuPerDay>(),
        ey.to::<AuPerDay>(),
        ez.to::<AuPerDay>(),
    )
}

// ── Public API ───────────────────────────────────────────────────────────

/// Sun barycentric position in Ecliptic J2000 (AU).
pub fn sun_barycentric(jd: JulianDate) -> Target<Position<Barycentric, Ecliptic, AstronomicalUnit>>
{
    let jd_tdb = JulianDate::tt_to_tdb(jd);
    let sun_icrf = data::SUN.position(jd_tdb);
    Target::new_static(icrf_km_to_ecliptic_au_position(sun_icrf), jd)
}

/// Earth barycentric position in Ecliptic J2000 (AU).
///
/// `Earth_bary = EMB − Moon_offset × FRAC_MOON`
pub fn earth_barycentric(
    jd: JulianDate,
) -> Target<Position<Barycentric, Ecliptic, AstronomicalUnit>> {
    let jd_tdb = JulianDate::tt_to_tdb(jd);
    let emb = data::EMB.position(jd_tdb);
    let moon_off = data::MOON.position(jd_tdb);
    let earth_icrf = emb - moon_off.scale(FRAC_MOON);
    Target::new_static(icrf_km_to_ecliptic_au_position(earth_icrf), jd)
}

/// Earth heliocentric position in Ecliptic J2000 (AU).
///
/// `Earth_helio = Earth_bary − Sun_bary`
pub fn earth_heliocentric(
    jd: JulianDate,
) -> Target<Position<Heliocentric, Ecliptic, AstronomicalUnit>> {
    let jd_tdb = JulianDate::tt_to_tdb(jd);
    let emb = data::EMB.position(jd_tdb);
    let moon_off = data::MOON.position(jd_tdb);
    let sun = data::SUN.position(jd_tdb);
    let earth_icrf = emb - moon_off.scale(FRAC_MOON) - sun;
    Target::new_static(icrf_km_to_ecliptic_au_position(earth_icrf), jd)
}

/// Earth barycentric velocity in Ecliptic J2000 (AU/day).
///
/// `v_Earth = v_EMB − v_Moon_offset × FRAC_MOON`
pub fn earth_barycentric_velocity(jd: JulianDate) -> Velocity<Ecliptic, AuPerDay> {
    let jd_tdb = JulianDate::tt_to_tdb(jd);
    let v_emb = data::EMB.velocity(jd_tdb);
    let v_moon_off = data::MOON.velocity(jd_tdb);
    let v_earth_icrf = v_emb - v_moon_off.scale(FRAC_MOON);
    icrf_to_ecliptic_velocity(v_earth_icrf)
}

/// Moon geocentric position in Ecliptic J2000 (km).
///
/// `Moon_geo = Moon_offset × FRAC_EARTH`
pub fn moon_geocentric(jd: JulianDate) -> Position<Geocentric, Ecliptic, Kilometer> {
    let jd_tdb = JulianDate::tt_to_tdb(jd);
    let moon_off = data::MOON.position(jd_tdb);
    let moon_geo_icrf = moon_off.scale(FRAC_EARTH);
    let rot = icrf_to_ecliptic_rotation();
    let [ex, ey, ez] = rot * [moon_geo_icrf.x(), moon_geo_icrf.y(), moon_geo_icrf.z()];
    Position::new(ex, ey, ez)
}
