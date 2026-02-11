// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Body chain resolution for DE440 ephemeris.
//!
//! DE440 segments give:
//! - Sun  (NAIF 10) relative to SSB  (NAIF 0)   → Sun barycentric (ICRF, km)
//! - EMB  (NAIF 3)  relative to SSB  (NAIF 0)   → EMB barycentric (ICRF, km)
//! - Moon (NAIF 301) relative to EMB (NAIF 3)    → Moon–EMB offset  (ICRF, km)
//!
//! From these we derive:
//! - **Earth barycentric**   = EMB − Moon_offset × μ_Moon / (μ_Earth + μ_Moon)
//! - **Moon geocentric**     = Moon_offset × μ_Earth / (μ_Earth + μ_Moon)
//! - **Earth heliocentric**  = Earth_bary − Sun_bary
//!
//! All outputs are rotated from ICRF (equatorial) → Ecliptic J2000 and
//! converted from km → AU (positions) or km/day → AU/day (velocities).

use super::data;
use super::eval;

use crate::coordinates::{
    cartesian::{Position, Velocity},
    centers::{Barycentric, Geocentric, Heliocentric},
    frames::Ecliptic,
};
use crate::targets::Target;
use crate::time::JulianDate;
use qtty::{AstronomicalUnit, AstronomicalUnits, Day, Kilometer, Kilometers};

/// AU in km (IAU 2012 exact).
const KM_PER_AU: f64 = 149_597_870.7;

/// Earth/Moon mass ratio from DE440.
const EARTH_MOON_RATIO: f64 = 81.300_569_074_190_62;

/// μ_Earth / (μ_Earth + μ_Moon)
const FRAC_EARTH: f64 = EARTH_MOON_RATIO / (EARTH_MOON_RATIO + 1.0);

/// μ_Moon / (μ_Earth + μ_Moon)
const FRAC_MOON: f64 = 1.0 / (EARTH_MOON_RATIO + 1.0);

/// J2000 mean obliquity (IAU 2006): 84381.406″ → radians.
const OBLIQUITY_RAD: f64 = 84381.406 / 3600.0 * std::f64::consts::PI / 180.0;

// Pre-computed sin/cos of obliquity for ICRF→Ecliptic rotation.
// (Computed once at first use via std::sync::LazyLock for const-correctness.)
//
// The rotation from equatorial (ICRF) to ecliptic J2000 about the X axis by ε:
//   x_ecl =  x_eq
//   y_ecl =  y_eq · cos ε + z_eq · sin ε
//   z_ecl = -y_eq · sin ε + z_eq · cos ε

/// Rotate an [x, y, z] vector from ICRF (equatorial) to Ecliptic J2000.
#[inline]
fn icrf_to_ecliptic(v: [f64; 3]) -> [f64; 3] {
    let (sin_e, cos_e) = OBLIQUITY_RAD.sin_cos();
    [
        v[0],
        v[1] * cos_e + v[2] * sin_e,
        -v[1] * sin_e + v[2] * cos_e,
    ]
}

/// Convert TT Julian Date → TDB Julian Date (Fairhead & Bretagnon).
#[inline]
fn tt_to_tdb(jd_tt: f64) -> f64 {
    let j2000 = 2_451_545.0;
    let e = (357.53 + 0.985_600_28 * (jd_tt - j2000)).to_radians();
    let delta_days = (1.658e-3 * e.sin() + 1.4e-6 * (2.0 * e).sin()) / 86_400.0;
    jd_tt + delta_days
}

// ── Velocity type alias (same as in ephemeris/mod.rs) ────────────────────

type AuPerDay = qtty::Per<AstronomicalUnit, Day>;

// ── Segment accessors ────────────────────────────────────────────────────
//
// Each body module in `data` exposes: INIT, INTLEN, NCOEFF, RSIZE, N_RECORDS,
// and `fn record(i: usize) -> &'static [f64]`.

/// Evaluate Sun barycentric position (ICRF, km) at JD(TDB).
fn sun_bary_icrf(jd_tdb: f64) -> [f64; 3] {
    eval::position(
        jd_tdb,
        data::sun::INIT,
        data::sun::INTLEN,
        data::sun::NCOEFF,
        data::sun::RSIZE,
        data::sun::N_RECORDS,
        data::sun::record,
    )
}

/// Evaluate EMB barycentric position (ICRF, km) at JD(TDB).
fn emb_bary_icrf(jd_tdb: f64) -> [f64; 3] {
    eval::position(
        jd_tdb,
        data::emb::INIT,
        data::emb::INTLEN,
        data::emb::NCOEFF,
        data::emb::RSIZE,
        data::emb::N_RECORDS,
        data::emb::record,
    )
}

/// Evaluate Moon offset from EMB (ICRF, km) at JD(TDB).
fn moon_emb_icrf(jd_tdb: f64) -> [f64; 3] {
    eval::position(
        jd_tdb,
        data::moon::INIT,
        data::moon::INTLEN,
        data::moon::NCOEFF,
        data::moon::RSIZE,
        data::moon::N_RECORDS,
        data::moon::record,
    )
}

/// Evaluate EMB barycentric velocity (ICRF, km/day) at JD(TDB).
fn emb_bary_vel_icrf(jd_tdb: f64) -> [f64; 3] {
    eval::velocity(
        jd_tdb,
        data::emb::INIT,
        data::emb::INTLEN,
        data::emb::NCOEFF,
        data::emb::RSIZE,
        data::emb::N_RECORDS,
        data::emb::record,
    )
}

/// Evaluate Moon offset from EMB velocity (ICRF, km/day) at JD(TDB).
fn moon_emb_vel_icrf(jd_tdb: f64) -> [f64; 3] {
    eval::velocity(
        jd_tdb,
        data::moon::INIT,
        data::moon::INTLEN,
        data::moon::NCOEFF,
        data::moon::RSIZE,
        data::moon::N_RECORDS,
        data::moon::record,
    )
}

#[allow(dead_code)]
/// Evaluate Sun barycentric velocity (ICRF, km/day) at JD(TDB).
fn sun_bary_vel_icrf(jd_tdb: f64) -> [f64; 3] {
    eval::velocity(
        jd_tdb,
        data::sun::INIT,
        data::sun::INTLEN,
        data::sun::NCOEFF,
        data::sun::RSIZE,
        data::sun::N_RECORDS,
        data::sun::record,
    )
}

// ── Public API ───────────────────────────────────────────────────────────

/// Sun barycentric position in Ecliptic J2000 (AU).
pub fn sun_barycentric(jd: JulianDate) -> Target<Position<Barycentric, Ecliptic, AstronomicalUnit>>
{
    let jd_tdb = tt_to_tdb(jd.value());
    let pos_km = icrf_to_ecliptic(sun_bary_icrf(jd_tdb));
    Target::new_static(
        Position::new(
            AstronomicalUnits::new(pos_km[0] / KM_PER_AU),
            AstronomicalUnits::new(pos_km[1] / KM_PER_AU),
            AstronomicalUnits::new(pos_km[2] / KM_PER_AU),
        ),
        jd,
    )
}

/// Earth barycentric position in Ecliptic J2000 (AU).
///
/// Earth_bary = EMB_bary − Moon_offset × μ_Moon / (μ_Earth + μ_Moon)
pub fn earth_barycentric(
    jd: JulianDate,
) -> Target<Position<Barycentric, Ecliptic, AstronomicalUnit>> {
    let jd_tdb = tt_to_tdb(jd.value());
    let emb = emb_bary_icrf(jd_tdb);
    let moon_off = moon_emb_icrf(jd_tdb);
    let earth_icrf = [
        emb[0] - moon_off[0] * FRAC_MOON,
        emb[1] - moon_off[1] * FRAC_MOON,
        emb[2] - moon_off[2] * FRAC_MOON,
    ];
    let pos_ecl = icrf_to_ecliptic(earth_icrf);
    Target::new_static(
        Position::new(
            AstronomicalUnits::new(pos_ecl[0] / KM_PER_AU),
            AstronomicalUnits::new(pos_ecl[1] / KM_PER_AU),
            AstronomicalUnits::new(pos_ecl[2] / KM_PER_AU),
        ),
        jd,
    )
}

/// Earth heliocentric position in Ecliptic J2000 (AU).
///
/// Earth_helio = Earth_bary − Sun_bary
pub fn earth_heliocentric(
    jd: JulianDate,
) -> Target<Position<Heliocentric, Ecliptic, AstronomicalUnit>> {
    let jd_tdb = tt_to_tdb(jd.value());
    let emb = emb_bary_icrf(jd_tdb);
    let moon_off = moon_emb_icrf(jd_tdb);
    let sun = sun_bary_icrf(jd_tdb);
    let earth_icrf = [
        emb[0] - moon_off[0] * FRAC_MOON - sun[0],
        emb[1] - moon_off[1] * FRAC_MOON - sun[1],
        emb[2] - moon_off[2] * FRAC_MOON - sun[2],
    ];
    let pos_ecl = icrf_to_ecliptic(earth_icrf);
    Target::new_static(
        Position::new(
            AstronomicalUnits::new(pos_ecl[0] / KM_PER_AU),
            AstronomicalUnits::new(pos_ecl[1] / KM_PER_AU),
            AstronomicalUnits::new(pos_ecl[2] / KM_PER_AU),
        ),
        jd,
    )
}

/// Earth barycentric velocity in Ecliptic J2000 (AU/day).
///
/// v_Earth_bary = v_EMB_bary − v_Moon_offset × μ_Moon / (μ_Earth + μ_Moon)
pub fn earth_barycentric_velocity(jd: JulianDate) -> Velocity<Ecliptic, AuPerDay> {
    let jd_tdb = tt_to_tdb(jd.value());
    let v_emb = emb_bary_vel_icrf(jd_tdb);
    let v_moon_off = moon_emb_vel_icrf(jd_tdb);
    let v_earth_icrf = [
        v_emb[0] - v_moon_off[0] * FRAC_MOON,
        v_emb[1] - v_moon_off[1] * FRAC_MOON,
        v_emb[2] - v_moon_off[2] * FRAC_MOON,
    ];
    let v_ecl = icrf_to_ecliptic(v_earth_icrf);
    // eval::velocity returns km/day; convert to AU/day
    Velocity::new(
        qtty::velocity::Velocity::<AstronomicalUnit, Day>::new(v_ecl[0] / KM_PER_AU),
        qtty::velocity::Velocity::<AstronomicalUnit, Day>::new(v_ecl[1] / KM_PER_AU),
        qtty::velocity::Velocity::<AstronomicalUnit, Day>::new(v_ecl[2] / KM_PER_AU),
    )
}

/// Moon geocentric position in Ecliptic J2000 (km).
///
/// Moon_geo = Moon_offset × μ_Earth / (μ_Earth + μ_Moon)
///
/// Note: This gives the Moon's position relative to Earth's center,
/// expressed in ecliptic coordinates, matching the signature expected
/// by the `Ephemeris` trait.
pub fn moon_geocentric(jd: JulianDate) -> Position<Geocentric, Ecliptic, Kilometer> {
    let jd_tdb = tt_to_tdb(jd.value());
    let moon_off = moon_emb_icrf(jd_tdb);
    let moon_geo_icrf = [
        moon_off[0] * FRAC_EARTH,
        moon_off[1] * FRAC_EARTH,
        moon_off[2] * FRAC_EARTH,
    ];
    let pos_ecl = icrf_to_ecliptic(moon_geo_icrf);
    Position::new(
        Kilometers::new(pos_ecl[0]),
        Kilometers::new(pos_ecl[1]),
        Kilometers::new(pos_ecl[2]),
    )
}
