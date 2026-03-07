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

use super::eval::{DynSegmentDescriptor, SegmentDescriptor};

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

// ═══════════════════════════════════════════════════════════════════════════
// DynSegmentDescriptor variants — runtime-loaded data
// ═══════════════════════════════════════════════════════════════════════════

/// Sun barycentric position (runtime data).
#[inline]
pub fn dyn_sun_barycentric(
    jd: JulianDate,
    sun: &DynSegmentDescriptor,
) -> Position<Barycentric, EclipticMeanJ2000, AstronomicalUnit> {
    let jd_tdb = JulianDate::tt_to_tdb(jd);
    let sun_icrf = sun.position(jd_tdb);
    let sun_ecl_au = sun_icrf
        .to_frame::<EclipticMeanJ2000>(&JulianDate::J2000)
        .to_unit::<AstronomicalUnit>();
    Position::new(sun_ecl_au.x(), sun_ecl_au.y(), sun_ecl_au.z())
}

/// Earth barycentric position (runtime data).
#[inline]
pub fn dyn_earth_barycentric(
    jd: JulianDate,
    emb: &DynSegmentDescriptor,
    moon: &DynSegmentDescriptor,
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

/// Earth heliocentric position (runtime data).
#[inline]
pub fn dyn_earth_heliocentric(
    jd: JulianDate,
    sun: &DynSegmentDescriptor,
    emb: &DynSegmentDescriptor,
    moon: &DynSegmentDescriptor,
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

/// Earth barycentric velocity (runtime data).
#[inline]
pub fn dyn_earth_barycentric_velocity(
    jd: JulianDate,
    emb: &DynSegmentDescriptor,
    moon: &DynSegmentDescriptor,
) -> Velocity<EclipticMeanJ2000, AuPerDay> {
    let jd_tdb = JulianDate::tt_to_tdb(jd);
    let v_emb = emb.velocity(jd_tdb);
    let v_moon_off = moon.velocity(jd_tdb);
    let v_earth_icrf = v_emb - v_moon_off.scale(FRAC_MOON);
    v_earth_icrf
        .to_frame::<EclipticMeanJ2000>(&JulianDate::J2000)
        .to_unit::<AuPerDay>()
}

/// Moon geocentric position (runtime data).
#[inline]
pub fn dyn_moon_geocentric(
    jd: JulianDate,
    moon: &DynSegmentDescriptor,
) -> Position<Geocentric, EclipticMeanJ2000, Kilometer> {
    let jd_tdb = JulianDate::tt_to_tdb(jd);
    let moon_off = moon.position(jd_tdb);
    let moon_geo_icrf = moon_off.scale(FRAC_EARTH);
    let moon_geo_ecl = moon_geo_icrf.to_frame::<EclipticMeanJ2000>(&JulianDate::J2000);
    Position::new(moon_geo_ecl.x(), moon_geo_ecl.y(), moon_geo_ecl.z())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::calculus::jpl::eval::{DynSegmentDescriptor, SegmentDescriptor};

    const SECONDS_PER_DAY: f64 = 86400.0;
    const JD_J2000: f64 = 2451545.0;

    // ── Shared test data for SegmentDescriptor (compile-time variant) ─────

    /// Static record for sun-like segment: position ~0.5 AU at tau=0.
    static SUN_RECORD: [f64; 8] = [
        500.0 * 86400.0, // mid (seconds past J2000)
        500.0 * 86400.0, // radius = half-interval
        7.5e7,
        0.0, // cx ~0.5 AU in km
        3.0e7,
        0.0, // cy
        1.0e7,
        0.0, // cz
    ];

    /// Static record for EMB-like segment: ~1 AU from SSB.
    static EMB_RECORD: [f64; 8] = [
        500.0 * 86400.0,
        500.0 * 86400.0,
        1.5e8,
        0.0, // ~1 AU
        0.0,
        0.0,
        0.0,
        0.0,
    ];

    /// Static record for Moon offset (~Earth-Moon distance).
    static MOON_RECORD: [f64; 8] = [
        500.0 * 86400.0,
        500.0 * 86400.0,
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
        use qtty::Seconds;
        SegmentDescriptor {
            init: Seconds::new(0.0),
            intlen: Seconds::new(1000.0 * SECONDS_PER_DAY),
            ncoeff: 2,
            n_records: 1,
            record_fn,
        }
    }

    fn jd_test_static() -> JulianDate {
        JulianDate::new(JD_J2000 + 500.0)
    }

    // ── SegmentDescriptor (compile-time) tests ────────────────────────────

    #[test]
    fn sun_barycentric_is_finite() {
        let sun = make_static_seg(sun_record_fn);
        let pos = sun_barycentric(jd_test_static(), &sun);
        assert!(pos.x().value().is_finite());
        assert!(pos.y().value().is_finite());
        assert!(pos.z().value().is_finite());
    }

    #[test]
    fn earth_barycentric_is_finite() {
        let emb = make_static_seg(emb_record_fn);
        let moon = make_static_seg(moon_record_fn);
        let pos = earth_barycentric(jd_test_static(), &emb, &moon);
        assert!(pos.x().value().is_finite());
        assert!(pos.y().value().is_finite());
        assert!(pos.z().value().is_finite());
    }

    #[test]
    fn earth_heliocentric_is_finite() {
        let sun = make_static_seg(sun_record_fn);
        let emb = make_static_seg(emb_record_fn);
        let moon = make_static_seg(moon_record_fn);
        let pos = earth_heliocentric(jd_test_static(), &sun, &emb, &moon);
        assert!(pos.x().value().is_finite());
        assert!(pos.y().value().is_finite());
        assert!(pos.z().value().is_finite());
    }

    #[test]
    fn earth_barycentric_velocity_is_finite() {
        let emb = make_static_seg(emb_record_fn);
        let moon = make_static_seg(moon_record_fn);
        let vel = earth_barycentric_velocity(jd_test_static(), &emb, &moon);
        assert!(vel.x().value().is_finite());
        assert!(vel.y().value().is_finite());
        assert!(vel.z().value().is_finite());
    }

    #[test]
    fn moon_geocentric_is_finite() {
        let moon = make_static_seg(moon_record_fn);
        let pos = moon_geocentric(jd_test_static(), &moon);
        assert!(pos.x().value().is_finite());
        assert!(pos.y().value().is_finite());
        assert!(pos.z().value().is_finite());
    }

    #[test]
    fn moon_geocentric_magnitude_scaled_by_frac_earth() {
        let r_km = 384400.0;
        let moon = make_static_seg(moon_record_fn); // position = (r_km, 0, 0) at tau=0
        let pos = moon_geocentric(jd_test_static(), &moon);
        let mag =
            (pos.x().value().powi(2) + pos.y().value().powi(2) + pos.z().value().powi(2)).sqrt();
        let expected = r_km * FRAC_EARTH;
        assert!(
            (mag - expected).abs() / expected < 0.01,
            "mag={mag}, expected={expected}"
        );
    }

    // ── DynSegmentDescriptor tests ────────────────────────────────────────

    /// Create a synthetic DynSegmentDescriptor representing a constant position.
    ///
    /// The segment spans 1000 days from J2000. The position at tau=0 (J2000+500d)
    /// is (x_km, y_km, z_km) in ICRF km. The velocity is zero.
    fn make_seg(x_km: f64, y_km: f64, z_km: f64) -> DynSegmentDescriptor {
        use qtty::Seconds;
        let ncoeff = 2usize;
        let rsize = 2 + 3 * ncoeff; // 8
        let intlen_secs = 1000.0 * SECONDS_PER_DAY;
        let mid = intlen_secs / 2.0;
        let radius = intlen_secs / 2.0;
        // Record: [mid, radius, cx0, cx1, cy0, cy1, cz0, cz1]
        let data = vec![mid, radius, x_km, 0.0, y_km, 0.0, z_km, 0.0];
        DynSegmentDescriptor {
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
        JulianDate::new(JD_J2000 + 500.0)
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
        assert!(pos.x().value().is_finite());
        assert!(pos.y().value().is_finite());
        assert!(pos.z().value().is_finite());
    }

    // ── dyn_earth_barycentric ─────────────────────────────────────────────

    #[test]
    fn dyn_earth_barycentric_is_finite() {
        let emb = make_seg(1.5e8, 0.0, 0.0);
        let moon = make_seg(3.84e5, 0.0, 0.0);
        let pos = dyn_earth_barycentric(jd_test(), &emb, &moon);
        assert!(pos.x().value().is_finite());
        assert!(pos.y().value().is_finite());
        assert!(pos.z().value().is_finite());
    }

    // ── dyn_earth_heliocentric ─────────────────────────────────────────────

    #[test]
    fn dyn_earth_heliocentric_is_finite() {
        let sun = make_seg(1.0e8, 0.0, 0.0);
        let emb = make_seg(1.5e8, 0.0, 0.0);
        let moon = make_seg(3.84e5, 0.0, 0.0);
        let pos = dyn_earth_heliocentric(jd_test(), &sun, &emb, &moon);
        assert!(pos.x().value().is_finite());
        assert!(pos.y().value().is_finite());
        assert!(pos.z().value().is_finite());
    }

    // ── dyn_earth_barycentric_velocity ────────────────────────────────────

    #[test]
    fn dyn_earth_barycentric_velocity_is_finite() {
        let emb = make_seg(1.5e8, 0.0, 0.0);
        let moon = make_seg(3.84e5, 0.0, 0.0);
        let vel = dyn_earth_barycentric_velocity(jd_test(), &emb, &moon);
        assert!(vel.x().value().is_finite());
        assert!(vel.y().value().is_finite());
        assert!(vel.z().value().is_finite());
    }

    // ── dyn_moon_geocentric ───────────────────────────────────────────────

    #[test]
    fn dyn_moon_geocentric_is_finite() {
        let moon = make_seg(3.84e5, 1.0e4, 2.0e3);
        let pos = dyn_moon_geocentric(jd_test(), &moon);
        assert!(pos.x().value().is_finite());
        assert!(pos.y().value().is_finite());
        assert!(pos.z().value().is_finite());
    }

    #[test]
    fn dyn_moon_geocentric_scaled_by_frac_earth() {
        // With moon_off = (R, 0, 0) at tau=0, geocentric = moon_off * FRAC_EARTH
        // After frame transform, x should be close to R * FRAC_EARTH (ignoring rotation)
        let r_km = 384400.0;
        let moon = make_seg(r_km, 0.0, 0.0);
        let pos = dyn_moon_geocentric(jd_test(), &moon);
        // The frame rotation can change x/y/z but the magnitude should be ~r_km * FRAC_EARTH
        let magnitude =
            (pos.x().value().powi(2) + pos.y().value().powi(2) + pos.z().value().powi(2)).sqrt();
        let expected = r_km * FRAC_EARTH;
        assert!(
            (magnitude - expected).abs() / expected < 0.01,
            "magnitude={magnitude}, expected={expected}"
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
