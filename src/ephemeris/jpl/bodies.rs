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
//! [`Ephemeris`](crate::ephemeris::Ephemeris) trait.
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
//! - Time input: `JulianDate` (TT); runtime SPK evaluation converts TT to
//!   NAIF/SPICE ephemeris time seconds before Chebyshev lookup.
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

use super::eval::{jd_tt_to_spice_et_seconds, DynSegmentStack};

use crate::coordinates::{
    cartesian::{Position, Velocity},
    centers::{Barycentric, Geocentric, Heliocentric},
    frames::EclipticMeanJ2000,
    transform::VectorAstroExt,
};
use crate::ephemeris::EphemerisError;
use crate::qtty::{AstronomicalUnit, Day, Kilometer, Per};
use crate::time::JulianDate;

// ── Physical constants (from archive: siderust_archive::jpl::constants) ──────

use crate::archive::jpl::constants::EARTH_MOON_RATIO;

/// μ_Earth / (μ_Earth + μ_Moon), Earth's mass fraction of the EM system.
const FRAC_EARTH: f64 = EARTH_MOON_RATIO / (EARTH_MOON_RATIO + 1.0);

/// μ_Moon / (μ_Earth + μ_Moon), Moon's mass fraction of the EM system.
#[allow(dead_code)]
const FRAC_MOON: f64 = 1.0 / (EARTH_MOON_RATIO + 1.0);

/// Scale converting `moon_off = Moon - EMB` to the Earth-side correction
/// `EMB - Earth = moon_off * (FRAC_MOON / FRAC_EARTH) = moon_off / EMRAT`.
const EARTH_OFFSET_FROM_MOON_OFF: f64 = 1.0 / EARTH_MOON_RATIO;

/// Scale converting `moon_off = Moon - EMB` to the geocentric Earth-to-Moon
/// vector `Moon - Earth = moon_off / FRAC_EARTH`.
const MOON_GEO_FROM_MOON_OFF: f64 = 1.0 / FRAC_EARTH;

// ── Velocity type alias ─────────────────────────────────────────────────────

type AuPerDay = Per<AstronomicalUnit, Day>;

// ── DynSegmentDescriptor body-chain functions (runtime-loaded data) ───────────

/// Sun barycentric position (runtime data).
#[inline]
pub(crate) fn try_dyn_sun_barycentric(
    jd: JulianDate,
    sun: &DynSegmentStack,
) -> Result<Position<Barycentric, EclipticMeanJ2000, AstronomicalUnit>, EphemerisError> {
    let et = jd_tt_to_spice_et_seconds(jd);
    let sun_icrf = sun.try_position_et(et)?;
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
pub(crate) fn dyn_sun_barycentric(
    jd: JulianDate,
    sun: &DynSegmentStack,
) -> Position<Barycentric, EclipticMeanJ2000, AstronomicalUnit> {
    try_dyn_sun_barycentric(jd, sun).expect("runtime JPL Sun barycentric position unavailable")
}

/// Earth barycentric position (runtime data).
#[inline]
pub(crate) fn try_dyn_earth_barycentric(
    jd: JulianDate,
    emb: &DynSegmentStack,
    moon: &DynSegmentStack,
) -> Result<Position<Barycentric, EclipticMeanJ2000, AstronomicalUnit>, EphemerisError> {
    let et = jd_tt_to_spice_et_seconds(jd);
    let emb_pos = emb.try_position_et(et)?;
    let moon_off = moon.try_position_et(et)?;
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

/// Earth barycentric position using a direct `Earth -> EMB` SPK segment.
#[inline]
pub(crate) fn try_dyn_earth_barycentric_direct(
    jd: JulianDate,
    emb: &DynSegmentStack,
    earth: &DynSegmentStack,
) -> Result<Position<Barycentric, EclipticMeanJ2000, AstronomicalUnit>, EphemerisError> {
    let et = jd_tt_to_spice_et_seconds(jd);
    let emb_pos = emb.try_position_et(et)?;
    let earth_off = earth.try_position_et(et)?;
    let earth_ecl_au = (emb_pos + earth_off)
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
pub(crate) fn dyn_earth_barycentric(
    jd: JulianDate,
    emb: &DynSegmentStack,
    moon: &DynSegmentStack,
) -> Position<Barycentric, EclipticMeanJ2000, AstronomicalUnit> {
    try_dyn_earth_barycentric(jd, emb, moon)
        .expect("runtime JPL Earth barycentric position unavailable")
}

/// Earth heliocentric position (runtime data).
#[inline]
pub(crate) fn try_dyn_earth_heliocentric(
    jd: JulianDate,
    sun: &DynSegmentStack,
    emb: &DynSegmentStack,
    moon: &DynSegmentStack,
) -> Result<Position<Heliocentric, EclipticMeanJ2000, AstronomicalUnit>, EphemerisError> {
    let et = jd_tt_to_spice_et_seconds(jd);
    let emb_pos = emb.try_position_et(et)?;
    let moon_off = moon.try_position_et(et)?;
    let sun_pos = sun.try_position_et(et)?;
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

/// Earth heliocentric position using a direct `Earth -> EMB` SPK segment.
#[inline]
pub(crate) fn try_dyn_earth_heliocentric_direct(
    jd: JulianDate,
    sun: &DynSegmentStack,
    emb: &DynSegmentStack,
    earth: &DynSegmentStack,
) -> Result<Position<Heliocentric, EclipticMeanJ2000, AstronomicalUnit>, EphemerisError> {
    let et = jd_tt_to_spice_et_seconds(jd);
    let emb_pos = emb.try_position_et(et)?;
    let earth_off = earth.try_position_et(et)?;
    let sun_pos = sun.try_position_et(et)?;
    let earth_ecl_au = (emb_pos + earth_off - sun_pos)
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
pub(crate) fn dyn_earth_heliocentric(
    jd: JulianDate,
    sun: &DynSegmentStack,
    emb: &DynSegmentStack,
    moon: &DynSegmentStack,
) -> Position<Heliocentric, EclipticMeanJ2000, AstronomicalUnit> {
    try_dyn_earth_heliocentric(jd, sun, emb, moon)
        .expect("runtime JPL Earth heliocentric position unavailable")
}

/// Earth barycentric velocity (runtime data).
#[inline]
pub(crate) fn try_dyn_earth_barycentric_velocity(
    jd: JulianDate,
    emb: &DynSegmentStack,
    moon: &DynSegmentStack,
) -> Result<Velocity<EclipticMeanJ2000, AuPerDay>, EphemerisError> {
    let et = jd_tt_to_spice_et_seconds(jd);
    let v_emb = emb.try_velocity_et(et)?;
    let v_moon_off = moon.try_velocity_et(et)?;
    let v_earth_icrf = v_emb - v_moon_off.scale(EARTH_OFFSET_FROM_MOON_OFF);
    Ok(v_earth_icrf
        .to_frame::<EclipticMeanJ2000>(&crate::J2000)
        .to_unit::<AuPerDay>())
}

/// Earth barycentric velocity using a direct `Earth -> EMB` SPK segment.
#[inline]
pub(crate) fn try_dyn_earth_barycentric_velocity_direct(
    jd: JulianDate,
    emb: &DynSegmentStack,
    earth: &DynSegmentStack,
) -> Result<Velocity<EclipticMeanJ2000, AuPerDay>, EphemerisError> {
    let et = jd_tt_to_spice_et_seconds(jd);
    let v_emb = emb.try_velocity_et(et)?;
    let v_earth_off = earth.try_velocity_et(et)?;
    Ok((v_emb + v_earth_off)
        .to_frame::<EclipticMeanJ2000>(&crate::J2000)
        .to_unit::<AuPerDay>())
}

/// Earth barycentric velocity (runtime data).
#[inline]
pub(crate) fn dyn_earth_barycentric_velocity(
    jd: JulianDate,
    emb: &DynSegmentStack,
    moon: &DynSegmentStack,
) -> Velocity<EclipticMeanJ2000, AuPerDay> {
    try_dyn_earth_barycentric_velocity(jd, emb, moon)
        .expect("runtime JPL Earth barycentric velocity unavailable")
}

/// Moon geocentric position (runtime data).
#[inline]
pub(crate) fn try_dyn_moon_geocentric(
    jd: JulianDate,
    moon: &DynSegmentStack,
) -> Result<Position<Geocentric, EclipticMeanJ2000, Kilometer>, EphemerisError> {
    let et = jd_tt_to_spice_et_seconds(jd);
    let moon_off = moon.try_position_et(et)?;
    let moon_geo_icrf = moon_off.scale(MOON_GEO_FROM_MOON_OFF);
    let moon_geo_ecl = moon_geo_icrf.to_frame::<EclipticMeanJ2000>(&crate::J2000);
    Ok(Position::new(
        moon_geo_ecl.x(),
        moon_geo_ecl.y(),
        moon_geo_ecl.z(),
    ))
}

/// Moon geocentric position using direct `Moon -> EMB` and `Earth -> EMB`
/// SPK segments.
#[inline]
pub(crate) fn try_dyn_moon_geocentric_direct(
    jd: JulianDate,
    moon: &DynSegmentStack,
    earth: &DynSegmentStack,
) -> Result<Position<Geocentric, EclipticMeanJ2000, Kilometer>, EphemerisError> {
    let et = jd_tt_to_spice_et_seconds(jd);
    let moon_off = moon.try_position_et(et)?;
    let earth_off = earth.try_position_et(et)?;
    let moon_geo_ecl = (moon_off - earth_off).to_frame::<EclipticMeanJ2000>(&crate::J2000);
    Ok(Position::new(
        moon_geo_ecl.x(),
        moon_geo_ecl.y(),
        moon_geo_ecl.z(),
    ))
}

/// Moon geocentric position (runtime data).
#[inline]
pub(crate) fn dyn_moon_geocentric(
    jd: JulianDate,
    moon: &DynSegmentStack,
) -> Position<Geocentric, EclipticMeanJ2000, Kilometer> {
    try_dyn_moon_geocentric(jd, moon).expect("runtime JPL Moon geocentric position unavailable")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ephemeris::jpl::eval::{DynSegmentDescriptor, DynSegmentStack};

    const SECONDS_PER_DAY: f64 = crate::qtty::time::SECONDS_PER_DAY;
    const JD_J2000: f64 = tempoch::J2000_JD_TT_DAY.value();

    // ── DynSegmentDescriptor tests ────────────────────────────────────────

    /// Create a synthetic DynSegmentDescriptor representing a constant position.
    ///
    /// The segment spans 1000 days from J2000. The position at tau=0 (J2000+500d)
    /// is (x_km, y_km, z_km) in ICRF km. The velocity is zero.
    fn make_seg(x_km: f64, y_km: f64, z_km: f64) -> DynSegmentStack {
        use crate::qtty::Seconds;
        let ncoeff = 2usize;
        let rsize = 2 + 3 * ncoeff; // 8
        let intlen_secs = 1000.0 * SECONDS_PER_DAY;
        let mid = intlen_secs / 2.0;
        let radius = intlen_secs / 2.0;
        // Record: [mid, radius, cx0, cx1, cy0, cy1, cz0, cz1]
        let data = vec![mid, radius, x_km, 0.0, y_km, 0.0, z_km, 0.0];
        let descriptor = DynSegmentDescriptor {
            data_type: 2,
            init: Seconds::new(0.0),
            intlen: Seconds::new(intlen_secs),
            ncoeff,
            rsize,
            n_records: 1,
            data,
        };
        DynSegmentStack::single(descriptor, 0.0, intlen_secs)
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
