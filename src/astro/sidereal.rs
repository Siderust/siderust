// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Sidereal Time Module
//!
//! **Sidereal time** is the hour angle of the vernal equinox: a clock that tells us
//! "how far the Earth has rotated relative to the stars" instead of the Sun.  One
//! **mean sidereal day** is ≈ 23 h 56 m 4.09 s, i.e. about 0.99727 solar days.
//! Astronomers rely on it to aim equatorial‑mounted telescopes, reduce star‐track
//! images, or convert between Earth‐fixed and inertial coordinate frames.
//!
//! This module provides IAU 2006/2000A compliant sidereal time functions:
//!
//! | Function | Output | Notes |
//! |----------|--------|-------|
//! | [`gmst_iau2006`] | Greenwich Mean Sidereal Time (radians) | IAU 2006 ERA-based formula |
//! | [`gast_iau2006`] | Greenwich Apparent Sidereal Time (radians) | Includes equation of equinoxes |
//!
//! The implementation follows the IAU 2006 standard using the Earth Rotation Angle (ERA).
//! Accuracy is better than ±0.1″ for dates within ±100 years of J2000.
//!
//! ## Example
//! ```rust
//! use chrono::prelude::*;
//! use siderust::time::JulianDate;
//! use siderust::astro::sidereal::gmst_iau2006;
//! use qtty::*;
//!
//! let jd = JulianDate::from_utc(Utc::now());
//! let gmst = gmst_iau2006(jd, jd); // jd_ut1 ≈ jd_tt for most applications
//! let lon_madrid = Degrees::new(-3.7038).to::<Radian>();
//! let lst = gmst + lon_madrid;
//! println!("GMST = {:.4}°,  LST = {:.4}°", gmst.to::<Degree>(), lst.to::<Degree>());
//! ```

use crate::astro::era::earth_rotation_angle;
use crate::time::JulianDate;
use qtty::*;
use std::f64::consts::TAU;

/// Mean sidereal day length ≈ 0.9972696 solar days (23 h 56 m 4.09 s).
pub use qtty::time::SIDEREAL_DAY;

// ════════════════════════════════════════════════════════════════════════
// IAU 2006 ERA-based sidereal time
// ════════════════════════════════════════════════════════════════════════

/// Arcseconds-to-radians conversion factor.
const AS2RAD: f64 = std::f64::consts::PI / (180.0 * 3600.0);

/// Greenwich Mean Sidereal Time — IAU 2006 (ERA-based).
///
/// ```text
/// GMST(UT1, TT) = ERA(UT1) + polynomial(t)
/// ```
///
/// where the polynomial captures the accumulated precession of the equinox
/// relative to the CIO. The coefficients are from Capitaine et al. (2003),
/// adopted by IAU 2006 Resolution B1.
///
/// `jd_ut1`: Julian Date on the UT1 time scale (for ERA).
/// `jd_tt`:  Julian Date on the TT time scale (for the polynomial).
///
/// Returns GMST in **radians**, normalized to [0, 2π).
///
/// ## References
/// * Capitaine, Wallace & Chapront (2003), A&A 412, 567
/// * SOFA routine `iauGmst06`
#[inline]
pub fn gmst_iau2006(jd_ut1: JulianDate, jd_tt: JulianDate) -> Radians {
    let era = earth_rotation_angle(jd_ut1);
    let t = jd_tt.julian_centuries().value();

    // Polynomial: accumulated precession of equinox in arcseconds
    // Coefficients from Capitaine et al. (2003), eq. 42
    let poly_as = 0.014_506 + 4_612.156_534 * t + 1.391_581_7 * t.powi(2)
        - 0.000_000_44 * t.powi(3)
        - 0.000_029_956 * t.powi(4)
        - 0.000_000_036_8 * t.powi(5);

    let gmst = era.value() + poly_as * AS2RAD;
    Radians::new(gmst.rem_euclid(TAU))
}

/// Greenwich Apparent Sidereal Time — IAU 2006/2000A.
///
/// ```text
/// GAST = GMST + equation_of_the_equinoxes
/// ```
///
/// This adds the nutation correction (equation of the equinoxes) to the
/// mean sidereal time. The equation of the equinoxes is Δψ·cos(ε).
///
/// `jd_ut1`: Julian Date on the UT1 time scale.
/// `jd_tt`:  Julian Date on the TT time scale.
/// `dpsi`:   nutation in longitude (Δψ) in radians.
/// `true_obliquity`: ε_A + Δε in radians.
///
/// Returns GAST in **radians**, normalized to [0, 2π).
///
/// ## References
/// * SOFA routine `iauGst06`
#[inline]
pub fn gast_iau2006(
    jd_ut1: JulianDate,
    jd_tt: JulianDate,
    dpsi: Radians,
    true_obliquity: Radians,
) -> Radians {
    let gmst = gmst_iau2006(jd_ut1, jd_tt);
    let ee = crate::astro::era::equation_of_the_equinoxes(dpsi, true_obliquity);
    Radians::new((gmst.value() + ee.value()).rem_euclid(TAU))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gmst_iau2006_at_j2000() {
        // At J2000.0, GMST ≈ ERA ≈ 280.46° (polynomial term is tiny at t=0)
        let jd = JulianDate::J2000;
        let gmst = gmst_iau2006(jd, jd);
        let gmst_deg = gmst.to::<Degree>();
        assert!(
            (gmst_deg.value() - 280.46_f64).abs() < 0.1,
            "GMST at J2000 = {}°, expected ≈ 280.46°",
            gmst_deg
        );
    }

    #[test]
    fn gmst_iau2006_range() {
        for jd in [2_451_545.0, 2_460_000.5, 2_440_000.0] {
            let jd = JulianDate::new(jd);
            let gmst = gmst_iau2006(jd, jd);
            assert!(gmst >= Radians::new(0.0), "GMST should be ≥ 0");
            assert!(gmst < Radians::new(TAU), "GMST should be < 2π");
        }
    }

    #[test]
    fn gast_iau2006_close_to_gmst() {
        // GAST = GMST + Δψ·cos(ε). For small nutation (~17″ max), |GAST−GMST| < 20″.
        let jd = JulianDate::new(2_460_000.5);
        let nutation = crate::astro::nutation::nutation_iau2000b(jd);
        let true_obliquity = nutation.true_obliquity();
        let gast = gast_iau2006(jd, jd, nutation.dpsi, true_obliquity);
        let gmst = gmst_iau2006(jd, jd);

        let diff_as = (gast - gmst).value().abs() * 206_264.806;
        assert!(
            diff_as < 20.0,
            "GAST–GMST = {:.3}″, max nutation effect should be < 20″",
            diff_as
        );
    }
}
