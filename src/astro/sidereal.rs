// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Sidereal Time Module
//!
//! IAU 2006/2000A compliant Greenwich Mean and Apparent Sidereal Time
//! functions, suitable for converting between Earth-fixed and inertial
//! frames and for pointing equatorially mounted telescopes.
//!
//! ## Scientific scope
//!
//! Sidereal time is the hour angle of the vernal equinox: a clock that
//! tracks how far the Earth has rotated relative to the stars rather than
//! the Sun. One mean sidereal day is ≈ 23 h 56 m 4.09 s (≈ 0.99727 solar
//! days). Greenwich Mean Sidereal Time (GMST) measures rotation against the
//! mean equinox; Greenwich Apparent Sidereal Time (GAST) additionally
//! includes the equation of the equinoxes (the projection of nutation onto
//! the equator). Local Sidereal Time (LST) is GMST/GAST plus the observer's
//! geodetic longitude.
//!
//! ## Technical scope
//!
//! [`gmst_iau2006`] implements the ERA-based IAU 2006 formulation
//! `GMST(UT1, TT) = ERA(UT1) + polynomial(t)`, with the polynomial
//! coefficients of Capitaine, Wallace & Chapront (2003) adopted by IAU 2006
//! Resolution B1. The two time arguments are kept distinct (UT1 for ERA,
//! TT for the polynomial) to match SOFA semantics. [`gast_iau2006`] adds
//! the equation of the equinoxes built from the nutation model and CIO
//! locator. Accuracy is better than ±0.1″ for dates within ±100 yr of
//! J2000.
//!
//! ## Example
//!
//! ```rust
//! use chrono::prelude::*;
//! use siderust::time::JulianDate;
//! use siderust::astro::sidereal::gmst_iau2006;
//! use siderust::qtty::*;
//!
//! let jd = JulianDate::from_chrono(Utc::now());
//! let gmst = gmst_iau2006(jd, jd); // jd_ut1 ≈ jd_tt for most applications
//! let lon_madrid = Degrees::new(-3.7038).to::<Radian>();
//! let lst = gmst + lon_madrid;
//! println!("GMST = {:.4}°,  LST = {:.4}°", gmst.to::<Degree>(), lst.to::<Degree>());
//! ```
//!
//! ## References
//!
//! * Capitaine, N., Wallace, P. T., Chapront, J. (2003), A&A 412, 567
//! * IERS Conventions (2010), §5.5.7
//! * SOFA routines `iauGmst06`, `iauGst06a`

use crate::astro::era::earth_rotation_angle;
use crate::qtty::*;
use crate::time::JulianDate;

/// Mean sidereal day length ≈ 0.9972696 solar days (23 h 56 m 4.09 s).
pub use crate::qtty::time::SIDEREAL_DAY;

// ════════════════════════════════════════════════════════════════════════
// IAU 2006 ERA-based sidereal time
// ════════════════════════════════════════════════════════════════════════

/// Greenwich Mean Sidereal Time, IAU 2006 (ERA-based).
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
    let t = (jd_tt.raw().value() - 2_451_545.0) / 36525.0;

    // Polynomial: accumulated precession of equinox in arcseconds
    // Coefficients from Capitaine et al. (2003), eq. 42
    let poly_as = 0.014_506 + 4_612.156_534 * t + 1.391_581_7 * t.powi(2)
        - 0.000_000_44 * t.powi(3)
        - 0.000_029_956 * t.powi(4)
        - 0.000_000_036_8 * t.powi(5);

    let gmst = era + Arcseconds::new(poly_as).to::<Radian>();
    gmst.wrap_pos()
}

/// Greenwich Apparent Sidereal Time, IAU 2006/2000A.
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
    (gmst + ee).wrap_pos()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::TAU;

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
            let jd = JulianDate::from_raw_unchecked(qtty::Day::new(jd));
            let gmst = gmst_iau2006(jd, jd);
            assert!(gmst >= Radians::new(0.0), "GMST should be ≥ 0");
            assert!(gmst < Radians::new(TAU), "GMST should be < 2π");
        }
    }

    #[test]
    fn gast_iau2006_close_to_gmst() {
        // GAST = GMST + Δψ·cos(ε). For small nutation (~17″ max), |GAST−GMST| < 20″.
        let jd = JulianDate::from_raw_unchecked(qtty::Day::new(2_460_000.5));
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
