// SPDX-License-Identifier: AGPL-3.0-only
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
//! the IAU 2000 equation of the equinoxes, including the complementary
//! periodic terms required by modern UT1 conventions.
//!
//! ## Example
//!
//! ```rust
//! use chrono::prelude::*;
//! use siderust::astro::sidereal::gmst_iau2006;
//! use siderust::qtty::*;
//! use siderust::time::{JulianDate, JD, TT, Time, UTC};
//!
//! let jd: JulianDate = Time::<UTC>::from_chrono(Utc::now()).to::<TT>().to::<JD>();
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
    let t = (jd_tt.raw().value() - 2_451_545.0_f64) / 36_525.0_f64;

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
/// This adds the IAU 2000 equation of the equinoxes to mean sidereal time:
/// the dominant nutation term and the complementary periodic terms required
/// by modern UT1 conventions.
///
/// `jd_ut1`: Julian Date on the UT1 time scale.
/// `jd_tt`:  Julian Date on the TT time scale.
/// `dpsi`:   nutation in longitude (Δψ) in radians.
/// `mean_obliquity`: ε_A in radians.
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
    mean_obliquity: Radians,
) -> Radians {
    let gmst = gmst_iau2006(jd_ut1, jd_tt);
    let ee = crate::astro::era::equation_of_the_equinoxes_iau2000(jd_tt, dpsi, mean_obliquity);
    (gmst + ee).wrap_pos()
}

/// Greenwich Apparent Sidereal Time, IAU 2006/2000A (high-level).
///
/// Computes nutation and mean obliquity internally from `jd_tt` using
/// [`nutation_iau2000b`](crate::astro::nutation::nutation_iau2000b), then
/// evaluates [`gast_iau2006`]. Prefer this at call sites that do not already
/// hold nutation angles; use [`gast_iau2006`] when supplying precomputed
/// nutation for batch transforms.
#[inline]
pub fn gast_iau2006a(jd_ut1: JulianDate, jd_tt: JulianDate) -> Radians {
    let nut = crate::astro::nutation::nutation_iau2000b(jd_tt);
    gast_iau2006(jd_ut1, jd_tt, nut.dpsi, nut.mean_obliquity)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::TAU;

    #[test]
    fn gmst_iau2006_at_j2000() {
        // At J2000.0, GMST ≈ ERA ≈ 280.46° (polynomial term is tiny at t=0)
        let jd = crate::J2000;
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
            let jd = crate::time::JulianDate::new(jd);
            let gmst = gmst_iau2006(jd, jd);
            assert!(gmst >= Radians::new(0.0), "GMST should be ≥ 0");
            assert!(gmst < Radians::new(TAU), "GMST should be < 2π");
        }
    }

    #[test]
    fn gast_iau2006_close_to_gmst() {
        // GAST differs from GMST by the equation of the equinoxes.
        let jd = crate::time::JulianDate::new(2_460_000.5);
        let nutation = crate::astro::nutation::nutation_iau2000b(jd);
        let gast = gast_iau2006(jd, jd, nutation.dpsi, nutation.mean_obliquity);
        let gmst = gmst_iau2006(jd, jd);

        let diff_as = (gast - gmst).value().abs() * 206_264.806;
        assert!(
            diff_as < 20.0,
            "GAST–GMST = {:.3}″, max nutation effect should be < 20″",
            diff_as
        );
    }

    #[test]
    fn gast_iau2006_matches_sofa_gst06a_reference() {
        use crate::astro::nutation::{Iau2006A, NutationModel};

        let jd = crate::time::JulianDate::new(2_453_736.5);
        let nut = <Iau2006A as NutationModel>::nutation(jd);
        let gast = gast_iau2006(jd, jd, nut.dpsi, nut.mean_obliquity);

        // SOFA/ERFA `gst06a(2400000.5, 53736.0, 2400000.5, 53736.0)`.
        let sofa = 1.754_166_137_675_019_2_f64;
        assert!(
            (gast.value() - sofa).abs() < 1e-12,
            "GAST = {:.16e}, SOFA gst06a = {:.16e}",
            gast.value(),
            sofa
        );
    }

    #[test]
    fn gast_iau2006a_matches_low_level_gast() {
        let jd = crate::time::JulianDate::new(2_460_000.5);
        let nut = crate::astro::nutation::nutation_iau2000b(jd);
        let low = gast_iau2006(jd, jd, nut.dpsi, nut.mean_obliquity);
        let high = gast_iau2006a(jd, jd);
        assert!((low.value() - high.value()).abs() < 1.0e-15);
    }
}
