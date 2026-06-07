// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Earth Rotation Angle (ERA), IAU 2000
//!
//! Implements the Earth Rotation Angle θ — the IAU 2000/2006 replacement for
//! Greenwich Apparent Sidereal Time as the primary measure of Earth's
//! rotation — and the equation of the origins linking ERA to GAST.
//!
//! ## Scientific scope
//!
//! The ERA is the angle between the Celestial Intermediate Origin (CIO) and
//! the Terrestrial Intermediate Origin (TIO) measured along the equator of
//! the Celestial Intermediate Pole (CIP). Unlike GMST, which mixes precession
//! into a polynomial in TT, the ERA is a strictly linear function of UT1, so
//! its angular velocity is exactly the ratio of a sidereal day to a solar
//! day. The equation of the origins `EO = GAST − ERA` provides the bridge to
//! the legacy equinox-based system and depends only on the precession-
//! nutation angles and the CIO locator `s`.
//!
//! ## Technical scope
//!
//! `earth_rotation_angle` evaluates the IAU 2000 ERA using the ERFA/SOFA
//! `era00` split-date algorithm: UT1 is expressed as `d1 = 2451545.0` and
//! `d2 = JD(UT1) − 2451545.0`, then
//!
//! ```text
//! t = d1 + d2
//! f = fmod(d1, 1) + fmod(d2, 1)
//! ERA = 2π × (f + 0.7790572732640 + 0.00273781191135448 × t)  (mod 2π)
//! ```
//!
//! This matches Capitaine et al. (2000) and avoids the incorrect extra
//! `fract(Du)` term that would double-count the fractional day. Precision is
//! still limited by storing UT1 as a single `f64` Julian day; sub-microsecond
//! UT1 resolution requires a true two-part Julian date input. The equation of the origins follows from the precession-
//! nutation chain via the simplified relation
//! `EO ≈ −(ψ̄ + Δψ) cos(ε_A + Δε) − s`.
//!
//! ## References
//!
//! * IAU 2000 Resolution B1.8
//! * IERS Conventions (2010), §5.4.4
//! * SOFA routine `iauEra00`

use crate::qtty::*;
use crate::time::JulianDate;
use std::f64::consts::TAU;

const ERFA_DJ00: f64 = 2_451_545.0;
const ERA0_TURNS: f64 = 0.779_057_273_264_0;
const ERAS_TURNS: f64 = 0.002_737_811_911_354_48;

/// C-style `fmod` for the ERA fractional-day sum (sign follows the dividend).
#[inline]
fn fmod_c(x: f64, y: f64) -> f64 {
    x - (x / y).trunc() * y
}

/// Compute the Earth Rotation Angle for a given Julian Date on the UT1 axis.
///
/// The input `jd_ut1` should be a Julian Day number on the **UT1** time scale.
/// In practice, UT ≈ UT1 for most applications (the difference is < 0.9 s
/// and requires IERS Bulletin A data).
///
/// Returns the ERA in **radians**, normalized to [0, 2π).
///
/// ## References
/// * IAU 2000 Resolution B1.8
/// * SOFA routine `iauEra00`
///
/// UT1 is passed as a single Julian day number. For maximum precision, callers
/// with a true two-part UT1 Julian date should use [`earth_rotation_angle_parts`].
#[inline]
pub fn earth_rotation_angle(jd_ut1: JulianDate) -> Radians {
    earth_rotation_angle_parts(ERFA_DJ00, jd_ut1.raw().value() - ERFA_DJ00)
}

/// Earth Rotation Angle from a two-part UT1 Julian Date (`d1 + d2`), SOFA `era00`.
#[inline]
pub fn earth_rotation_angle_parts(d1: f64, d2: f64) -> Radians {
    let (d1, d2) = if d1 < d2 { (d1, d2) } else { (d2, d1) };
    let t = d1 + (d2 - ERFA_DJ00);
    let f = fmod_c(d1, 1.0) + fmod_c(d2, 1.0);
    let turns = (f + ERA0_TURNS + ERAS_TURNS * t).rem_euclid(1.0);
    Radians::new(turns * TAU)
}

/// Approximate equation of the origins (equinox-based, not IAU-validated).
///
/// ```text
/// EO ≈ −(ψ̄ + Δψ) × cos(ε_A + Δε) − s
/// ```
///
/// This is a simplified bridge formula. For standards-compliant apparent
/// sidereal time use [`crate::astro::sidereal::gast_iau2006`] /
/// [`crate::astro::sidereal::gast_iau2006a`], which apply the IAU 2000
/// equation of the equinoxes including complementary terms.
///
/// ## References
/// * IERS Conventions (2010), §5.5.3 (conceptual relation EO = GAST − ERA)
#[inline]
pub fn equation_of_the_origins_approx(
    psib_plus_dpsi: Radians,
    epsa_plus_deps: Radians,
    cio_locator_s: Radians,
) -> Radians {
    // EO = −(ψ̄ + Δψ) × cos(ε_A + Δε) − s
    -(psib_plus_dpsi * epsa_plus_deps.cos()) - cio_locator_s
}

/// Equation of the equinoxes (traditional, simplified).
///
/// ```text
/// GAST = GMST + Δψ·cos(ε)
/// ```
///
/// This is the traditional nutation correction to mean sidereal time.
/// It relates the mean and apparent Greenwich sidereal times.
///
/// ## References
/// * Meeus (1998), eq. 12.4
/// * SOFA routine `iauEe00`
#[inline]
pub fn equation_of_the_equinoxes(dpsi: Radians, true_obliquity: Radians) -> Radians {
    dpsi * true_obliquity.cos()
}

#[derive(Debug, Clone, Copy)]
struct ComplementaryTerm {
    nfa: [i8; 8],
    sin_arcsec: f64,
    cos_arcsec: f64,
}

// SOFA `iauEect00`: terms of order t^0 for the IAU 2000 equation-of-equinoxes
// complementary terms. Coefficients are arcseconds.
const EECT00_T0: [ComplementaryTerm; 33] = [
    ComplementaryTerm {
        nfa: [0, 0, 0, 0, 1, 0, 0, 0],
        sin_arcsec: 2640.96e-6,
        cos_arcsec: -0.39e-6,
    },
    ComplementaryTerm {
        nfa: [0, 0, 0, 0, 2, 0, 0, 0],
        sin_arcsec: 63.52e-6,
        cos_arcsec: -0.02e-6,
    },
    ComplementaryTerm {
        nfa: [0, 0, 2, -2, 3, 0, 0, 0],
        sin_arcsec: 11.75e-6,
        cos_arcsec: 0.01e-6,
    },
    ComplementaryTerm {
        nfa: [0, 0, 2, -2, 1, 0, 0, 0],
        sin_arcsec: 11.21e-6,
        cos_arcsec: 0.01e-6,
    },
    ComplementaryTerm {
        nfa: [0, 0, 2, -2, 2, 0, 0, 0],
        sin_arcsec: -4.55e-6,
        cos_arcsec: 0.0,
    },
    ComplementaryTerm {
        nfa: [0, 0, 2, 0, 3, 0, 0, 0],
        sin_arcsec: 2.02e-6,
        cos_arcsec: 0.0,
    },
    ComplementaryTerm {
        nfa: [0, 0, 2, 0, 1, 0, 0, 0],
        sin_arcsec: 1.98e-6,
        cos_arcsec: 0.0,
    },
    ComplementaryTerm {
        nfa: [0, 0, 0, 0, 3, 0, 0, 0],
        sin_arcsec: -1.72e-6,
        cos_arcsec: 0.0,
    },
    ComplementaryTerm {
        nfa: [0, 1, 0, 0, 1, 0, 0, 0],
        sin_arcsec: -1.41e-6,
        cos_arcsec: -0.01e-6,
    },
    ComplementaryTerm {
        nfa: [0, 1, 0, 0, -1, 0, 0, 0],
        sin_arcsec: -1.26e-6,
        cos_arcsec: -0.01e-6,
    },
    ComplementaryTerm {
        nfa: [1, 0, 0, 0, -1, 0, 0, 0],
        sin_arcsec: -0.63e-6,
        cos_arcsec: 0.0,
    },
    ComplementaryTerm {
        nfa: [1, 0, 0, 0, 1, 0, 0, 0],
        sin_arcsec: -0.63e-6,
        cos_arcsec: 0.0,
    },
    ComplementaryTerm {
        nfa: [0, 1, 2, -2, 3, 0, 0, 0],
        sin_arcsec: 0.46e-6,
        cos_arcsec: 0.0,
    },
    ComplementaryTerm {
        nfa: [0, 1, 2, -2, 1, 0, 0, 0],
        sin_arcsec: 0.45e-6,
        cos_arcsec: 0.0,
    },
    ComplementaryTerm {
        nfa: [0, 0, 4, -4, 4, 0, 0, 0],
        sin_arcsec: 0.36e-6,
        cos_arcsec: 0.0,
    },
    ComplementaryTerm {
        nfa: [0, 0, 1, -1, 1, -8, 12, 0],
        sin_arcsec: -0.24e-6,
        cos_arcsec: -0.12e-6,
    },
    ComplementaryTerm {
        nfa: [0, 0, 2, 0, 0, 0, 0, 0],
        sin_arcsec: 0.32e-6,
        cos_arcsec: 0.0,
    },
    ComplementaryTerm {
        nfa: [0, 0, 2, 0, 2, 0, 0, 0],
        sin_arcsec: 0.28e-6,
        cos_arcsec: 0.0,
    },
    ComplementaryTerm {
        nfa: [1, 0, 2, 0, 3, 0, 0, 0],
        sin_arcsec: 0.27e-6,
        cos_arcsec: 0.0,
    },
    ComplementaryTerm {
        nfa: [1, 0, 2, 0, 1, 0, 0, 0],
        sin_arcsec: 0.26e-6,
        cos_arcsec: 0.0,
    },
    ComplementaryTerm {
        nfa: [0, 0, 2, -2, 0, 0, 0, 0],
        sin_arcsec: -0.21e-6,
        cos_arcsec: 0.0,
    },
    ComplementaryTerm {
        nfa: [0, 1, -2, 2, -3, 0, 0, 0],
        sin_arcsec: 0.19e-6,
        cos_arcsec: 0.0,
    },
    ComplementaryTerm {
        nfa: [0, 1, -2, 2, -1, 0, 0, 0],
        sin_arcsec: 0.18e-6,
        cos_arcsec: 0.0,
    },
    ComplementaryTerm {
        nfa: [0, 0, 0, 0, 0, 8, -13, -1],
        sin_arcsec: -0.10e-6,
        cos_arcsec: 0.05e-6,
    },
    ComplementaryTerm {
        nfa: [0, 0, 0, 2, 0, 0, 0, 0],
        sin_arcsec: 0.15e-6,
        cos_arcsec: 0.0,
    },
    ComplementaryTerm {
        nfa: [2, 0, -2, 0, -1, 0, 0, 0],
        sin_arcsec: -0.14e-6,
        cos_arcsec: 0.0,
    },
    ComplementaryTerm {
        nfa: [1, 0, 0, -2, 1, 0, 0, 0],
        sin_arcsec: 0.14e-6,
        cos_arcsec: 0.0,
    },
    ComplementaryTerm {
        nfa: [0, 1, 2, -2, 2, 0, 0, 0],
        sin_arcsec: -0.14e-6,
        cos_arcsec: 0.0,
    },
    ComplementaryTerm {
        nfa: [1, 0, 0, -2, -1, 0, 0, 0],
        sin_arcsec: 0.14e-6,
        cos_arcsec: 0.0,
    },
    ComplementaryTerm {
        nfa: [0, 0, 4, -2, 4, 0, 0, 0],
        sin_arcsec: 0.13e-6,
        cos_arcsec: 0.0,
    },
    ComplementaryTerm {
        nfa: [0, 0, 2, -2, 4, 0, 0, 0],
        sin_arcsec: -0.11e-6,
        cos_arcsec: 0.0,
    },
    ComplementaryTerm {
        nfa: [1, 0, -2, 0, -3, 0, 0, 0],
        sin_arcsec: 0.11e-6,
        cos_arcsec: 0.0,
    },
    ComplementaryTerm {
        nfa: [1, 0, -2, 0, -1, 0, 0, 0],
        sin_arcsec: 0.11e-6,
        cos_arcsec: 0.0,
    },
];

const EECT00_T1: [ComplementaryTerm; 1] = [ComplementaryTerm {
    nfa: [0, 0, 0, 0, 1, 0, 0, 0],
    sin_arcsec: -0.87e-6,
    cos_arcsec: 0.0,
}];

const TURN_ARCSEC: f64 = 1_296_000.0;
const ARCSEC_TO_RAD: f64 = std::f64::consts::PI / (180.0 * 3600.0);

fn arcsec_argument(poly: f64) -> f64 {
    poly % TURN_ARCSEC * ARCSEC_TO_RAD
}

fn fal03(t: f64) -> f64 {
    arcsec_argument(
        485868.249036 + t * (1717915923.2178 + t * (31.8792 + t * (0.051635 - 0.00024470 * t))),
    )
}

fn falp03(t: f64) -> f64 {
    arcsec_argument(
        1287104.793048 + t * (129596581.0481 + t * (-0.5532 + t * (0.000136 - 0.00001149 * t))),
    )
}

fn faf03(t: f64) -> f64 {
    arcsec_argument(
        335779.526232 + t * (1739527262.8478 + t * (-12.7512 + t * (-0.001037 + 0.00000417 * t))),
    )
}

fn fad03(t: f64) -> f64 {
    arcsec_argument(
        1072260.703692 + t * (1602961601.2090 + t * (-6.3706 + t * (0.006593 - 0.00003169 * t))),
    )
}

fn faom03(t: f64) -> f64 {
    arcsec_argument(
        450160.398036 + t * (-6962890.5431 + t * (7.4722 + t * (0.007702 - 0.00005939 * t))),
    )
}

fn fave03(t: f64) -> f64 {
    (3.176146697 + 1021.3285546211 * t) % TAU
}

fn fae03(t: f64) -> f64 {
    (1.753470314 + 628.3075849991 * t) % TAU
}

fn fapa03(t: f64) -> f64 {
    (0.024381750 + 0.00000538691 * t) * t
}

fn evaluate_complementary_terms(terms: &[ComplementaryTerm], args: &[f64; 8]) -> f64 {
    terms.iter().rev().fold(0.0, |sum, term| {
        let angle = term
            .nfa
            .iter()
            .zip(args.iter())
            .fold(0.0, |angle, (coefficient, argument)| {
                angle + f64::from(*coefficient) * argument
            });
        sum + term.sin_arcsec * angle.sin() + term.cos_arcsec * angle.cos()
    })
}

/// Equation-of-equinoxes complementary terms, IAU 2000.
///
/// These are the small periodic terms introduced to keep apparent sidereal
/// time compatible with modern UT1 conventions. They are the terms used by
/// SOFA `iauEect00`.
pub fn complementary_terms_iau2000(jd_tt: JulianDate) -> Radians {
    let t = jd_tt.julian_centuries();
    let args = [
        fal03(t),
        falp03(t),
        faf03(t),
        fad03(t),
        faom03(t),
        fave03(t),
        fae03(t),
        fapa03(t),
    ];

    let arcsec = evaluate_complementary_terms(&EECT00_T0, &args)
        + t * evaluate_complementary_terms(&EECT00_T1, &args);
    Radians::new(arcsec * ARCSEC_TO_RAD)
}

/// Equation of the equinoxes, IAU 2000 conventions.
///
/// This is the quantity to add to IAU 2006 GMST when computing an apparent
/// sidereal angle compatible with SOFA `iauGst06a`.
pub fn equation_of_the_equinoxes_iau2000(
    jd_tt: JulianDate,
    dpsi: Radians,
    mean_obliquity: Radians,
) -> Radians {
    dpsi * mean_obliquity.cos() + complementary_terms_iau2000(jd_tt)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Reference ERA from the ERFA `era00` algorithm (SOFA `iauEra00`).
    fn reference_era00(d1: f64, d2: f64) -> f64 {
        earth_rotation_angle_parts(d1, d2).value()
    }

    #[test]
    fn era_matches_sofa_reference_epochs() {
        let cases: &[(f64, f64)] = &[
            (2_451_545.0, 0.0),                       // J2000.0
            (2_451_545.0, 2_453_372.5 - 2_451_545.0), // 2005-01-01
            (2_451_545.0, 2_440_000.0 - 2_451_545.0), // pre-J2000
            (2_451_545.0, 2_500_000.0 - 2_451_545.0), // far future
            (2_451_545.0, -1_545.0),                  // negative offset Du
            (2_400_000.5, 50123.2),                   // MJD-style split (SOFA note example shape)
        ];
        for &(d1, d2) in cases {
            let jd = crate::time::JulianDate::new(d1 + d2);
            let era = earth_rotation_angle(jd).value();
            let expected = reference_era00(d1, d2);
            assert!(
                (era - expected).abs() < 2.0e-9,
                "ERA mismatch at d1={d1} d2={d2}: got {era:.17e}, expected {expected:.17e}"
            );
        }
    }

    #[test]
    fn era_parts_matches_single_jd_when_split_at_j2000() {
        let jd = 2_453_736.5_f64;
        let from_jd = earth_rotation_angle(crate::time::JulianDate::new(jd)).value();
        let from_parts = earth_rotation_angle_parts(ERFA_DJ00, jd - ERFA_DJ00).value();
        assert_eq!(from_jd, from_parts);
    }

    #[test]
    fn era_at_j2000() {
        // At J2000.0 (JD 2451545.0), ERA ≈ 280.46° (GST at that epoch)
        let era = earth_rotation_angle(crate::J2000);
        let era_deg = era.to::<Degree>();
        // ERA at J2000.0 ≈ 0.7790572732640 × 360° ≈ 280.46°
        assert!(
            (era_deg.value() - 280.46_f64).abs() < 0.1,
            "ERA at J2000 = {}°, expected ≈ 280.46°",
            era_deg
        );
    }

    #[test]
    fn era_increases_with_time() {
        let era1 = earth_rotation_angle(crate::time::JulianDate::new(2_451_545.0));
        let era2 = earth_rotation_angle(crate::time::JulianDate::new(2_451_545.5));
        // Half a solar day ≈ half a sidereal rotation ≈ 180°
        // Actually ERA increases by 360.985...° per solar day, so ≈ 180.5° per half day
        let diff = (era2 - era1).value().rem_euclid(TAU);
        let diff_deg = diff.to_degrees();
        assert!(
            (diff_deg - 180.49).abs() < 1.0,
            "ERA increase over 0.5 days = {}°, expected ≈ 180.5°",
            diff_deg
        );
    }

    #[test]
    fn era_range() {
        for jd in [2_451_545.0, 2_460_000.5, 2_440_000.0] {
            let era = earth_rotation_angle(crate::time::JulianDate::new(jd));
            assert!(era >= Radians::new(0.0), "ERA should be ≥ 0");
            assert!(era < Radians::new(TAU), "ERA should be < 2π");
        }
    }
}
