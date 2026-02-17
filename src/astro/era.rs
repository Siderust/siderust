// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Earth Rotation Angle (ERA) — IAU 2000
//!
//! The **Earth Rotation Angle** (ERA, θ) is the angle between the Celestial
//! Intermediate Origin (CIO) and the Terrestrial Intermediate Origin (TIO),
//! measured along the equator of the Celestial Intermediate Pole (CIP).
//!
//! ERA replaces Greenwich Apparent Sidereal Time (GAST) in the IAU 2000/2006
//! framework as the primary measure of Earth's rotation.
//!
//! ## Defining relation
//!
//! ```text
//! ERA = 2π × (0.7790572732640 + 1.00273781191135448 × Du)
//! ```
//!
//! where **Du** is the Julian UT1 date − 2451545.0 (i.e., days since
//! J2000.0 on the UT1 axis).
//!
//! This is a simple linear function of UT1, with no polynomial terms.
//! The angular velocity is exactly the ratio of a sidereal day to a solar day.
//!
//! ## Equation of the Origins
//!
//! The **equation of the origins** (EO) connects ERA to GAST:
//!
//! ```text
//! GAST = ERA + EO
//! ```
//!
//! where EO is a function of the precession-nutation angles.
//!
//! ## References
//!
//! * IAU 2000 Resolution B1.8
//! * IERS Conventions (2010), §5.4.4
//! * SOFA/ERFA routine `iauEra00`

use crate::time::JulianDate;
use qtty::*;
use std::f64::consts::TAU;

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
#[inline]
pub fn earth_rotation_angle(jd_ut1: JulianDate) -> Radians {
    // Du = JD(UT1) − J2000.0
    let du = (jd_ut1 - JulianDate::J2000).value();

    // ERA = 2π × (0.7790572732640 + 1.00273781191135448 × Du)
    // To maintain precision, split the fractional and integer parts:
    let frac = du.fract();
    let era = TAU * (0.779_057_273_264_0 + 0.002_737_811_911_354_48 * du + frac);

    Radians::new(era.rem_euclid(TAU))
}

/// Equation of the origins: connects ERA to Greenwich Apparent Sidereal Time.
///
/// ```text
/// EO = GAST − ERA
/// ```
///
/// This is computed from the CIO locator `s` and the precession-nutation
/// angles. In the equinox-based framework:
///
/// ```text
/// EO ≈ −(ψ̄ + Δψ) × cos(ε_A + Δε) − s
/// ```
///
/// where the CIO locator `s` is small (~miliarcseconds).
///
/// For most practical purposes, the dominant term is −Δψ·cos(ε).
///
/// ## References
/// * IERS Conventions (2010), §5.5.3
/// * SOFA routine `iauEe06a`
#[inline]
pub fn equation_of_the_origins(
    psib_plus_dpsi: Radians,
    epsa_plus_deps: Radians,
    cio_locator_s: Radians,
) -> Radians {
    // EO = −(ψ̄ + Δψ) × cos(ε_A + Δε) − s
    Radians::new(-(psib_plus_dpsi.value()) * epsa_plus_deps.cos() - cio_locator_s.value())
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
    Radians::new(dpsi.value() * true_obliquity.cos())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn era_at_j2000() {
        // At J2000.0 (JD 2451545.0), ERA ≈ 280.46° (GST at that epoch)
        let era = earth_rotation_angle(JulianDate::J2000);
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
        let era1 = earth_rotation_angle(JulianDate::new(2_451_545.0));
        let era2 = earth_rotation_angle(JulianDate::new(2_451_545.5));
        // Half a solar day ≈ half a sidereal rotation ≈ 180°
        // Actually ERA increases by 360.985...° per solar day, so ≈ 180.5° per half day
        let diff = ((era2 - era1).value()).rem_euclid(TAU);
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
            let era = earth_rotation_angle(JulianDate::new(jd));
            assert!(era >= Radians::new(0.0), "ERA should be ≥ 0");
            assert!(era < Radians::new(TAU), "ERA should be < 2π");
        }
    }
}
