// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # IAU 2006 Precession Model (Capitaine et al. 2003, P03)
//!
//! This module implements the **IAU 2006 precession** of the Earth's mean equator
//! and equinox using the **Fukushima-Williams** four-angle parameterisation.
//!
//! ## Method
//!
//! The precession matrix is constructed from four angles (all in arcseconds,
//! expressed as polynomials in Julian centuries `t` from J2000.0 TT):
//!
//! * **γ̄** (gamb): ecliptic pole x, GCRS → mean ecliptic
//! * **φ̄** (phib): ecliptic pole y, including frame bias
//! * **ψ̄** (psib): equator pole, including frame bias
//! * **ε_A** (epsa): obliquity of the ecliptic
//!
//! The precession matrix is:
//!
//! ```text
//! P = R₁(−ε_A) · R₃(−ψ̄) · R₁(φ̄) · R₃(γ̄)
//! ```
//!
//! This naturally incorporates the frame bias between GCRS and the J2000 mean
//! equator/equinox.
//!
//! ## References
//!
//! * IAU 2006 Resolution B1 (precession)
//! * Hilton et al. (2006), Celestial Mechanics 94, 351–367
//! * Capitaine et al. (2003), A&A 412, 567–586
//! * SOFA/ERFA: `iauPfw06`, `iauObl06`, `iauPmat06`
//!
//! ## Accuracy
//!
//! Better than **0.1 mas** for epochs within ±200 years of J2000.0.
//! The polynomial series has no formal time limit but accuracy degrades
//! for epochs far from J2000.

use crate::time::JulianDate;
use affn::Rotation3;
use qtty::*;

// ═══════════════════════════════════════════════════════════════════════════
// Fukushima-Williams precession angles (SOFA iauPfw06)
// ═══════════════════════════════════════════════════════════════════════════

/// Fukushima-Williams precession angles for a given Julian Date (TT).
///
/// Returns `(gamb, phib, psib, epsa)` in **radians**, where:
/// * `gamb` (γ̄): F-W angle from GCRS to ecliptic pole (x)
/// * `phib` (φ̄): F-W angle from GCRS to ecliptic pole (y), with frame bias
/// * `psib` (ψ̄): F-W angle for equator pole, with frame bias
/// * `epsa` (ε_A): obliquity of the ecliptic (mean obliquity of date)
///
/// Polynomials from Hilton et al. (2006), Table 1.
/// Coefficients match SOFA `iauPfw06` / ERFA `eraPfw06`.
///
/// ## References
/// * IAU 2006 Resolution B1
/// * Hilton et al. (2006), Celestial Mechanics 94, 351–367
/// * SOFA routine `iauPfw06`
#[inline]
pub fn precession_fw_angles(jd: JulianDate) -> (Radians, Radians, Radians, Radians) {
    let t = jd.julian_centuries().value();
    let t2 = t * t;
    let t3 = t2 * t;
    let t4 = t3 * t;
    let t5 = t4 * t;

    // γ̄ (gamb) — arcseconds
    let gamb_as = -0.052_928
        + 10.556_378 * t
        + 0.493_2044 * t2
        - 0.000_312_38 * t3
        - 2.788e-6 * t4
        + 2.60e-8 * t5;

    // φ̄ (phib) — arcseconds
    let phib_as = 84381.412_819
        - 46.811_016 * t
        + 0.051_1268 * t2
        + 0.000_532_89 * t3
        - 4.40e-7 * t4
        - 1.76e-8 * t5;

    // ψ̄ (psib) — arcseconds
    let psib_as = -0.041_775
        + 5038.481_484 * t
        + 1.558_4175 * t2
        - 0.000_185_22 * t3
        - 2.6452e-5 * t4
        - 1.48e-8 * t5;

    // ε_A (epsa) — mean obliquity of date, arcseconds
    // IAU 2006 obliquity (Hilton et al. 2006 / SOFA iauObl06)
    let epsa_as = 84381.406
        - 46.836_769 * t
        - 0.000_1831 * t2
        + 0.002_003_40 * t3
        - 5.76e-7 * t4
        - 4.34e-8 * t5;

    let as_to_rad = |a: f64| Radians::new(a.to_radians() / 3600.0);
    (
        as_to_rad(gamb_as),
        as_to_rad(phib_as),
        as_to_rad(psib_as),
        as_to_rad(epsa_as),
    )
}

/// Mean obliquity of the ecliptic (IAU 2006).
///
/// Returns ε_A in **radians** for the given Julian Date (TT).
///
/// ```text
/// ε_A = 84381.406″ − 46.836769″·t − 0.0001831″·t² + 0.00200340″·t³
///       − 0.000000576″·t⁴ − 0.0000000434″·t⁵
/// ```
///
/// ## References
/// * IAU 2006 Resolution B1
/// * SOFA routine `iauObl06`
#[inline]
pub fn mean_obliquity_iau2006(jd: JulianDate) -> Radians {
    let t = jd.julian_centuries().value();
    let t2 = t * t;
    let t3 = t2 * t;
    let t4 = t3 * t;
    let t5 = t4 * t;

    let epsa_as = 84381.406
        - 46.836_769 * t
        - 0.000_1831 * t2
        + 0.002_003_40 * t3
        - 5.76e-7 * t4
        - 4.34e-8 * t5;

    Radians::new(epsa_as.to_radians() / 3600.0)
}

/// Mean obliquity of the ecliptic at J2000.0 (IAU 2006): 84381.406″.
///
/// ## References
/// * IAU 2006 Resolution B1
pub const J2000_MEAN_OBLIQUITY_ARCSEC: f64 = 84381.406;

// ═══════════════════════════════════════════════════════════════════════════
// Rotation matrix construction
// ═══════════════════════════════════════════════════════════════════════════

#[inline]
fn rotation_x(angle: Radians) -> Rotation3 {
    let (s, c) = angle.sin_cos();
    Rotation3::from_matrix([[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]])
}

#[inline]
fn rotation_z(angle: Radians) -> Rotation3 {
    let (s, c) = angle.sin_cos();
    Rotation3::from_matrix([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
}

/// Construct the Fukushima-Williams precession matrix from four angles.
///
/// ```text
/// P = R₁(−ε_A) · R₃(−ψ̄) · R₁(φ̄) · R₃(γ̄)
/// ```
///
/// This matrix transforms vectors from the GCRS (≈ ICRS) to the mean
/// equator and equinox of date. Frame bias is included.
///
/// ## References
/// * SOFA routine `iauFw2m`
#[inline]
pub fn fw_matrix(gamb: Radians, phib: Radians, psib: Radians, epsa: Radians) -> Rotation3 {
    rotation_x(-epsa) * rotation_z(-psib) * rotation_x(phib) * rotation_z(gamb)
}

/// IAU 2006 precession matrix from GCRS to mean equator/equinox of `jd`.
///
/// This matrix includes the frame bias, so it transforms directly from
/// GCRS (≈ ICRS) to the mean equatorial frame of date.
///
/// The `JulianDate` input is interpreted as TT.
///
/// ## References
/// * IAU 2006 Resolution B1
/// * SOFA routine `iauPmat06`
pub fn precession_matrix_iau2006(jd: JulianDate) -> Rotation3 {
    let (gamb, phib, psib, epsa) = precession_fw_angles(jd);
    fw_matrix(gamb, phib, psib, epsa)
}

/// Precession-nutation matrix (NPB) from GCRS to true equator/equinox of date.
///
/// Given nutation corrections Δψ and Δε (in radians), this constructs:
///
/// ```text
/// NPB = R₁(−(ε_A + Δε)) · R₃(−(ψ̄ + Δψ)) · R₁(φ̄) · R₃(γ̄)
/// ```
///
/// This is the complete frame rotation from GCRS to the true (apparent)
/// equatorial frame, including frame bias, precession, and nutation.
///
/// ## References
/// * IERS Conventions (2010), §5.6
/// * SOFA routine `iauPnm06a`
pub fn precession_nutation_matrix(jd: JulianDate, dpsi: Radians, deps: Radians) -> Rotation3 {
    let (gamb, phib, psib, epsa) = precession_fw_angles(jd);
    fw_matrix(gamb, phib, psib + dpsi, epsa + deps)
}

#[cfg(test)]
mod tests {
    use super::*;
    use qtty::{Degrees, Radians};

    #[test]
    fn mean_obliquity_at_j2000() {
        let eps = mean_obliquity_iau2006(JulianDate::J2000);
        // IAU 2006: 84381.406″ = 23.4392911111...°
        let expected_deg = 84381.406 / 3600.0;
        assert!(
            (eps.to::<Degree>().value() - expected_deg).abs() < 1e-10,
            "obliquity at J2000 = {}°, expected {}°",
            eps.to::<Degree>(),
            expected_deg
        );
    }

    #[test]
    fn fw_angles_at_j2000_are_approximately_identity() {
        // At J2000, t=0: gamb≈-0.053″, phib≈84381.413″, psib≈-0.042″, epsa=84381.406″
        // The precession matrix at J2000 should be close to identity (with frame bias).
        let mat = precession_matrix_iau2006(JulianDate::J2000);
        let m = mat.as_matrix();
        // Diagonal should be very close to 1
        for i in 0..3 {
            assert!(
                (m[i][i] - 1.0).abs() < 1e-7,
                "diagonal[{}] = {}, expected ~1.0",
                i,
                m[i][i]
            );
        }
    }

    #[test]
    fn precession_matrix_j2025_reasonable() {
        // JD of approximately 2025-01-01
        let jd = JulianDate::new(2_460_676.5);
        let mat = precession_matrix_iau2006(jd);
        let m = mat.as_matrix();

        // Precession over 25 years is small: off-diagonal < 0.01
        for i in 0..3 {
            for j in 0..3 {
                if i == j {
                    assert!(
                        (m[i][j] - 1.0).abs() < 0.001,
                        "diagonal[{}][{}] too far from 1: {}",
                        i,
                        j,
                        m[i][j]
                    );
                } else {
                    assert!(
                        m[i][j].abs() < 0.01,
                        "off-diagonal[{}][{}] too large: {}",
                        i,
                        j,
                        m[i][j]
                    );
                }
            }
        }
    }

    #[test]
    fn precession_nutation_matrix_includes_corrections() {
        let jd = JulianDate::new(2_460_000.5);
        let mat_prec = precession_matrix_iau2006(jd);
        let mat_pn =
            precession_nutation_matrix(jd, Radians::new(1e-5), Radians::new(1e-5));

        // The matrices should differ slightly due to nutation
        let m1 = mat_prec.as_matrix();
        let m2 = mat_pn.as_matrix();
        let mut max_diff = 0.0f64;
        for i in 0..3 {
            for j in 0..3 {
                max_diff = max_diff.max((m1[i][j] - m2[i][j]).abs());
            }
        }
        assert!(
            max_diff > 1e-8,
            "nutation should cause detectable difference"
        );
        assert!(
            max_diff < 1e-3,
            "nutation difference too large: {}",
            max_diff
        );
    }

    #[test]
    fn mean_obliquity_decreases_with_time() {
        let eps_2000 = mean_obliquity_iau2006(JulianDate::J2000);
        let eps_2100 = mean_obliquity_iau2006(JulianDate::new(2_488_069.5));
        // Obliquity is currently decreasing at ~47″/century
        assert!(eps_2100 < eps_2000, "obliquity should decrease over time");
        let diff_arcsec = (eps_2000 - eps_2100).to::<Degree>().value() * 3600.0;
        assert!(
            (diff_arcsec - 47.0_f64).abs() < 2.0,
            "obliquity change over 1 century ≈ 47″, got {}″",
            diff_arcsec
        );
    }
}
