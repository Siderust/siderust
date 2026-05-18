// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # IAU 2006 Precession
//!
//! Provides precession of the Earth's mean equator and equinox using the
//! IAU 2006 model (Capitaine et al. 2003, P03) in the Fukushima-Williams
//! four-angle parameterisation, plus the combined precession-nutation matrix.
//!
//! ## Scientific scope
//!
//! The Earth's rotation axis is not fixed in inertial space: torques from
//! the Sun and Moon on the equatorial bulge produce a slow conical motion
//! with a period of ≈ 25,770 years (the *Platonic year*), known as
//! lunisolar precession; planetary tides add a much smaller contribution.
//! In equatorial coordinates this drifts the celestial poles and equinox
//! at ≈ 50″/yr. Whenever stellar positions from different epochs are
//! compared, coordinates must be precessed to a common reference date.
//! Precession is the **secular** part of the orientation; nutation is the
//! superposed periodic wobble.
//!
//! ## Technical scope
//!
//! The precession matrix is constructed from four Fukushima-Williams angles
//! `(γ̄, φ̄, ψ̄, ε_A)`, each evaluated as a 5th-order polynomial in Julian
//! centuries `t` from J2000.0 TT, and assembled as
//!
//! ```text
//! P = R₁(−ε_A) · R₃(−ψ̄) · R₁(φ̄) · R₃(γ̄),
//! ```
//!
//! which naturally incorporates the GCRS↔J2000 mean equator/equinox frame
//! bias. Public entry points include [`precession_fw_angles`], [`fw_matrix`],
//! [`precession_matrix_iau2006`], [`precession_nutation_matrix`] and
//! [`mean_obliquity_iau2006`]. Accuracy is better than 0.1 mas within
//! ±200 yr of J2000.0, degrading gracefully outside that window.
//!
//! ## References
//!
//! * IAU 2006 Resolution B1
//! * Capitaine, N., Wallace, P. T., Chapront, J. (2003), A&A 412, 567
//! * IERS Conventions (2010), §5.6
//! * SOFA routine `iauPfw06`

use crate::qtty::*;
use crate::time::JulianDate;
use affn::Rotation3;

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
/// Coefficients match SOFA `iauPfw06`.
///
/// ## References
/// * IAU 2006 Resolution B1
/// * Hilton et al. (2006), Celestial Mechanics 94, 351–367
/// * SOFA routine `iauPfw06`
#[inline]
pub fn precession_fw_angles(jd: JulianDate) -> (Radians, Radians, Radians, Radians) {
    let t = (jd.raw().value() - 2_451_545.0_f64) / 36_525.0_f64;
    let t2 = t * t;
    let t3 = t2 * t;
    let t4 = t3 * t;
    let t5 = t4 * t;

    // γ̄ (gamb), arcseconds
    let gamb_as =
        -0.052_928 + 10.556_378 * t + 0.493_204_4 * t2 - 0.000_312_38 * t3 - 2.788e-6 * t4
            + 2.60e-8 * t5;

    // φ̄ (phib), arcseconds
    let phib_as = 84_381.412_819 - 46.811_016 * t + 0.051_126_8 * t2 + 0.000_532_89 * t3
        - 4.40e-7 * t4
        - 1.76e-8 * t5;

    // ψ̄ (psib), arcseconds
    let psib_as = -0.041_775 + 5_038.481_484 * t + 1.558_417_5 * t2
        - 0.000_185_22 * t3
        - 2.6452e-5 * t4
        - 1.48e-8 * t5;

    // ε_A (epsa), mean obliquity of date, arcseconds
    // IAU 2006 obliquity (Hilton et al. 2006 / SOFA iauObl06)
    let epsa_as = 84381.406 - 46.836_769 * t - 0.000_183_1 * t2 + 0.002_003_40 * t3
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
    let t = (jd.raw().value() - 2_451_545.0_f64) / 36_525.0_f64;
    let t2 = t * t;
    let t3 = t2 * t;
    let t4 = t3 * t;
    let t5 = t4 * t;

    let epsa_as = 84381.406 - 46.836_769 * t - 0.000_183_1 * t2 + 0.002_003_40 * t3
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

/// Construct the Fukushima-Williams precession matrix from four angles.
///
/// The SOFA formula is:
/// ```text
/// P = R₁(−ε_A) · R₃(−ψ̄) · R₁(φ̄) · R₃(γ̄)
/// ```
/// where R₁, R₃ are the SOFA elementary rotation functions.
///
/// SOFA's Rx/Rz use the **opposite sign convention** from the standard
/// math rotation matrices: `Rx_SOFA(θ) = Rx_standard(−θ)`.
///
/// Translated to **standard** (siderust) convention:
/// ```text
/// P = Rx(ε_A) · Rz(ψ̄) · Rx(−φ̄) · Rz(−γ̄)
/// ```
///
/// This matrix transforms vectors from the GCRS (≈ ICRS) to the mean
/// equator and equinox of date. Frame bias is included.
///
/// ## References
/// * SOFA routine `iauFw2m`
#[inline]
pub fn fw_matrix(gamb: Radians, phib: Radians, psib: Radians, epsa: Radians) -> Rotation3 {
    // Fused 4-rotation constructor: ~35% faster than sequential composition
    Rotation3::fused_rx_rz_rx_rz(epsa, psib, -phib, -gamb)
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

/// Rotation matrix from GCRS (≈ ICRS) equatorial to ecliptic-of-date.
///
/// Combines the IAU 2006 precession matrix with a rotation about the X-axis
/// by the **of-date** mean obliquity:
///
/// ```text
/// R = Rx(−ε_A) · P_iau2006
/// ```
///
/// This converts a Cartesian vector in GCRS equatorial coordinates to the
/// **mean** ecliptic of date (no nutation).  Ecliptic longitude is measured
/// from the mean equinox.
///
/// ## Parameters
/// * `jd`: Julian Date in TT scale.
pub fn gcrs_to_ecliptic_of_date_matrix(jd: JulianDate) -> Rotation3 {
    let prec = precession_matrix_iau2006(jd);
    let eps_a = mean_obliquity_iau2006(jd);
    Rotation3::rx(-eps_a) * prec
}

/// Rotation matrix from GCRS to the **true** ecliptic of date.
///
/// Combines the IAU 2006/2000B precession-nutation (NPB) matrix with a
/// rotation by the **true** obliquity (ε_A + Δε):
///
/// ```text
/// R = Rx(−ε_true) · NPB
/// ```
///
/// Ecliptic longitude is measured from the **true** equinox (precession +
/// nutation in longitude applied).
///
/// ## Parameters
/// * `jd`: Julian Date in TT scale.
pub fn gcrs_to_true_ecliptic_of_date_matrix(jd: JulianDate) -> Rotation3 {
    let nut = crate::astro::nutation::nutation_iau2000b(jd);
    let npb = precession_nutation_matrix(jd, nut.dpsi, nut.deps);
    let eps_true = nut.mean_obliquity + nut.deps;
    Rotation3::rx(-eps_true) * npb
}

/// Rotation matrix from the **true** ecliptic of date to GCRS.
///
/// Transpose (inverse) of [`gcrs_to_true_ecliptic_of_date_matrix`].
pub fn true_ecliptic_of_date_to_gcrs_matrix(jd: JulianDate) -> Rotation3 {
    gcrs_to_true_ecliptic_of_date_matrix(jd).transpose()
}

/// Rotation matrix from ecliptic-of-date to GCRS (≈ ICRS) equatorial.
///
/// This is the transpose (inverse) of [`gcrs_to_ecliptic_of_date_matrix`].
///
/// ## Parameters
/// * `jd`: Julian Date in TT scale.
pub fn ecliptic_of_date_to_gcrs_matrix(jd: JulianDate) -> Rotation3 {
    gcrs_to_ecliptic_of_date_matrix(jd).transpose()
}

/// Rotation matrix from mean equatorial of date to mean ecliptic of date.
///
/// This transformation applies only the obliquity rotation (without nutation),
/// converting from the mean equatorial frame (precessed but not nutated) to the
/// mean ecliptic frame of the same epoch.
///
/// ## Parameters
/// * `jd`: Julian Date in TT scale.
pub fn mean_equatorial_to_ecliptic_of_date_matrix(jd: JulianDate) -> Rotation3 {
    let eps_a = mean_obliquity_iau2006(jd);
    Rotation3::rx(-eps_a)
}

/// Rotation matrix from mean equatorial of date to **true** ecliptic of date.
///
/// Applies the IAU 2000B nutation matrix and rotates by the true obliquity:
///
/// ```text
/// R = Rx(−ε_true) · N
/// ```
///
/// where N is the equinox-based nutation matrix and ε_true = ε_A + Δε.
///
/// ## Parameters
/// * `jd`: Julian Date in TT scale.
pub fn mean_equatorial_to_true_ecliptic_of_date_matrix(jd: JulianDate) -> Rotation3 {
    let nut = crate::astro::nutation::nutation_iau2000b(jd);
    let eps = nut.mean_obliquity;
    let dpsi = nut.dpsi;
    let deps = nut.deps;
    let nutation_matrix = affn::Rotation3::fused_rx_rz_rx(eps + deps, dpsi, -eps);
    let eps_true = eps + deps;
    Rotation3::rx(-eps_true) * nutation_matrix
}

/// Rotation matrix from true ecliptic of date to mean equatorial of date.
///
/// Transpose (inverse) of [`mean_equatorial_to_true_ecliptic_of_date_matrix`].
pub fn true_ecliptic_of_date_to_mean_equatorial_matrix(jd: JulianDate) -> Rotation3 {
    mean_equatorial_to_true_ecliptic_of_date_matrix(jd).transpose()
}

/// Rotation matrix from ecliptic of date to mean equatorial of date.
///
/// This is the transpose (inverse) of [`mean_equatorial_to_ecliptic_of_date_matrix`].
///
/// ## Parameters
/// * `jd`: Julian Date in TT scale.
pub fn ecliptic_of_date_to_mean_equatorial_matrix(jd: JulianDate) -> Rotation3 {
    mean_equatorial_to_ecliptic_of_date_matrix(jd).transpose()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::qtty::Radians;

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
        for (i, row) in m.iter().enumerate().take(3) {
            assert!(
                (row[i] - 1.0).abs() < 1e-7,
                "diagonal[{}] = {}, expected ~1.0",
                i,
                row[i]
            );
        }
    }

    #[test]
    fn precession_matrix_j2025_reasonable() {
        // JD of approximately 2025-01-01
        let jd = JulianDate::from_raw_unchecked(qtty::Day::new(2_460_676.5));
        let mat = precession_matrix_iau2006(jd);
        let m = mat.as_matrix();

        // Precession over 25 years is small: off-diagonal < 0.01
        for (i, row) in m.iter().enumerate().take(3) {
            for (j, &val) in row.iter().enumerate().take(3) {
                if i == j {
                    assert!(
                        (val - 1.0).abs() < 0.001,
                        "diagonal[{}][{}] too far from 1: {}",
                        i,
                        j,
                        val
                    );
                } else {
                    assert!(
                        val.abs() < 0.01,
                        "off-diagonal[{}][{}] too large: {}",
                        i,
                        j,
                        val
                    );
                }
            }
        }
    }

    #[test]
    fn precession_nutation_matrix_includes_corrections() {
        let jd = JulianDate::from_raw_unchecked(qtty::Day::new(2_460_000.5));
        let mat_prec = precession_matrix_iau2006(jd);
        let mat_pn = precession_nutation_matrix(jd, Radians::new(1e-5), Radians::new(1e-5));

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
        let eps_2100 = mean_obliquity_iau2006(JulianDate::from_raw_unchecked(qtty::Day::new(2_488_069.5)));
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
