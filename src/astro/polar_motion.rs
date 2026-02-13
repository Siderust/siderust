// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Polar Motion — IERS Conventions
//!
//! Polar motion describes the deviation of the Earth's rotation axis from
//! its crust-fixed reference position. It is characterized by the pole
//! coordinates **(xₚ, yₚ)** published by the IERS.
//!
//! ## The W matrix
//!
//! The polar motion matrix **W** transforms from the Terrestrial Intermediate
//! Reference System (TIRS) to the International Terrestrial Reference System
//! (ITRS):
//!
//! ```text
//! W = R₃(−s') · R₂(xₚ) · R₁(yₚ)
//! ```
//!
//! where **s'** is the TIO locator (Terrestrial Intermediate Origin locator),
//! a very small quantity (~μas) that accounts for the motion of the TIO on
//! the Earth's surface.
//!
//! ## TIO locator s'
//!
//! ```text
//! s' ≈ −47 μas × t  (where t is Julian centuries from J2000)
//! ```
//!
//! This is negligible for most applications (< 0.01 mas over a century).
//!
//! ## References
//!
//! * IERS Conventions (2010), §5.4.2
//! * SOFA routines `iauSp00`, `iauPom00`

use crate::time::JulianDate;
use affn::Rotation3;
use qtty::*;

/// TIO locator s' (radians).
///
/// ```text
/// s' ≈ −47 μas × t
/// ```
///
/// where t is Julian centuries from J2000 on the TT scale.
///
/// ## References
/// * IERS Conventions (2010), eq. 5.13
/// * SOFA routine `iauSp00`
#[inline]
pub fn tio_locator_sp(jd_tt: JulianDate) -> f64 {
    let t = jd_tt.julian_centuries().value();
    let as2rad = std::f64::consts::PI / (180.0 * 3600.0);
    // −47 μas/century
    -47e-6 * as2rad * t
}

/// Polar motion matrix **W**.
///
/// Transforms from TIRS to ITRS (or equivalently, from the intermediate
/// frame to the Earth-fixed frame):
///
/// ```text
/// W = R₃(−s') · R₂(xₚ) · R₁(yₚ)
/// ```
///
/// `xp`, `yp`: pole coordinates in **radians**.
/// `sp`: TIO locator s' in **radians** (from [`tio_locator_sp`]).
///
/// ## References
/// * IERS Conventions (2010), §5.4.2
/// * SOFA routine `iauPom00`
pub fn polar_motion_matrix(xp: f64, yp: f64, sp: f64) -> Rotation3 {
    let (sxp, cxp) = xp.sin_cos();
    let (syp, cyp) = yp.sin_cos();
    let (ssp, csp) = sp.sin_cos();

    // R₃(−s')
    #[rustfmt::skip]
    let r3 = [
        [csp, ssp, 0.0],
        [-ssp, csp, 0.0],
        [0.0, 0.0, 1.0],
    ];

    // R₂(xₚ)
    #[rustfmt::skip]
    let r2 = [
        [cxp, 0.0, -sxp],
        [0.0, 1.0, 0.0],
        [sxp, 0.0, cxp],
    ];

    // R₁(yₚ)
    #[rustfmt::skip]
    let r1 = [
        [1.0, 0.0, 0.0],
        [0.0, cyp, syp],
        [0.0, -syp, cyp],
    ];

    // W = R₃(−s') · R₂(xₚ) · R₁(yₚ)
    // Compute R₂(xₚ) · R₁(yₚ) first, then R₃(−s') · result
    let mut tmp = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            tmp[i][j] = r2[i][0] * r1[0][j] + r2[i][1] * r1[1][j] + r2[i][2] * r1[2][j];
        }
    }

    let mut w = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            w[i][j] = r3[i][0] * tmp[0][j] + r3[i][1] * tmp[1][j] + r3[i][2] * tmp[2][j];
        }
    }

    Rotation3::from_matrix(w)
}

/// Convenience: compute W from pole coordinates (arcseconds) and Julian Date.
///
/// This handles the unit conversion and TIO locator computation.
#[inline]
pub fn polar_motion_matrix_from_arcsec(xp_as: f64, yp_as: f64, jd_tt: JulianDate) -> Rotation3 {
    let as2rad = std::f64::consts::PI / (180.0 * 3600.0);
    let sp = tio_locator_sp(jd_tt);
    polar_motion_matrix(xp_as * as2rad, yp_as * as2rad, sp)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tio_locator_at_j2000_is_zero() {
        let sp = tio_locator_sp(JulianDate::J2000);
        assert!(sp.abs() < 1e-15, "s' at J2000 should be ~0, got {}", sp);
    }

    #[test]
    fn tio_locator_is_small() {
        // After 25 years: s' ≈ −47e-6 × 0.25 arcsec ≈ −11.75 μas
        let jd = JulianDate::new(2_460_000.5); // ~2023
        let sp = tio_locator_sp(jd);
        let sp_uas = sp * 206_264_806_000.0; // rad → μas
        assert!(
            sp_uas.abs() < 50.0,
            "s' should be ~12 μas, got {} μas",
            sp_uas
        );
    }

    #[test]
    fn polar_motion_identity_with_zero_poles() {
        let w = polar_motion_matrix(0.0, 0.0, 0.0);
        let m = w.as_matrix();
        for i in 0..3 {
            assert!(
                (m[i][i] - 1.0).abs() < 1e-15,
                "W[{}][{}] = {}, expected 1",
                i, i, m[i][i]
            );
        }
    }

    #[test]
    fn polar_motion_is_proper_rotation() {
        let as2rad = std::f64::consts::PI / (180.0 * 3600.0);
        let w = polar_motion_matrix(0.1 * as2rad, 0.3 * as2rad, -1e-9);
        let m = w.as_matrix();

        // Check det ≈ 1
        let det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
        assert!(
            (det - 1.0).abs() < 1e-14,
            "det(W) = {}, expected ≈ 1",
            det
        );
    }

    #[test]
    fn polar_motion_small_angles() {
        // For small pole coordinates (~0.3″), the rotation should be tiny
        let w = polar_motion_matrix_from_arcsec(0.3, 0.2, JulianDate::new(2_460_000.5));
        let m = w.as_matrix();

        // Diagonal elements should be very close to 1
        for i in 0..3 {
            assert!(
                (m[i][i] - 1.0).abs() < 1e-10,
                "W[{}][{}] = {}, should be ≈ 1 for small poles",
                i, i, m[i][i]
            );
        }
    }
}
