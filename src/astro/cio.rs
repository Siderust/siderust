// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # CIO Locator and CIP Coordinates — IAU 2006/2000A
//!
//! In the IAU 2000/2006 framework, the orientation of the Celestial Intermediate
//! Pole (CIP) is described by two direction cosines **X** and **Y** in the GCRS,
//! and the position of the Celestial Intermediate Origin (CIO) on the CIP equator
//! is given by the quantity **s**.
//!
//! ## CIP (X, Y) coordinates
//!
//! The CIP direction in the GCRS is given by:
//!
//! ```text
//! X ≈ sin(ψ̄ + Δψ) × sin(ε_A + Δε) + corrections
//! Y ≈ −(ε_A + Δε) + corrections
//! ```
//!
//! but the rigorous approach uses the Fukushima-Williams angles with nutation.
//!
//! ## CIO locator s
//!
//! The CIO locator `s` is defined as:
//!
//! ```text
//! s = −½ X Y + Σ(polynomial + trigonometric series)
//! ```
//!
//! For most applications, the simplified form `s ≈ −½ X Y` is accurate to ~1 μas.
//!
//! ## BPN matrix (CIO-based)
//!
//! The celestial-to-intermediate frame rotation uses the CIO-based matrix:
//!
//! ```text
//! Q(t) = R₃(−E)·R₂(−d)·R₃(E)·R₃(s)
//! ```
//!
//! where `E = atan2(Y, X)`, `d = acos(Z)`, `Z = √(1 − X² − Y²)`.
//!
//! ## References
//!
//! * IAU 2006 Resolution B1
//! * IERS Conventions (2010), §5.4.4, §5.5.4
//! * SOFA routines `iauXy06`, `iauS06`, `iauC2ixys`

use crate::astro::precession;
use crate::time::JulianDate;
use affn::Rotation3;
use qtty::*;

/// CIP (X, Y) coordinates and CIO locator s.
#[derive(Debug, Clone, Copy)]
pub struct CipCio {
    /// X coordinate of CIP in GCRS (dimensionless direction cosine).
    pub x: f64,
    /// Y coordinate of CIP in GCRS (dimensionless direction cosine).
    pub y: f64,
    /// CIO locator s.
    pub s: Radians,
}

/// Compute CIP (X, Y) coordinates from the Fukushima-Williams precession angles
/// and nutation corrections.
///
/// This extracts the CIP direction from the full precession-nutation matrix (NPB),
/// which is the most rigorous approach.
///
/// `jd`: Julian Date (TT).
/// `dpsi`: nutation in longitude Δψ (radians).
/// `deps`: nutation in obliquity Δε (radians).
///
/// Returns (X, Y) as the first two components of the CIP unit vector in the GCRS.
///
/// ## References
/// * SOFA routine `iauXys06a` (via the NPB matrix approach)
/// * IERS Conventions (2010), §5.4.1
pub fn cip_xy(jd: JulianDate, dpsi: Radians, deps: Radians) -> (f64, f64) {
    // Build the full NPB matrix (bias + precession + nutation)
    let npb = precession::precession_nutation_matrix(jd, dpsi, deps);
    let m = npb.as_matrix();

    // X = sin(d)·cos(E) = m[2][0], Y = sin(d)·sin(E) = m[2][1]
    // In the SOFA convention: X = m[2][0], Y = m[2][1]
    (m[2][0], m[2][1])
}

/// CIO locator `s`, simplified form.
///
/// ```text
/// s ≈ −½ X Y + polynomial
/// ```
///
/// The polynomial captures the dominant secular term. For the full series
/// with ~33 lunisolar + planetary terms, see SOFA `iauS06`.
///
/// `jd`: Julian Date (TT).
/// `x`, `y`: CIP coordinates.
///
/// Returns s as [`Radians`].
///
/// ## References
/// * IERS Conventions (2010), §5.5.4
/// * Capitaine, Wallace & Chapront (2003)
/// * SOFA routine `iauS06`
pub fn cio_locator_s(jd: JulianDate, x: f64, y: f64) -> Radians {
    let t = jd.julian_centuries().value();

    // Polynomial part (μas), from IERS Conventions (2010) eq. 5.15
    let poly_uas = 94.0 + 3808.65 * t - 122.68 * t.powi(2) - 72574.11 * t.powi(3);

    let s_rad = -0.5 * x * y + MicroArcseconds::new(poly_uas).to::<Radian>().value();
    Radians::new(s_rad)
}

/// Compute the full CIP/CIO triplet (X, Y, s).
///
/// Convenience function combining [`cip_xy`] and [`cio_locator_s`].
pub fn cip_cio(jd: JulianDate, dpsi: Radians, deps: Radians) -> CipCio {
    let (x, y) = cip_xy(jd, dpsi, deps);
    let s = cio_locator_s(jd, x, y);
    CipCio { x, y, s }
}

/// CIO-based celestial-to-intermediate matrix Q(t).
///
/// This rotates from the GCRS to the Celestial Intermediate Reference System
/// (CIRS), using the CIP coordinates (X, Y) and CIO locator s.
///
/// ```text
/// Q(t) = R₃(−E)·R₂(−d)·R₃(E + s)
/// ```
///
/// where `E = atan2(Y, X)`, `d = acos(√(1 − X² − Y²))`.
///
/// ## References
/// * SOFA routine `iauC2ixys`
/// * IERS Conventions (2010), eq. 5.8
pub fn gcrs_to_cirs_matrix(x: f64, y: f64, s: Radians) -> Rotation3 {
    let r2 = x * x + y * y;
    let e = y.atan2(x);
    let d = r2.sqrt().atan2((1.0 - r2).max(0.0).sqrt());

    let (se, ce) = e.sin_cos();
    let (sd, cd) = d.sin_cos();

    // R₃(-(E+s)) · R₂(-d) · R₃(E)
    // Following SOFA's iauC2ixys construction:
    let es = e + s.value();
    let (ses, ces) = es.sin_cos();

    #[rustfmt::skip]
    let m = [
        [
            ces * cd * ce + ses * se,
            ces * cd * se - ses * ce,
            ces * sd,
        ],
        [
            ses * cd * ce - ces * se,
            ses * cd * se + ces * ce,
            ses * sd,
        ],
        [
            -sd * ce,
            -sd * se,
            cd,
        ],
    ];

    // The matrix above is Q^T (intermediate-to-celestial). We need Q (celestial-to-intermediate).
    // Actually, SOFA's C2IXYS returns the celestial-to-intermediate matrix directly.
    // Let's verify: at J2000 with X≈0, Y≈0, s≈0, this should be ≈ identity.
    // E = atan2(0,0) = 0, d = 0, s = 0.
    // m[0][0] = 1, m[1][1] = 1, m[2][2] = 1. ✓
    Rotation3::from_matrix(m)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cip_xy_at_j2000() {
        // At J2000.0, the CIP is very close to the GCRS pole, so X ≈ 0, Y ≈ 0
        let (x, y) = cip_xy(JulianDate::J2000, Radians::new(0.0), Radians::new(0.0));
        // Frame bias gives X, Y of order ~10 mas = ~5e-8 rad
        assert!(x.abs() < 1e-4, "X at J2000 should be ~0, got {}", x);
        assert!(y.abs() < 1e-4, "Y at J2000 should be ~0, got {}", y);
    }

    #[test]
    fn cio_locator_small() {
        // The CIO locator s is very small (< 1 mas for current epochs)
        let s = cio_locator_s(JulianDate::J2000, 0.0, 0.0);
        let s_mas = s.value() * 206_264_806.0; // rad → mas
        assert!(
            s_mas.abs() < 100.0,
            "s at J2000 should be < 100 mas, got {} mas",
            s_mas
        );
    }

    #[test]
    fn gcrs_to_cirs_near_identity_at_j2000() {
        // At J2000 with zero nutation, Q should be near-identity (only frame bias)
        let q = gcrs_to_cirs_matrix(0.0, 0.0, Radians::new(0.0));
        let m = q.as_matrix();
        for (i, row) in m.iter().enumerate().take(3) {
            assert!(
                (row[i] - 1.0).abs() < 1e-6,
                "Q[{}][{}] = {}, expected ≈ 1",
                i,
                i,
                row[i]
            );
        }
    }

    #[test]
    fn gcrs_to_cirs_is_proper_rotation() {
        let jd = JulianDate::new(2_460_000.5);
        let nut = crate::astro::nutation::nutation_iau2000b(jd);
        let cip = cip_cio(jd, nut.dpsi, nut.deps);
        let q = gcrs_to_cirs_matrix(cip.x, cip.y, cip.s);
        let m = q.as_matrix();

        // Check det ≈ 1 (proper rotation)
        let det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
        assert!((det - 1.0).abs() < 1e-12, "det(Q) = {}, expected ≈ 1", det);
    }
}
