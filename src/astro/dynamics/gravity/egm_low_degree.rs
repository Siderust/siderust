// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Built-in gravity field providers.
//!
//! ## Scope
//!
//! Provides two reference implementations of [`GravityFieldProvider`]:
//! [`TwoBodyEarth`] (degree 0, point-mass only) and [`LowDegreeEarth`]
//! (degree/order 4, based on EGM2008 normalised coefficients).
//!
//! ## Models
//!
//! | Type | Degree/Order | Description | Use case |
//! |------|-------------|-------------|----------|
//! | [`TwoBodyEarth`] | 0/0 | Two-body point-mass (no harmonics) | Baseline, validation |
//! | [`LowDegreeEarth`] | 4/4 | EGM2008 C/S through degree/order 4 | LEO primary gravity |
//!
//! ## Validity limits
//!
//! [`TwoBodyEarth`] is exact for a point-mass Earth.  [`LowDegreeEarth`]:
//! - Includes J2, J3, J4, and coupled m ≠ 0 harmonics.
//! - Ignores all coefficients beyond degree/order 4.
//! - Suitable for LEO; higher-degree fields may be needed for high-precision
//!   work or geosynchronous orbits.
//!
//! ## EGM2008 constants
//!
//! - GM = 398,600.4418 km³/s² (IERS 2010 value)
//! - R = 6,378.137 km (WGS-84 equatorial radius)
//!
//! ## References
//!
//! * Pavlis, N.K., et al. (2012), "The development and evaluation of the
//!   Earth Gravitational Model 2008 (EGM2008)", *J. Geophys. Res.*,
//!   117, B04406.
//! * IERS Conventions (2010).

use super::provider::GravityFieldProvider;
use crate::astro::dynamics::units::GravitationalParameter;
use crate::qtty::Kilometers;

// ---------------------------------------------------------------------------
// WGS-84 / EGM2008 constants
// ---------------------------------------------------------------------------

/// Earth GM (WGS-84), km³/s².
const GM_EARTH: f64 = 398_600.441_8;
/// Earth equatorial radius (WGS-84 / EGM2008), km.
const R_EARTH: f64 = 6_378.137;

// ---------------------------------------------------------------------------
// TwoBodyEarth
// ---------------------------------------------------------------------------

/// Trivial degree-0 Earth gravity field for two-body propagation.
///
/// All harmonic coefficients are zero except `C̄₀₀ = 1`.
#[derive(Debug, Clone, Copy)]
pub struct TwoBodyEarth;

impl GravityFieldProvider for TwoBodyEarth {
    fn gm(&self) -> GravitationalParameter {
        GravitationalParameter::new(GM_EARTH)
    }

    fn reference_radius(&self) -> Kilometers {
        Kilometers::new(R_EARTH)
    }

    fn max_degree(&self) -> usize {
        0
    }

    fn max_order(&self) -> usize {
        0
    }

    fn c_normalized(&self, n: usize, _m: usize) -> f64 {
        if n == 0 {
            1.0
        } else {
            0.0
        }
    }

    fn s_normalized(&self, _n: usize, _m: usize) -> f64 {
        0.0
    }
}

// ---------------------------------------------------------------------------
// LowDegreeEarth
// ---------------------------------------------------------------------------

/// Earth geopotential through degree/order 4 using EGM2008 normalised
/// Stokes coefficients.
///
/// The zonal `C̄₂₀` value is derived from the J₂ constant used in
/// [`J2`](crate::astro::dynamics::forces::J2) (`1.082 626 68 × 10⁻³`) so
/// that a [`Geopotential`](crate::astro::dynamics::forces::Geopotential)
/// truncated to `(2, 0)` reproduces the analytic J₂ model to machine
/// precision.
///
/// All other coefficients are from the public EGM2008 release.
#[derive(Debug, Clone, Copy)]
pub struct LowDegreeEarth;

/// EGM2008 normalised coefficients through degree 4.
///
/// Stored as `[(n, m, C̄_nm, S̄_nm)]`.
///
/// `C̄₂₀` is pinned to `−J2 / √5` where `J2 = 1.082_626_68e-3` to match
/// the [`J2`] force model exactly.
///
/// Other values: EGM2008 (Pavlis et al. 2012).
#[rustfmt::skip]
const EGM2008_COEFFS: &[(usize, usize, f64, f64)] = &[
    // degree 2
    (2, 0, -4.841_653_701_469_824e-4,  0.0),               // = -J2/√5, J2 = 1.082_626_68e-3
    (2, 1, -2.066_155_09e-10,         1.384_413_9e-9),
    (2, 2,  2.439_383_74e-6,         -1.400_273_8e-6),
    // degree 3
    (3, 0,  9.571_612_0e-7,           0.0),
    (3, 1,  2.030_462_0e-6,           2.484_348_0e-7),
    (3, 2,  9.047_878_0e-7,          -6.190_135_0e-7),
    (3, 3,  7.213_217_0e-7,           1.414_349_0e-6),
    // degree 4
    (4, 0,  5.399_737_0e-7,           0.0),
    (4, 1, -5.362_616_0e-7,          -4.735_903_0e-7),
    (4, 2,  3.505_016_0e-7,           6.624_480_0e-7),
    (4, 3,  9.908_551_0e-8,          -2.009_382_0e-7),
    (4, 4, -1.888_440_0e-8,           3.088_772_0e-8),
];

impl LowDegreeEarth {
    /// Look up a normalised coefficient pair from the embedded table.
    fn lookup(n: usize, m: usize) -> (f64, f64) {
        if n == 0 && m == 0 {
            return (1.0, 0.0);
        }
        for &(tn, tm, c, s) in EGM2008_COEFFS {
            if tn == n && tm == m {
                return (c, s);
            }
        }
        (0.0, 0.0)
    }
}

impl GravityFieldProvider for LowDegreeEarth {
    fn gm(&self) -> GravitationalParameter {
        GravitationalParameter::new(GM_EARTH)
    }

    fn reference_radius(&self) -> Kilometers {
        Kilometers::new(R_EARTH)
    }

    fn max_degree(&self) -> usize {
        4
    }

    fn max_order(&self) -> usize {
        4
    }

    fn c_normalized(&self, n: usize, m: usize) -> f64 {
        if n > self.max_degree() || m > n {
            return 0.0;
        }
        Self::lookup(n, m).0
    }

    fn s_normalized(&self, n: usize, m: usize) -> f64 {
        if n > self.max_degree() || m > n {
            return 0.0;
        }
        Self::lookup(n, m).1
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn two_body_earth_c00() {
        assert_eq!(TwoBodyEarth.c_normalized(0, 0), 1.0);
        assert_eq!(TwoBodyEarth.c_normalized(2, 0), 0.0);
        assert_eq!(TwoBodyEarth.s_normalized(2, 1), 0.0);
    }

    #[test]
    fn two_body_earth_gm_radius() {
        let g = TwoBodyEarth;
        assert!((g.gm().value() - 398_600.441_8).abs() < 1e-3);
        assert!((g.reference_radius().value() - 6_378.137).abs() < 1e-3);
        assert_eq!(g.max_degree(), 0);
    }

    #[test]
    fn low_degree_earth_c00() {
        assert_eq!(LowDegreeEarth.c_normalized(0, 0), 1.0);
    }

    #[test]
    fn low_degree_earth_c20_consistent_with_j2() {
        // C̄₂₀ × √5 should equal −J2 (= −1.082_626_68e-3)
        let c20 = LowDegreeEarth.c_normalized(2, 0);
        let j2_recovered = -c20 * 5.0_f64.sqrt();
        let j2_expected = 1.082_626_68e-3;
        assert!(
            (j2_recovered - j2_expected).abs() / j2_expected < 1e-7,
            "C̄₂₀ inconsistent with J2: got {j2_recovered:.9e}, want {j2_expected:.9e}"
        );
    }

    #[test]
    fn low_degree_earth_out_of_range_returns_zero() {
        assert_eq!(LowDegreeEarth.c_normalized(5, 0), 0.0);
        assert_eq!(LowDegreeEarth.s_normalized(2, 5), 0.0);
    }
}
