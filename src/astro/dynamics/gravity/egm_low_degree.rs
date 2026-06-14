// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Built-in Earth gravity-field providers.

use principia::{GravityFieldProvider, PrincipiaError};

use crate::astro::dynamics::units::GravitationalParameter;
use crate::qtty::Kilometers;

use crate::archive::gravity::tables::EGM2008_COEFFS;

/// Earth GM (WGS-84), km³/s².
const GM_EARTH_KM3_S2: f64 = 398_600.441_8;
/// Earth equatorial radius (WGS-84 / EGM2008), km.
const R_EARTH_KM: f64 = 6_378.137;

/// Trivial degree-0 Earth gravity field for two-body propagation.
#[derive(Debug, Clone, Copy)]
pub struct TwoBodyEarth;

impl GravityFieldProvider for TwoBodyEarth {
    fn mu(&self) -> GravitationalParameter {
        GravitationalParameter::new(GM_EARTH_KM3_S2)
    }

    fn reference_radius(&self) -> Kilometers {
        Kilometers::new(R_EARTH_KM)
    }

    fn max_degree(&self) -> usize {
        0
    }

    fn max_order(&self) -> usize {
        0
    }

    fn c_normalized(&self, n: usize, _m: usize) -> Result<f64, PrincipiaError> {
        Ok(if n == 0 { 1.0 } else { 0.0 })
    }

    fn s_normalized(&self, _n: usize, _m: usize) -> Result<f64, PrincipiaError> {
        Ok(0.0)
    }
}

/// Earth geopotential through degree/order 4 using EGM2008 normalized coefficients.
#[derive(Debug, Clone, Copy)]
pub struct LowDegreeEarth;

impl LowDegreeEarth {
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
    fn mu(&self) -> GravitationalParameter {
        GravitationalParameter::new(GM_EARTH_KM3_S2)
    }

    fn reference_radius(&self) -> Kilometers {
        Kilometers::new(R_EARTH_KM)
    }

    fn max_degree(&self) -> usize {
        4
    }

    fn max_order(&self) -> usize {
        4
    }

    fn c_normalized(&self, n: usize, m: usize) -> Result<f64, PrincipiaError> {
        Ok(if n > self.max_degree() || m > n {
            0.0
        } else {
            Self::lookup(n, m).0
        })
    }

    fn s_normalized(&self, n: usize, m: usize) -> Result<f64, PrincipiaError> {
        Ok(if n > self.max_degree() || m > n {
            0.0
        } else {
            Self::lookup(n, m).1
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn two_body_earth_c00() {
        assert_eq!(TwoBodyEarth.c_normalized(0, 0).unwrap(), 1.0);
        assert_eq!(TwoBodyEarth.c_normalized(2, 0).unwrap(), 0.0);
    }

    #[test]
    fn low_degree_earth_c20_consistent_with_j2() {
        let c20 = LowDegreeEarth.c_normalized(2, 0).unwrap();
        let j2_recovered = -c20 * 5.0_f64.sqrt();
        let j2_expected = 1.082_626_68e-3;
        assert!((j2_recovered - j2_expected).abs() / j2_expected < 1e-7);
    }
}
