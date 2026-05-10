// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Gravity field provider abstraction.
//!
//! ## Scientific scope
//!
//! High-fidelity orbit propagation and POD use spherical-harmonic
//! geopotential models (EGM2008, EIGEN-6C, …) expressed as normalised
//! Stokes coefficients `C_{n,m}` and `S_{n,m}` together with the field
//! constants `GM` and the equatorial reference radius `R`.
//!
//! This module defines:
//! - [`GravityConstants`] — `GM`, `R`, and the maximum degree stored.
//! - [`GravityFieldProvider`] — trait returning coefficients on demand.
//! - [`TwoBodyEarth`] — trivial degree-0 field for two-body propagation.
//!
//! Concrete implementations for full geopotential models (EGM2008, …)
//! live in downstream crates and plug in via `impl GravityFieldProvider`.

use crate::astro::dynamics::units::GravitationalParameter;
use crate::qtty::Kilometers;

/// Gravity-field constants.
#[derive(Debug, Clone, Copy)]
pub struct GravityConstants {
    /// `GM = G·M_central`.
    pub gm: GravitationalParameter,
    /// Equatorial reference radius of the field.
    pub radius: Kilometers,
    /// Maximum degree available in the model.
    pub max_degree: u16,
}

/// Provider returning normalised geopotential coefficients.
///
/// Implement this trait to plug in any gravity field model.
pub trait GravityFieldProvider: Send + Sync {
    /// Field constants: GM, equatorial radius, and maximum degree.
    fn constants(&self) -> GravityConstants;

    /// Normalised cosine coefficient `C_{n,m}`.
    fn c_nm(&self, n: u16, m: u16) -> f64;

    /// Normalised sine coefficient `S_{n,m}`.
    fn s_nm(&self, n: u16, m: u16) -> f64;
}

/// Trivial two-body Earth gravity field (degree 0 only).
///
/// Uses the standard WGS-84 values:
/// `GM = 398600.4418 km³/s²`, `R = 6378.137 km`.
#[derive(Debug, Clone, Copy)]
pub struct TwoBodyEarth;

impl GravityFieldProvider for TwoBodyEarth {
    fn constants(&self) -> GravityConstants {
        GravityConstants {
            gm: GravitationalParameter::new(398_600.441_8),
            radius: Kilometers::new(6_378.137),
            max_degree: 0,
        }
    }

    fn c_nm(&self, n: u16, _m: u16) -> f64 {
        if n == 0 { 1.0 } else { 0.0 }
    }

    fn s_nm(&self, _n: u16, _m: u16) -> f64 {
        0.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn two_body_earth_constants() {
        let g = TwoBodyEarth.constants();
        assert!((g.gm.value() - 398_600.441_8).abs() < 1e-3);
        assert!((g.radius.value() - 6_378.137).abs() < 1e-3);
        assert_eq!(g.max_degree, 0);
    }

    #[test]
    fn two_body_earth_coefficients() {
        assert_eq!(TwoBodyEarth.c_nm(0, 0), 1.0);
        assert_eq!(TwoBodyEarth.c_nm(2, 0), 0.0);
        assert_eq!(TwoBodyEarth.s_nm(2, 1), 0.0);
    }
}
