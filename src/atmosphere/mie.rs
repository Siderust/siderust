// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Mie / aerosol optical depth — power-law parameterization.

use crate::ext_qtty::length::Nanometers;

/// Aerosol optical-depth model parameters.
///
/// The optical depth is `τ_M(λ) = τ₀ · (λ / λ_ref)^α`. This is the
/// Patat (2011) power-law commonly used in modern Cerro Paranal sky
/// brightness models. Negative α encodes the usual aerosol behaviour
/// (more extinction at shorter wavelengths).
///
/// Reference: Patat, F. (2011), "The dancing sky: 6 years of night-sky
/// observations at Cerro Paranal", *Astronomy & Astrophysics* 527, A91.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MieParams {
    /// Optical depth at the reference wavelength.
    pub tau0: f64,
    /// Wavelength exponent `α` (typically negative).
    pub alpha: f64,
    /// Reference wavelength `λ_ref`.
    pub lambda_ref: Nanometers,
}

impl MieParams {
    /// Cerro Paranal default used by `darknsb` / NSB
    /// (`τ₀ = 0.05`, `α = -1.38`, `λ_ref = 550 nm`).
    pub const PARANAL: MieParams = MieParams {
        tau0: 0.05,
        alpha: -1.38,
        lambda_ref: Nanometers::new(550.0),
    };
}

/// Mie / aerosol optical depth at the given wavelength using the
/// power-law parameterization stored in `params`.
#[inline]
pub fn mie_optical_depth(params: &MieParams, wavelength: Nanometers) -> f64 {
    params.tau0 * (wavelength.value() / params.lambda_ref.value()).powf(params.alpha)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn at_lambda_ref_returns_tau0() {
        let p = MieParams::PARANAL;
        assert_eq!(mie_optical_depth(&p, Nanometers::new(550.0)), p.tau0);
    }

    #[test]
    fn shorter_wavelength_higher_tau() {
        let p = MieParams::PARANAL;
        let blue = mie_optical_depth(&p, Nanometers::new(400.0));
        let red = mie_optical_depth(&p, Nanometers::new(700.0));
        assert!(blue > p.tau0);
        assert!(red < p.tau0);
    }
}
