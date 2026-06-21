// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

/// Pixel value for a stellar surface-brightness map.
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct StellarSurfaceBrightness {
    /// Integrated photon radiance in photons cm^-2 ns^-1 sr^-1.
    pub integrated_ph_cm2_ns_sr: f64,
    /// B-band tenth-magnitude-star surface brightness per square degree.
    pub b_s10: f64,
    /// V-band tenth-magnitude-star surface brightness per square degree.
    pub v_s10: f64,
}

impl StellarSurfaceBrightness {
    /// Return an empty pixel value.
    #[must_use]
    pub fn zero() -> Self {
        Self {
            integrated_ph_cm2_ns_sr: 0.0,
            b_s10: 0.0,
            v_s10: 0.0,
        }
    }

    /// Return true when all values are finite and non-negative.
    #[must_use]
    pub fn is_valid(self) -> bool {
        self.integrated_ph_cm2_ns_sr.is_finite()
            && self.b_s10.is_finite()
            && self.v_s10.is_finite()
            && self.integrated_ph_cm2_ns_sr >= 0.0
            && self.b_s10 >= 0.0
            && self.v_s10 >= 0.0
    }
}
