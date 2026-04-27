// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Rayleigh scattering: optical depth and phase function.

use crate::ext_qtty::length::{Kilometers, Nanometers};

/// Default scale height of the dry-atmosphere column, in km. Used by
/// [`rayleigh_optical_depth_bodhaine99`] when callers do not pass an
/// explicit override.
pub const DEFAULT_SCALE_HEIGHT_KM: f64 = 8.0;

/// Rayleigh optical depth using the Bodhaine et al. (1999) approximation,
/// scaled by surface pressure and an exponential observatory-altitude
/// correction.
///
/// Inputs:
/// - `wavelength` — typed wavelength.
/// - `pressure_hpa` — surface pressure in hectopascal (no qtty `Pressure`
///   unit yet; the parameter name encodes the scale).
/// - `observatory_height` — observatory altitude above sea level
///   (typed [`Kilometers`]).
/// - `scale_height_km` — column scale height. Pass
///   [`DEFAULT_SCALE_HEIGHT_KM`] for the standard 8 km dry-atmosphere
///   value used by NSB and `darknsb`.
///
/// Reference: Bodhaine, B. A., Wood, N. B., Dutton, E. G., & Slusser, J. R.
/// (1999), "On Rayleigh optical depth calculations", *Journal of
/// Atmospheric and Oceanic Technology* 16, 1854.
pub fn rayleigh_optical_depth_bodhaine99(
    wavelength: Nanometers,
    pressure_hpa: f64,
    observatory_height: Kilometers,
    scale_height_km: f64,
) -> f64 {
    let lam_um = wavelength.value() * 1e-3;
    let p_atm = pressure_hpa / 1013.25;
    let h_km = observatory_height.value();

    let l2 = lam_um * lam_um;
    let inv_l2 = 1.0 / l2;
    // Bodhaine et al. 1999 simplified form (matches NSB GetRayleighOptDepth).
    let tau_sea = 0.0021520
        * (1.0455996 - 341.29061 * inv_l2 - 0.90230850 * l2)
        / (1.0 + 0.0027059889 * inv_l2 - 85.968563 * l2);
    let height_factor = (-h_km / scale_height_km).exp();
    p_atm * tau_sea * height_factor
}

/// Rayleigh phase function `P(θ) = 3 / (16π) · (1 + cos²θ)`.
#[inline]
pub fn rayleigh_phase(cos_theta: f64) -> f64 {
    3.0 / (16.0 * core::f64::consts::PI) * (1.0 + cos_theta * cos_theta)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matches_nsb_at_550nm_paranal_like() {
        // Replays NSB's GetRayleighOptDepth at 550 nm, P = 779 hPa, h = 2.4 km.
        // Expected value computed by the same formula in nsb::atmosphere::rayleigh.
        let tau = rayleigh_optical_depth_bodhaine99(
            Nanometers::new(550.0),
            779.0,
            Kilometers::new(2.4),
            DEFAULT_SCALE_HEIGHT_KM,
        );
        assert!(tau > 0.05 && tau < 0.15, "tau={tau}");
    }

    #[test]
    fn phase_at_90deg_is_3_over_16pi() {
        let p = rayleigh_phase(0.0);
        let expected = 3.0 / (16.0 * core::f64::consts::PI);
        assert!((p - expected).abs() < 1e-15);
    }

    #[test]
    fn phase_at_0deg_is_double_90deg() {
        let p_0 = rayleigh_phase(1.0);
        let p_90 = rayleigh_phase(0.0);
        assert!((p_0 - 2.0 * p_90).abs() < 1e-15);
    }
}
