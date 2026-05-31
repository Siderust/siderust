// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Rayleigh scattering: optical depth and phase function
//!
//! ## Scientific scope
//!
//! Rayleigh scattering is the elastic scattering of light by atmospheric
//! molecules whose size is much smaller than the wavelength. It dominates
//! the molecular component of broadband optical extinction and gives the
//! sky its blue colour. This module provides:
//!
//! - the **Bodhaine et al. 1999** column-integrated Rayleigh optical depth
//!   parameterisation, scaled by surface pressure and corrected for
//!   observer altitude via an exponential atmosphere of fixed scale height;
//! - the **Rayleigh phase function** `P(θ) = (3 / 16π)·(1 + cos²θ)`.
//!
//! ## Technical scope
//!
//! - Wavelength is taken as a typed [`Nanometers`].
//! - Surface pressure is taken as a typed [`Hectopascals`].
//! - Observer altitude and scale height are typed [`Kilometers`].
//! - Optical depth is returned as [`OpticalDepths`].
//! - The phase function takes a typed scattering angle [`Radians`] and
//!   returns a dimensionless [`ScatteringFactor`]
//!   quantity.
//!
//! Internally the kernel still operates on `f64` for performance; the
//! conversion seam is one named binding per function.
//!
//! ## References
//!
//! - Bodhaine, B. A., Wood, N. B., Dutton, E. G., & Slusser, J. R. (1999).
//!   "On Rayleigh optical depth calculations". *J. Atmos. Oceanic Technol.*
//!   16, 1854.

use crate::atmosphere::ScatteringFactor;
use crate::ext_qtty::pressure::Hectopascals;
use crate::ext_qtty::Quantity;
use crate::qtty::unit::Micrometer;
use crate::qtty::{Kilometers, Nanometers, OpticalDepths, Radians};

/// Default scale height of the dry-atmosphere column (`8 km`).
///
/// Used by [`rayleigh_optical_depth_bodhaine99`] when callers do not pass
/// an explicit override; matches the standard NSB / `darknsb` value.
pub const DEFAULT_SCALE_HEIGHT: Kilometers = Kilometers::new(8.0);

const SEA_LEVEL_PRESSURE: Hectopascals = Hectopascals::new(1013.25);

/// Rayleigh optical depth using the Bodhaine et al. (1999) approximation,
/// scaled by surface pressure and an exponential observatory-altitude
/// correction.
///
/// # Inputs
///
/// - `wavelength` — typed wavelength.
/// - `surface_pressure` — surface pressure (typed [`Hectopascals`]).
/// - `observer_altitude` — observatory altitude above sea level
///   (typed [`Kilometers`]).
/// - `scale_height` — atmospheric column scale height (typed
///   [`Kilometers`]). Pass [`DEFAULT_SCALE_HEIGHT`] for the standard 8 km
///   value used by NSB and `darknsb`.
///
/// # Returns
///
/// Vertical Rayleigh optical depth at the supplied wavelength, as a typed
/// [`OpticalDepths`] value.
///
/// # References
///
/// Bodhaine, B. A., Wood, N. B., Dutton, E. G., & Slusser, J. R. (1999).
pub fn rayleigh_optical_depth_bodhaine99(
    wavelength: Nanometers,
    surface_pressure: Hectopascals,
    observer_altitude: Kilometers,
    scale_height: Kilometers,
) -> OpticalDepths {
    let lam_um = wavelength.to::<Micrometer>().value();
    let p_atm: f64 = surface_pressure / SEA_LEVEL_PRESSURE;
    let h_km = observer_altitude.value();
    let scale_km = scale_height.value();

    let l2 = lam_um * lam_um;
    let inv_l2 = 1.0 / l2;
    // Bodhaine et al. 1999 simplified form (matches NSB GetRayleighOptDepth).
    let tau_sea = 0.0021520 * (1.0455996 - 341.29061 * inv_l2 - 0.90230850 * l2)
        / (1.0 + 0.0027059889 * inv_l2 - 85.968563 * l2);
    let height_factor = (-h_km / scale_km).exp();
    OpticalDepths::new(p_atm * tau_sea * height_factor)
}

/// Rayleigh phase function `P(θ) = 3 / (16π) · (1 + cos²θ)`.
///
/// Returns a dimensionless [`ScatteringFactor`] quantity at the typed
/// scattering angle `theta`.
#[inline]
pub fn rayleigh_phase(theta: Radians) -> Quantity<ScatteringFactor> {
    let c = theta.cos();
    Quantity::<ScatteringFactor>::new(3.0 / (16.0 * core::f64::consts::PI) * (1.0 + c * c))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matches_nsb_at_550nm_paranal_like() {
        let tau = rayleigh_optical_depth_bodhaine99(
            Nanometers::new(550.0),
            Hectopascals::new(779.0),
            Kilometers::new(2.4),
            DEFAULT_SCALE_HEIGHT,
        );
        let v = tau.value();
        assert!(v > 0.05 && v < 0.15, "tau={v}");
    }

    #[test]
    fn phase_at_90deg_is_3_over_16pi() {
        let p = rayleigh_phase(Radians::new(core::f64::consts::FRAC_PI_2));
        let expected = 3.0 / (16.0 * core::f64::consts::PI);
        assert!((p.value() - expected).abs() < 1e-15);
    }

    #[test]
    fn phase_at_0deg_is_double_90deg() {
        let p_0 = rayleigh_phase(Radians::new(0.0));
        let p_90 = rayleigh_phase(Radians::new(core::f64::consts::FRAC_PI_2));
        assert!((p_0.value() - 2.0 * p_90.value()).abs() < 1e-15);
    }
}
