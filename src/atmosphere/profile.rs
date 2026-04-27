// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Vertically-integrated atmosphere profile: Rayleigh + Mie optical depth.

use crate::ext_qtty::length::{Kilometers, Nanometers};
use crate::ext_qtty::pressure::Hectopascals;
use crate::atmosphere::extinction::transmission as ext_transmission;
use crate::atmosphere::mie::{mie_optical_depth, MieParams};
use crate::atmosphere::rayleigh::{
    rayleigh_optical_depth_bodhaine99, DEFAULT_SCALE_HEIGHT_KM,
};

/// Vertically-integrated atmosphere model used to compute the total
/// dry/aerosol optical depth at a given wavelength for a fixed observer
/// altitude and surface pressure.
///
/// # Ozone
/// Ozone transmittance is **not** part of this aggregate τ.
/// Apply it separately as a multiplicative transmittance via
/// `siderust::atmosphere::ozone::transmission_table()` for parity
/// with the NSB / Patat 2008 et al. workflow, where ozone extinction
/// is handled as an independent spectral factor.
#[derive(Clone, Debug)]
pub struct AtmosphereProfile {
    /// Observer surface pressure.
    pub surface_pressure: Hectopascals,
    /// Observer altitude above sea level.
    pub observer_altitude: Kilometers,
    /// Rayleigh atmosphere scale height (King-correction reference), in km.
    pub rayleigh_scale_height_km: f64,
    /// Mie aerosol parameters (τ₀, α, λ_ref).
    pub mie_params: MieParams,
}

impl AtmosphereProfile {
    /// Cerro Paranal / La Palma default profile.
    ///
    /// - Surface pressure: 744 hPa (standard Paranal value).
    /// - Altitude: 2.635 km AMSL — geodetic height of Cerro Paranal
    ///   (`siderust::observatories::EL_PARANAL`, 2635 m).
    /// - Rayleigh scale height: [`DEFAULT_SCALE_HEIGHT_KM`] (8.0 km).
    /// - Mie aerosols: [`MieParams::PARANAL`].
    pub const PARANAL: AtmosphereProfile = AtmosphereProfile {
        surface_pressure: Hectopascals::new(744.0),
        observer_altitude: Kilometers::new(2.635),
        rayleigh_scale_height_km: DEFAULT_SCALE_HEIGHT_KM,
        mie_params: MieParams::PARANAL,
    };

    /// Total optical depth at wavelength `lambda` (Rayleigh + Mie).
    ///
    /// Ozone is **not** included; apply `ozone::transmission_table()`
    /// separately as a multiplicative transmittance.
    pub fn optical_depth(&self, lambda: Nanometers) -> f64 {
        let tau_r = rayleigh_optical_depth_bodhaine99(
            lambda,
            self.surface_pressure.value(),
            self.observer_altitude,
            self.rayleigh_scale_height_km,
        );
        let tau_m = mie_optical_depth(&self.mie_params, lambda);
        tau_r + tau_m
    }

    /// Multiplicative transmission `exp(-airmass · τ)` at wavelength `lambda`.
    ///
    /// Convenience wrapper around [`crate::atmosphere::extinction::transmission`].
    /// Ozone is **not** included in τ; see [`Self::optical_depth`].
    pub fn transmission(&self, lambda: Nanometers, airmass: f64) -> f64 {
        ext_transmission(self.optical_depth(lambda), airmass)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Verify that `AtmosphereProfile::PARANAL.optical_depth(550 nm)` is
    /// bit-identical to the value produced by calling the underlying
    /// Rayleigh and Mie kernels directly with the same inputs — the same
    /// chain that NSB `extinction::optical_depth(550.0, 744.0, 2.635)`
    /// would execute.
    #[test]
    fn paranal_550nm_matches_manual_chain() {
        let profile_tau = AtmosphereProfile::PARANAL.optical_depth(Nanometers::new(550.0));

        let tau_r = rayleigh_optical_depth_bodhaine99(
            Nanometers::new(550.0),
            744.0,
            Kilometers::new(2.635),
            DEFAULT_SCALE_HEIGHT_KM,
        );
        let tau_m = mie_optical_depth(&MieParams::PARANAL, Nanometers::new(550.0));

        assert_eq!(profile_tau, tau_r + tau_m,
            "AtmosphereProfile must be bit-identical to the manual chain");
    }

    /// Custom profile: verify optical_depth equals the sum of Rayleigh + Mie
    /// parts computed independently.
    #[test]
    fn custom_profile_is_sum_of_parts() {
        use crate::ext_qtty::pressure::Hectopascals;

        let custom_mie = MieParams {
            tau0: 0.03,
            alpha: -1.2,
            lambda_ref: Nanometers::new(500.0),
        };
        let profile = AtmosphereProfile {
            surface_pressure: Hectopascals::new(900.0),
            observer_altitude: Kilometers::new(0.5),
            rayleigh_scale_height_km: 8.0,
            mie_params: custom_mie,
        };
        let lambda = Nanometers::new(450.0);
        let tau = profile.optical_depth(lambda);

        let tau_r = rayleigh_optical_depth_bodhaine99(
            lambda,
            900.0,
            Kilometers::new(0.5),
            8.0,
        );
        let tau_m = mie_optical_depth(&custom_mie, lambda);

        assert_eq!(tau, tau_r + tau_m,
            "optical_depth must equal tau_r + tau_m for custom profile");
    }

    #[test]
    fn transmission_decreases_with_airmass() {
        let t1 = AtmosphereProfile::PARANAL.transmission(Nanometers::new(550.0), 1.0);
        let t2 = AtmosphereProfile::PARANAL.transmission(Nanometers::new(550.0), 2.0);
        assert!(t2 < t1, "transmission must decrease with airmass");
    }
}
