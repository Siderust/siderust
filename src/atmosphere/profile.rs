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
    /// Cerro Paranal default profile.
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

    /// Roque de los Muchachos Observatory, La Palma (ORM) profile.
    ///
    /// - Surface pressure: 744 hPa (source: `NSB_Utils.py:59`, `la_palma_pres = 744 hPa`;
    ///   consistent with `siderust::observatories::ROQUE_DE_LOS_MUCHACHOS`).
    /// - Altitude: 2.396 km AMSL (WGS84 height of Roque de los Muchachos, 2396 m).
    /// - Rayleigh scale height: [`DEFAULT_SCALE_HEIGHT_KM`] (8.0 km).
    /// - Mie aerosols: [`MieParams::LA_PALMA`] (τ₀ = 0.02 at 550 nm).
    ///
    /// See [`crate::observatories::ROQUE_DE_LOS_MUCHACHOS`] for the observatory entry.
    pub const LA_PALMA: AtmosphereProfile = AtmosphereProfile {
        surface_pressure: Hectopascals::new(744.0),
        observer_altitude: Kilometers::new(2.396),
        rayleigh_scale_height_km: DEFAULT_SCALE_HEIGHT_KM,
        mie_params: MieParams::LA_PALMA,
    };

    /// Mauna Kea Observatories (CFHT / Subaru) profile.
    ///
    /// - Surface pressure: 615 hPa. At 4207 m AMSL the ISA standard-atmosphere
    ///   pressure is ≈ 600–620 hPa; 615 hPa is consistent with CFHT site
    ///   documentation and Krisciunas (1990, PASP 102, 1235).
    /// - Altitude: 4.207 km AMSL (WGS84 height, `siderust::observatories::MAUNA_KEA`).
    /// - Rayleigh scale height: [`DEFAULT_SCALE_HEIGHT_KM`] (8.0 km).
    /// - Mie aerosols: [`MieParams::MAUNA_KEA`] (τ₀ = 0.03 at 550 nm).
    ///
    /// **Note:** The lower pressure at Mauna Kea reduces the Rayleigh optical depth
    /// relative to sea-level sites by the ratio 615/1013 ≈ 0.61, making the
    /// atmosphere unusually transparent at near-UV/blue wavelengths.
    ///
    /// See [`crate::observatories::MAUNA_KEA`] for the observatory entry.
    pub const MAUNA_KEA: AtmosphereProfile = AtmosphereProfile {
        surface_pressure: Hectopascals::new(615.0),
        observer_altitude: Kilometers::new(4.207),
        rayleigh_scale_height_km: DEFAULT_SCALE_HEIGHT_KM,
        mie_params: MieParams::MAUNA_KEA,
    };

    /// La Silla Observatory (ESO, Chile) profile.
    ///
    /// - Surface pressure: 760 hPa. At 2400 m AMSL the ISA standard-atmosphere
    ///   pressure is ≈ 756 hPa; 760 hPa is consistent with ESO La Silla
    ///   atmospheric monitoring and Burki et al. (1995, A&AS 112, 383).
    /// - Altitude: 2.400 km AMSL (WGS84 height,
    ///   `siderust::observatories::LA_SILLA_OBSERVATORY`).
    /// - Rayleigh scale height: [`DEFAULT_SCALE_HEIGHT_KM`] (8.0 km).
    /// - Mie aerosols: [`MieParams::LA_SILLA`] (τ₀ = 0.04 at 550 nm).
    ///
    /// See [`crate::observatories::LA_SILLA_OBSERVATORY`] for the observatory entry.
    pub const LA_SILLA: AtmosphereProfile = AtmosphereProfile {
        surface_pressure: Hectopascals::new(760.0),
        observer_altitude: Kilometers::new(2.400),
        rayleigh_scale_height_km: DEFAULT_SCALE_HEIGHT_KM,
        mie_params: MieParams::LA_SILLA,
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

    /// Site profiles must differ from PARANAL on at least one parameter.
    #[test]
    fn site_profiles_differ_from_paranal() {
        let p = AtmosphereProfile::PARANAL;
        let lp = AtmosphereProfile::LA_PALMA;
        let mk = AtmosphereProfile::MAUNA_KEA;
        let ls = AtmosphereProfile::LA_SILLA;

        // τ₀ differs for all three sites.
        assert_ne!(lp.mie_params.tau0, p.mie_params.tau0, "La Palma τ₀ must differ");
        assert_ne!(mk.mie_params.tau0, p.mie_params.tau0, "Mauna Kea τ₀ must differ");
        assert_ne!(ls.mie_params.tau0, p.mie_params.tau0, "La Silla τ₀ must differ");

        // Mauna Kea also has different pressure and altitude.
        assert_ne!(mk.surface_pressure.value(), p.surface_pressure.value(),
            "Mauna Kea pressure must differ from Paranal");
        assert_ne!(mk.observer_altitude.value(), p.observer_altitude.value(),
            "Mauna Kea altitude must differ from Paranal");
    }

    /// Optical depth at 550 nm / airmass 2 (altitude ≈ 30°) must be
    /// positive, finite, and ordered roughly by aerosol loading.
    #[test]
    fn site_profiles_optical_depth_ordering_at_30deg() {
        let lambda = Nanometers::new(550.0);
        let airmass = 2.0; // ~30° elevation

        let t_lp = AtmosphereProfile::LA_PALMA.transmission(lambda, airmass);
        let t_mk = AtmosphereProfile::MAUNA_KEA.transmission(lambda, airmass);
        let t_ls = AtmosphereProfile::LA_SILLA.transmission(lambda, airmass);
        let t_pn = AtmosphereProfile::PARANAL.transmission(lambda, airmass);

        // All transmissions must be in (0, 1).
        for (name, t) in [("La Palma", t_lp), ("Mauna Kea", t_mk),
                          ("La Silla", t_ls), ("Paranal", t_pn)] {
            assert!(t > 0.0 && t < 1.0 && t.is_finite(),
                "{name} transmission {t:.4} must be in (0,1)");
        }

        // Mauna Kea has lower pressure so its Rayleigh OD is substantially
        // reduced; it should therefore have higher overall transmission than
        // La Silla and Paranal despite being at a different absolute scale.
        assert!(t_mk > t_pn,
            "Mauna Kea should transmit more than Paranal (lower pressure + cleaner aerosol)");
        assert!(t_ls > t_pn,
            "La Silla should transmit more than Paranal (cleaner aerosol)");
        assert!(t_lp > t_pn,
            "La Palma should transmit more than Paranal (cleaner aerosol)");
    }
}
