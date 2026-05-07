// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Site-level atmospheric profiles (Rayleigh + Mie)
//!
//! ## Scientific scope
//!
//! Bundles the few parameters needed to evaluate a *clear-sky* atmospheric
//! optical depth at a specific observing site:
//!
//! - **Surface pressure** at the observatory altitude (Rayleigh column).
//! - **Observer altitude** above sea level (Rayleigh exponential decay).
//! - **Rayleigh scale height** (typically ≈ 8 km).
//! - **Mie aerosol parameters** ([`MieParams`]) – reference depth `τ₀`,
//!   reference wavelength `λ₀`, Ångström exponent `α`.
//!
//! The total Rayleigh + Mie optical depth at wavelength `λ` is the simple
//! sum of the two terms — a physically transparent first-order model
//! adequate for sky-brightness, photometric-extinction and exposure-time
//! work in clear conditions. Ozone, water vapour and trace species are
//! handled separately (see [`crate::atmosphere::ozone`]).
//!
//! ## Technical scope
//!
//! - All public fields and constructors take/return typed [`qtty`](crate::qtty)
//!   quantities — no raw `f64` units cross this module's boundary.
//! - Four ready-to-use observatory presets
//!   ([`AtmosphereProfile::PARANAL`], [`AtmosphereProfile::LA_PALMA`],
//!   [`AtmosphereProfile::MAUNA_KEA`], [`AtmosphereProfile::LA_SILLA`])
//!   are `const`-constructed from typed literals.
//! - [`AtmosphereProfile::optical_depth`] returns an [`OpticalDepths`]
//!   built from the typed sum of Rayleigh ([`rayleigh_optical_depth_bodhaine99`])
//!   and Mie ([`mie_optical_depth`]) contributions.
//! - [`AtmosphereProfile::transmission`] reuses
//!   [`crate::atmosphere::extinction::transmission`] and returns a
//!   typed [`Transmittances`].
//!
//! ## References
//!
//! - Bodhaine, B. A., Wood, N. B., Dutton, E. G., Slusser, J. R. (1999).
//!   "On Rayleigh Optical Depth Calculations". *J. Atmos. Ocean. Tech.*
//!   16, 1854.
//! - Patat, F., et al. (2011). "Optical atmospheric extinction over Cerro
//!   Paranal". *A&A* 527, A91.

use crate::atmosphere::extinction::transmission;
use crate::atmosphere::mie::{mie_optical_depth, MieParams};
use crate::atmosphere::rayleigh::{rayleigh_optical_depth_bodhaine99, DEFAULT_SCALE_HEIGHT};
use crate::atmosphere::Transmittances;
use crate::qtty::{Airmasses, Hectopascals, Kilometers, Nanometers, OpticalDepths};

/// Site-specific atmospheric profile combining Rayleigh + Mie inputs.
///
/// Construct one of the presets ([`EL_PARANAL`](Self::EL_PARANAL),
/// [`ROQUE_DE_LOS_MUCHACHOS`](Self::ROQUE_DE_LOS_MUCHACHOS),
/// [`MAUNA_KEA`](Self::MAUNA_KEA), [`LA_SILLA`](Self::LA_SILLA)) or build
/// a custom profile literal.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AtmosphereProfile {
    /// Atmospheric pressure at the observer (column reference for Rayleigh).
    pub surface_pressure: Hectopascals,
    /// Observer altitude above mean sea level.
    pub observer_altitude: Kilometers,
    /// Rayleigh exponential scale height (≈ 8 km in standard atmospheres).
    pub rayleigh_scale_height: Kilometers,
    /// Aerosol model parameters (`τ₀`, `λ₀`, `α`).
    pub mie_params: MieParams,
}

impl AtmosphereProfile {
    /// VLT — Cerro Paranal, Chile (≈ 2635 m).
    ///
    /// See [`crate::observatories::EL_PARANAL`].
    pub const EL_PARANAL: AtmosphereProfile = AtmosphereProfile {
        surface_pressure: Hectopascals::new(744.0),
        observer_altitude: Kilometers::new(2.635),
        rayleigh_scale_height: DEFAULT_SCALE_HEIGHT,
        mie_params: MieParams::PARANAL,
    };

    /// Roque de los Muchachos — La Palma, Spain (≈ 2396 m).
    ///
    /// See [`crate::observatories::ROQUE_DE_LOS_MUCHACHOS`].
    pub const ROQUE_DE_LOS_MUCHACHOS: AtmosphereProfile = AtmosphereProfile {
        surface_pressure: Hectopascals::new(770.0),
        observer_altitude: Kilometers::new(2.396),
        rayleigh_scale_height: DEFAULT_SCALE_HEIGHT,
        mie_params: MieParams::LA_PALMA,
    };

    /// Mauna Kea, Hawaiʻi, USA (≈ 4205 m).
    ///
    /// See [`crate::observatories::MAUNA_KEA`].
    pub const MAUNA_KEA: AtmosphereProfile = AtmosphereProfile {
        surface_pressure: Hectopascals::new(615.0),
        observer_altitude: Kilometers::new(4.205),
        rayleigh_scale_height: DEFAULT_SCALE_HEIGHT,
        mie_params: MieParams::MAUNA_KEA,
    };

    /// La Silla, Chile (≈ 2400 m).
    ///
    /// See [`crate::observatories::LA_SILLA`].
    pub const LA_SILLA: AtmosphereProfile = AtmosphereProfile {
        surface_pressure: Hectopascals::new(770.0),
        observer_altitude: Kilometers::new(2.400),
        rayleigh_scale_height: DEFAULT_SCALE_HEIGHT,
        mie_params: MieParams::LA_SILLA,
    };

    /// Total clear-sky optical depth at wavelength `lambda`:
    /// `τ_total(λ) = τ_Rayleigh(λ) + τ_Mie(λ)`.
    pub fn optical_depth(&self, lambda: Nanometers) -> OpticalDepths {
        let tau_r = rayleigh_optical_depth_bodhaine99(
            lambda,
            self.surface_pressure,
            self.observer_altitude,
            self.rayleigh_scale_height,
        );
        let tau_m = mie_optical_depth(&self.mie_params, lambda);
        tau_r + tau_m
    }

    /// Transmittance along a slant path of the given `airmass` at wavelength
    /// `lambda`: `T = exp(-τ_total · X)`.
    pub fn transmission(&self, lambda: Nanometers, airmass: Airmasses) -> Transmittances {
        transmission(self.optical_depth(lambda), airmass)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn presets_have_typed_units() {
        let p = AtmosphereProfile::EL_PARANAL;
        assert_eq!(p.surface_pressure.value(), 744.0);
        assert_eq!(p.observer_altitude.value(), 2.635);
        assert_eq!(p.rayleigh_scale_height.value(), 8.0);
    }

    #[test]
    fn optical_depth_is_positive_and_decreases_with_wavelength() {
        let p = AtmosphereProfile::EL_PARANAL;
        let blue = p.optical_depth(Nanometers::new(400.0));
        let red = p.optical_depth(Nanometers::new(700.0));
        assert!(blue.value() > 0.0);
        assert!(red.value() > 0.0);
        assert!(blue.value() > red.value(), "Rayleigh + Mie both fall with λ");
    }

    #[test]
    fn transmission_decreases_with_airmass() {
        let p = AtmosphereProfile::EL_PARANAL;
        let lambda = Nanometers::new(550.0);
        let t1 = p.transmission(lambda, Airmasses::new(1.0));
        let t2 = p.transmission(lambda, Airmasses::new(2.0));
        assert!(t1.value() > t2.value());
        assert!(t1.value() > 0.0 && t1.value() < 1.0);
    }

    #[test]
    fn higher_altitude_lowers_rayleigh_contribution() {
        // Mauna Kea (higher altitude, lower surface pressure) should give
        // less Rayleigh extinction at the same λ than Paranal.
        let lambda = Nanometers::new(500.0);
        let mk = AtmosphereProfile::MAUNA_KEA.optical_depth(lambda);
        let pa = AtmosphereProfile::EL_PARANAL.optical_depth(lambda);
        assert!(
            mk.value() < pa.value(),
            "Mauna Kea τ ({}) should be smaller than Paranal τ ({})",
            mk.value(),
            pa.value(),
        );
    }
}
