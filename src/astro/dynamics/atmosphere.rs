// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Atmospheric density models for satellite drag.
//!
//! ## Scientific scope
//!
//! Satellite drag acceleration requires the atmospheric mass density `ρ`
//! (kg/m³) at the spacecraft's geodetic altitude.  This module provides a
//! trait-based abstraction over density sources so that force-model
//! implementations remain independent of the atmosphere backend.
//!
//! The default [`ExponentialAtmosphere`] is a single-layer exponential
//! profile — cheap, analytic, adequate for LEO regression tests and
//! short-arc propagation, but not suitable for operational POD.  More
//! realistic backends (NRLMSISE-00, DTM2000) can be plugged in by
//! implementing [`DensityProvider`] in downstream crates.
//!
//! ## Technical scope
//!
//! The density interface accepts a typed [`Kilometers`] altitude, consistent
//! with the [`OrbitState`] position units.  The public boundary of the
//! *force model* is fully typed; the density provider is part of that
//! boundary.
//!
//! [`OrbitState`]: crate::astro::dynamics::state::OrbitState

use crate::qtty::Kilometers;

/// Provider returning atmospheric mass density at a geodetic altitude.
///
/// Implement this trait to plug in any density model (exponential, constant,
/// NRLMSISE-00, …) into the drag force model.
pub trait DensityProvider: Send + Sync {
    /// Mass density `ρ` in **kg/m³** at the given geodetic altitude.
    fn density_kg_m3(&self, altitude: Kilometers) -> f64;
}

/// Single-layer exponential atmosphere.
///
/// ```text
/// ρ(h) = ρ₀ · exp(−(h − h₀) / H)
/// ```
///
/// where:
/// - `ρ₀` is the reference density (kg/m³) at reference altitude `h₀` (km),
/// - `H` is the scale height (km).
///
/// This profile is representative of a narrow altitude band; for wider ranges
/// a table-driven multi-layer model is more accurate.
#[derive(Debug, Clone, Copy)]
pub struct ExponentialAtmosphere {
    /// Reference density at `h0`, kg/m³.
    pub rho0_kg_m3: f64,
    /// Reference altitude, km.
    pub h0_km: f64,
    /// Scale height, km.
    pub scale_height_km: f64,
}

impl ExponentialAtmosphere {
    /// Approximate USSA-like values representative of ~500 km altitude.
    ///
    /// These numbers are indicative; production runs should use calibrated
    /// tables or a full atmosphere model.
    pub const LEO_500KM: Self = Self {
        rho0_kg_m3: 6.967e-13,
        h0_km: 500.0,
        scale_height_km: 63.822,
    };
}

impl DensityProvider for ExponentialAtmosphere {
    #[inline]
    fn density_kg_m3(&self, altitude: Kilometers) -> f64 {
        self.rho0_kg_m3 * (-(altitude.value() - self.h0_km) / self.scale_height_km).exp()
    }
}

/// Constant-density atmosphere.  Useful only for synthetic tests.
#[derive(Debug, Clone, Copy)]
pub struct ConstantDensity {
    /// Density returned for every altitude, kg/m³.
    pub rho: f64,
}

impl DensityProvider for ConstantDensity {
    #[inline]
    fn density_kg_m3(&self, _altitude: Kilometers) -> f64 {
        self.rho
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exponential_decreases_with_altitude() {
        use crate::qtty::Kilometers;
        let atm = ExponentialAtmosphere::LEO_500KM;
        assert!(atm.density_kg_m3(Kilometers::new(400.0)) > atm.density_kg_m3(Kilometers::new(500.0)));
        assert!(atm.density_kg_m3(Kilometers::new(500.0)) > atm.density_kg_m3(Kilometers::new(600.0)));
    }

    #[test]
    fn constant_density_independent_of_altitude() {
        use crate::qtty::Kilometers;
        let atm = ConstantDensity { rho: 1.23e-12 };
        assert_eq!(atm.density_kg_m3(Kilometers::new(0.0)), 1.23e-12);
        assert_eq!(atm.density_kg_m3(Kilometers::new(1000.0)), 1.23e-12);
    }
}
