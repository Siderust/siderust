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
//! All public inputs and outputs are typed [`qtty`] quantities.  The density
//! interface accepts a typed [`Kilometers`] altitude (consistent with
//! [`OrbitState`] position units) and returns a typed
//! [`KilogramsPerCubicMeter`].  The [`DensityProvider`] is part of the typed
//! public boundary of the force model.
//!
//! [`OrbitState`]: crate::astro::dynamics::state::OrbitState

use crate::qtty::{KilogramsPerCubicMeter, Kilometers};

/// Provider returning atmospheric mass density at a geodetic altitude.
///
/// Implement this trait to plug in any density model (exponential, constant,
/// NRLMSISE-00, …) into the drag force model.
pub trait DensityProvider: Send + Sync {
    /// Mass density `ρ` at the given geodetic `altitude`.
    fn density(&self, altitude: Kilometers) -> KilogramsPerCubicMeter;
}

/// Single-layer exponential atmosphere.
///
/// ```text
/// ρ(h) = ρ₀ · exp(−(h − h₀) / H)
/// ```
///
/// where:
/// - `ρ₀` is the reference density at reference altitude `h₀`,
/// - `H` is the scale height.
///
/// This profile is representative of a narrow altitude band; for wider ranges
/// a table-driven multi-layer model is more accurate.
#[derive(Debug, Clone, Copy)]
pub struct ExponentialAtmosphere {
    /// Reference density at `h0`.
    pub rho0: KilogramsPerCubicMeter,
    /// Reference altitude.
    pub h0: Kilometers,
    /// Atmospheric scale height.
    pub scale_height: Kilometers,
}

impl ExponentialAtmosphere {
    /// Approximate USSA-like values representative of ~500 km altitude.
    ///
    /// These numbers are indicative; production runs should use calibrated
    /// tables or a full atmosphere model.
    pub const LEO_500KM: Self = Self {
        rho0: KilogramsPerCubicMeter::new(6.967e-13),
        h0: Kilometers::new(500.0),
        scale_height: Kilometers::new(63.822),
    };
}

impl DensityProvider for ExponentialAtmosphere {
    #[inline]
    fn density(&self, altitude: Kilometers) -> KilogramsPerCubicMeter {
        let exponent = -(altitude.value() - self.h0.value()) / self.scale_height.value();
        KilogramsPerCubicMeter::new(self.rho0.value() * exponent.exp())
    }
}

/// Constant-density atmosphere.  Useful only for synthetic tests.
#[derive(Debug, Clone, Copy)]
pub struct ConstantDensity {
    /// Density returned for every altitude.
    pub rho: KilogramsPerCubicMeter,
}

impl DensityProvider for ConstantDensity {
    #[inline]
    fn density(&self, _altitude: Kilometers) -> KilogramsPerCubicMeter {
        self.rho
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exponential_decreases_with_altitude() {
        let atm = ExponentialAtmosphere::LEO_500KM;
        assert!(
            atm.density(Kilometers::new(400.0)).value()
                > atm.density(Kilometers::new(500.0)).value()
        );
        assert!(
            atm.density(Kilometers::new(500.0)).value()
                > atm.density(Kilometers::new(600.0)).value()
        );
    }

    #[test]
    fn constant_density_independent_of_altitude() {
        let rho = KilogramsPerCubicMeter::new(1.23e-12);
        let atm = ConstantDensity { rho };
        assert_eq!(atm.density(Kilometers::new(0.0)).value(), 1.23e-12);
        assert_eq!(atm.density(Kilometers::new(1_000.0)).value(), 1.23e-12);
    }
}
