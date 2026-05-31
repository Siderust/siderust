// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Beer–Lambert atmospheric transmission
//!
//! ## Scientific scope
//!
//! Transmission `T = exp(−X · τ)` is the fraction of incident photons that
//! survive a slant path through an absorbing/scattering medium of vertical
//! optical depth `τ` traversed at airmass `X`. Combined with `airmass`
//! formulas and an `AtmosphereProfile` it yields the multiplicative
//! attenuation between the top of the atmosphere and an observer's
//! pupil.
//!
//! ## Technical scope
//!
//! Inputs are typed [`OpticalDepths`] and [`Airmasses`]; the result is a
//! typed [`Transmittances`] value in `[0, 1]`.
//!
//! ## References
//!
//! - Bouguer, P. (1729). *Essai d'optique sur la gradation de la lumière*.
//! - Lambert, J. H. (1760). *Photometria*.

use crate::atmosphere::Transmittances;
use crate::qtty::{Airmasses, OpticalDepths};

/// Slant-path transmission `T = exp(-X · τ)`.
///
/// `airmass` is the dimensionless geometric path length (see
/// [`mod@crate::atmosphere::airmass`]) and `tau` is the per-unit-airmass
/// vertical optical depth (e.g. Rayleigh + Mie + ozone).
#[inline]
pub fn transmission(tau: OpticalDepths, airmass: Airmasses) -> Transmittances {
    Transmittances::new((-tau.value() * airmass.value()).exp())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn zero_tau_means_unit_transmission() {
        assert_eq!(
            transmission(OpticalDepths::new(0.0), Airmasses::new(1.5)).value(),
            1.0
        );
    }

    #[test]
    fn monotonically_decreases_in_tau() {
        let t1 = transmission(OpticalDepths::new(0.1), Airmasses::new(1.0)).value();
        let t2 = transmission(OpticalDepths::new(0.2), Airmasses::new(1.0)).value();
        assert!(t2 < t1);
    }

    #[test]
    fn monotonically_decreases_in_airmass() {
        let t1 = transmission(OpticalDepths::new(0.1), Airmasses::new(1.0)).value();
        let t2 = transmission(OpticalDepths::new(0.1), Airmasses::new(2.0)).value();
        assert!(t2 < t1);
    }
}
