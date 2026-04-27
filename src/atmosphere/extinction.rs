// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Beer-Lambert atmospheric transmission.

/// Slant-path transmission `T = exp(-X · τ)`.
///
/// `airmass` is the dimensionless geometric path length (see
/// [`crate::atmosphere::airmass`]) and `tau` is the per-unit-airmass
/// optical depth (e.g. Rayleigh + Mie + ozone).
#[inline]
pub fn transmission(tau: f64, airmass: f64) -> f64 {
    (-tau * airmass).exp()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn zero_tau_means_unit_transmission() {
        assert_eq!(transmission(0.0, 1.5), 1.0);
    }

    #[test]
    fn monotonically_decreases_in_tau() {
        let t1 = transmission(0.1, 1.0);
        let t2 = transmission(0.2, 1.0);
        assert!(t2 < t1);
    }

    #[test]
    fn monotonically_decreases_in_airmass() {
        let t1 = transmission(0.1, 1.0);
        let t2 = transmission(0.1, 2.0);
        assert!(t2 < t1);
    }
}
