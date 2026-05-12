// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! [`GravityFieldProvider`] trait and the companion [`GravityConstants`] struct.

use crate::astro::dynamics::units::GravitationalParameter;
use crate::qtty::Kilometers;

/// Gravity-field constants (legacy helper, still used by old call sites).
#[derive(Debug, Clone, Copy)]
pub struct GravityConstants {
    /// `GM = G·M_central`.
    pub gm: GravitationalParameter,
    /// Equatorial reference radius of the field.
    pub radius: Kilometers,
    /// Maximum degree available in the model.
    pub max_degree: u16,
}

/// Provider of normalised spherical-harmonic geopotential coefficients.
///
/// ## API contract
///
/// * `c_normalized(n, m)` / `s_normalized(n, m)` return the **fully
///   normalised** (4π-normalised) Stokes coefficients `C̄ₙₘ` and `S̄ₙₘ`.
///   `C̄₀₀ = 1`, `S̄ₙₘ = 0` for all `m = 0`.
/// * Calls with `n > max_degree()` or `m > max_order()` return `0.0`.
/// * `gm()` and `reference_radius()` must be consistent with the field model.
///
/// ## Normalisation convention
///
/// The normalisation factor is
/// ```text
/// N_{n,m} = sqrt( (2n+1)(2 - δ_{m,0}) (n-m)! / (n+m)! )
/// ```
/// so that `C_{n,m}^{unnorm} = C̄_{n,m} · N_{n,m}`.
pub trait GravityFieldProvider: Send + Sync {
    /// `GM = G·M` (km³/s²).
    fn gm(&self) -> GravitationalParameter;

    /// Equatorial reference radius (km).
    fn reference_radius(&self) -> Kilometers;

    /// Maximum spherical-harmonic degree stored in this model.
    fn max_degree(&self) -> usize;

    /// Maximum spherical-harmonic order stored in this model.
    fn max_order(&self) -> usize;

    /// Normalised cosine coefficient `C̄_{n,m}`.
    fn c_normalized(&self, n: usize, m: usize) -> f64;

    /// Normalised sine coefficient `S̄_{n,m}`.
    fn s_normalized(&self, n: usize, m: usize) -> f64;

    // ------------------------------------------------------------------
    // Legacy helpers — kept so call sites that use the old API continue
    // to compile without modification.
    // ------------------------------------------------------------------

    /// Field constants packed into a legacy [`GravityConstants`].
    fn constants(&self) -> GravityConstants {
        GravityConstants {
            gm: self.gm(),
            radius: self.reference_radius(),
            max_degree: self.max_degree() as u16,
        }
    }

    /// Normalised cosine coefficient — legacy `u16` index variant.
    fn c_nm(&self, n: u16, m: u16) -> f64 {
        self.c_normalized(n as usize, m as usize)
    }

    /// Normalised sine coefficient — legacy `u16` index variant.
    fn s_nm(&self, n: u16, m: u16) -> f64 {
        self.s_normalized(n as usize, m as usize)
    }
}
