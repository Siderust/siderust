// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! [`GravityFieldProvider`] trait and coefficient-access interface.
//!
//! ## Scope
//!
//! Provides the [`GravityFieldProvider`] trait — an abstraction over normalised
//! spherical-harmonic geopotential coefficient sources — and the companion
//! [`GravityConstants`] helper struct.
//!
//! ## Equations
//!
//! The gravitational potential in spherical coordinates is:
//!
//! ```text
//! U = (GM/r) Σ_{n=0}^{N} (R/r)^n  Σ_{m=0}^{n}
//!       [ C̄_{nm} P̄_{nm}(sin φ) cos(mλ)
//!       + S̄_{nm} P̄_{nm}(sin φ) sin(mλ) ]
//! ```
//!
//! where C̄, S̄ are the **fully normalised** Stokes coefficients.
//!
//! ## Normalisation
//!
//! The 4π-normalisation factor is:
//! ```text
//! N_{n,m} = sqrt( (2n+1)(2 - δ_{m,0}) (n-m)! / (n+m)! )
//! C_{n,m}^{unnorm} = C̄_{n,m} · N_{n,m}
//! ```
//!
//! ## Provider contract
//!
//! - `c_normalized(n, m)` / `s_normalized(n, m)` return fully-normalised coefficients.
//! - `C̄₀₀ = 1`, `S̄_{nm} = 0` for all m = 0.
//! - Out-of-range calls return 0.0 (silent clipping).
//! - `gm()` and `reference_radius()` must be self-consistent.
//!
//! ## References
//!
//! * Montenbruck & Gill, *Satellite Orbits* (2001), §3.2.
//! * Vallado, *Fundamentals of Astrodynamics and Applications* (2013), §8.6.
//! * IERS Conventions (2010).

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
