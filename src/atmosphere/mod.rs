// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Atmospheric optics primitives.
//!
//! Available behind the `atmosphere` feature. This module collects
//! citation-backed, typed-input helpers for:
//!
//! - [`airmass`] — the geometric path-length multiplier through a
//!   plane-parallel or curved-Earth atmosphere, with compile-time formula
//!   selectors covering [`PlaneParallel`](airmass::PlaneParallel),
//!   [`Young1994`](airmass::Young1994),
//!   [`Rozenberg1966`](airmass::Rozenberg1966), and
//!   [`KrisciunasSchaefer1991`](airmass::KrisciunasSchaefer1991).
//! - [`rayleigh`] — Rayleigh optical depth (Bodhaine et al. 1999) and
//!   the Rayleigh phase function.
//! - [`mie`] — Mie / aerosol optical depth via the Patat 2011 power-law
//!   parameterization, with a `MieParams::PARANAL` preset.
//! - [`airglow`] — Van Rhijn airglow emission-layer geometry.
//! - [`scattering`] — phase-function traits and Rayleigh/tabulated helpers.
//! - [`extinction`] — Beer-Lambert transmission `exp(-airmass · τ)`.
//! - [`profile`] — [`AtmosphereProfile`]: first-class aggregate of
//!   Rayleigh + Mie optical depth for a fixed observer site. Ozone is
//!   kept separate; apply `ozone::transmission_table()` independently.

pub mod airglow;
pub mod airmass;
pub mod extinction;
pub mod mie;
#[cfg(feature = "spectra")]
pub mod ozone;
pub mod profile;
pub mod rayleigh;
pub mod scattering;

pub use airglow::{van_rhijn_factor, van_rhijn_factor_with_radius};
pub use airmass::{
    airmass, AirmassFormula, DefaultAirmassFormula, Formula, KrisciunasSchaefer1991,
    PlaneParallel, Rozenberg1966, Young1994,
};
pub use extinction::transmission;
pub use mie::{mie_optical_depth, MieParams};
#[cfg(feature = "spectra")]
pub use ozone::{transmission_table, Transmittance};
pub use profile::AtmosphereProfile;
pub use rayleigh::{rayleigh_optical_depth_bodhaine99, rayleigh_phase};
#[cfg(feature = "tables")]
pub use scattering::TabulatedPhaseFunction;
pub use scattering::{PhaseFunction, RayleighPhaseFunction, ScatteringFactor};
