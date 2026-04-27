// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Atmospheric optics primitives.
//!
//! Available behind the `atmosphere` feature. This module collects
//! citation-backed, typed-input helpers for:
//!
//! - [`airmass`] — the geometric path-length multiplier through a
//!   plane-parallel or curved-Earth atmosphere, with named formulas
//!   covering [Plane-parallel](airmass::AirmassFormula::PlaneParallel),
//!   [Young 1994](airmass::AirmassFormula::Young1994),
//!   [Rozenberg 1966](airmass::AirmassFormula::Rozenberg1966), and
//!   [Krisciunas & Schaefer 1991](airmass::AirmassFormula::KrisciunasSchaefer1991).
//! - [`rayleigh`] — Rayleigh optical depth (Bodhaine et al. 1999) and
//!   the Rayleigh phase function.
//! - [`mie`] — Mie / aerosol optical depth via the Patat 2011 power-law
//!   parameterization, with a `MieParams::PARANAL` preset.
//! - [`extinction`] — Beer-Lambert transmission `exp(-airmass · τ)`.
//! - [`profile`] — [`AtmosphereProfile`]: first-class aggregate of
//!   Rayleigh + Mie optical depth for a fixed observer site. Ozone is
//!   kept separate; apply `ozone::transmission_table()` independently.

pub mod airmass;
pub mod extinction;
pub mod mie;
pub mod profile;
pub mod rayleigh;
#[cfg(feature = "spectra")]
pub mod ozone;

pub use airmass::{airmass, AirmassFormula};
pub use extinction::transmission;
pub use mie::{mie_optical_depth, MieParams};
pub use profile::AtmosphereProfile;
pub use rayleigh::{rayleigh_optical_depth_bodhaine99, rayleigh_phase};
#[cfg(feature = "spectra")]
pub use ozone::{transmission_table, Transmittance};
