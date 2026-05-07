// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Atmospheric optics primitives
//!
//! ## Scientific scope
//!
//! Citation-backed, typed-input helpers for the line-of-sight effects of
//! Earth's atmosphere on a celestial source: airmass, Rayleigh / Mie
//! scattering optical depth, ozone band absorption, airglow geometry, and
//! Beer–Lambert extinction. These primitives compose the pipeline that
//! converts a physical irradiance at the top of the atmosphere into an
//! observed irradiance at a ground-based site.
//!
//! ## Technical scope
//!
//! Every public entry point in this module exchanges typed [`qtty`]
//! quantities — angles in [`Radians`](crate::qtty::Radians), wavelengths in
//! [`Nanometers`](crate::qtty::Nanometers), pressures in
//! [`Hectopascals`](crate::qtty::Hectopascals), heights in
//! [`Kilometers`](crate::qtty::Kilometers), airmass as
//! [`Airmasses`](crate::qtty::Airmasses), optical depth as
//! [`OpticalDepths`](crate::qtty::OpticalDepths), and transmittance as
//! [`Transmittance`]. The conversion seam between typed inputs and the raw
//! `f64` numerical kernels is a single named binding inside each function;
//! no `f64` leaks to the public API.
//!
//! Submodules:
//! - [`airmass`] — geometric path-length multipliers
//!   ([`PlaneParallel`](airmass::PlaneParallel),
//!   [`Young1994`](airmass::Young1994),
//!   [`Rozenberg1966`](airmass::Rozenberg1966),
//!   [`KrisciunasSchaefer1991`](airmass::KrisciunasSchaefer1991)).
//! - [`rayleigh`] — Rayleigh optical depth (Bodhaine et al. 1999) and
//!   the Rayleigh phase function.
//! - [`mie`] — Mie / aerosol optical depth via the Patat 2011 power-law,
//!   with site presets including [`MieParams::PARANAL`](mie::MieParams).
//! - [`airglow`] — Van Rhijn airglow emission-layer geometry.
//! - [`scattering`] — phase-function traits and Rayleigh / tabulated
//!   helpers; defines the dimensionless [`ScatteringFactor`] unit.
//! - [`extinction`] — Beer–Lambert transmission `exp(-X · τ)`.
//! - [`profile`] — [`AtmosphereProfile`]: site-bound aggregate of
//!   Rayleigh + Mie optical depth. Ozone is kept separate; apply
//!   `ozone::transmission_table()` independently.
//!
//! ## References
//!
//! - Bodhaine, B. A., Wood, N. B., Dutton, E. G., & Slusser, J. R. (1999).
//!   "On Rayleigh optical depth calculations". *J. Atmos. Oceanic Technol.*
//!   16, 1854.
//! - Patat, F., et al. (2011). "Optical atmospheric extinction over Cerro
//!   Paranal". *A&A* 527, A91.
//! - Krisciunas, K., & Schaefer, B. E. (1991). "A model of the brightness
//!   of moonlight". *PASP* 103, 1033.

pub mod airglow;
pub mod airmass;
pub mod extinction;
pub mod mie;
#[cfg(feature = "spectra")]
pub mod ozone;
pub mod profile;
pub mod rayleigh;
pub mod scattering;

pub use crate::qtty::{Transmittance, Transmittances};
pub use airglow::{van_rhijn_factor, van_rhijn_factor_with_radius};
pub use airmass::{
    airmass, AirmassFormula, DefaultAirmassFormula, Formula, KrisciunasSchaefer1991, PlaneParallel,
    Rozenberg1966, Young1994,
};
pub use extinction::transmission;
pub use mie::{mie_optical_depth, MieParams};
#[cfg(feature = "spectra")]
pub use ozone::transmission_table;
pub use profile::AtmosphereProfile;
pub use rayleigh::{rayleigh_optical_depth_bodhaine99, rayleigh_phase, DEFAULT_SCALE_HEIGHT};
#[cfg(feature = "tables")]
pub use scattering::TabulatedPhaseFunction;
pub use scattering::{PhaseFunction, RayleighPhaseFunction, ScatteringFactor};
