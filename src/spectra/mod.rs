// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Typed sampled spectra
//!
//! ## Scientific scope
//!
//! Many photometric and atmospheric quantities are intrinsically functions
//! of wavelength — filter transmission curves, stellar irradiance, ozone
//! transmission, sky brightness — but in practice they are almost always
//! distributed and consumed as discretely sampled tables. This module
//! defines the typed container for that pattern: a strictly-monotonic
//! vector of wavelengths paired with a same-length vector of spectral
//! values, treated as the piecewise-linear interpolant of a continuous
//! spectral function over the optical / near-infrared regime.
//!
//! Validity is bounded by the sampled domain (extrapolation is governed
//! by the configured [`OutOfRange`] policy, not by an extrapolating
//! model) and by the linear-interpolation assumption between samples,
//! which is the de-facto convention in synthetic photometry and the
//! one used by Bessell (1990), the SVO Filter Profile Service, and
//! `numpy.interp`.
//!
//! ## Technical scope
//!
//! This module provides:
//!
//! - [`SampledSpectrum<X, Y, S>`] — typed container; `X` and `Y` are
//!   `qtty::Unit` markers (e.g. [`Nanometer`](crate::qtty::Nanometer) on
//!   the x-axis; transmittance, irradiance, throughput on the y-axis),
//!   `S` defaults to `f64`.
//! - [`Interpolation`] — interpolation policy (only `Linear` is
//!   implemented today).
//! - [`OutOfRange`] — out-of-range policy (clamp / zero / error),
//!   re-exported from [`crate::interp`].
//! - [`Provenance`] / [`DataSource`] — re-exported provenance metadata.
//! - [`SpectrumError`] — error taxonomy.
//! - Submodules: [`algo`] (untyped `&[f64]` kernels for bit-for-bit
//!   parity with reference pipelines), [`integrate`] (free-function
//!   re-exports), [`loaders`] (file-format ingestion), [`passbands`]
//!   (curated filter datasets).
//!
//! Integration over a typed spectrum returns a typed
//! `Quantity<Prod<Y, X>>`; dimensional analysis flows through
//! interpolation and integration unchanged.
//!
//! ```ignore
//! use siderust::spectra::{SampledSpectrum, Interpolation, OutOfRange};
//! use siderust::qtty::{Nanometer, Nanometers};
//!
//! let xs = vec![Nanometers::new(300.0), Nanometers::new(400.0), Nanometers::new(500.0)];
//! let ys = vec![0.1_f64, 0.5, 0.9];
//! let s = SampledSpectrum::<Nanometer, _>::from_typed_xs_raw_ys(
//!     xs, ys, Interpolation::Linear, OutOfRange::ClampToEndpoints, None,
//! ).unwrap();
//! ```
//!
//! ## References
//!
//! - Bessell, M. S. (1990). "UBVRI Passbands". *Publications of the
//!   Astronomical Society of the Pacific* **102**, 1181.
//!   doi:10.1086/132749.
//! - Press, W. H., Teukolsky, S. A., Vetterling, W. T., Flannery, B. P.
//!   (1992). *Numerical Recipes in C*, 2nd ed., §3.1 (linear
//!   interpolation). Cambridge University Press.
//! - NumPy developers. *numpy.interp* documentation
//!   (boundary-handling semantics).

pub mod algo;
pub mod error;
pub mod integrate;
pub mod interp;
pub mod loaders;
pub mod passbands;
pub mod provenance;
pub mod sampled;

pub use error::SpectrumError;
pub use interp::{Interpolation, OutOfRange};
pub use provenance::{DataSource, Provenance};
pub use sampled::SampledSpectrum;
