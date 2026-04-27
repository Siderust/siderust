// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Generic typed sampled spectra.
//!
//! This module provides a small, reusable building block for **1-D sampled
//! spectra** — a strictly-monotonic vector of x-values (typically wavelength)
//! paired with a same-length vector of y-values (irradiance, transmittance,
//! response, …) and an optional [`Provenance`] record.
//!
//! Both axes are typed via the `qtty` unit system, so dimensional analysis
//! flows through interpolation and integration naturally:
//!
//! ```ignore
//! use siderust::spectra::{SampledSpectrum, Interpolation, OutOfRange};
//! use siderust::qtty::{Nanometer, Nanometers};
//!
//! // x: nanometres, y: dimensionless transmittance
//! let xs = vec![Nanometers::new(300.0), Nanometers::new(400.0), Nanometers::new(500.0)];
//! let ys = vec![0.1_f64, 0.5, 0.9];
//! let s = SampledSpectrum::<Nanometer, _>::from_typed_xs_raw_ys(
//!     xs, ys, Interpolation::Linear, OutOfRange::ClampToEndpoints, None,
//! ).unwrap();
//! ```
//!
//! ## Design choices
//!
//! - **Interpolation policies are explicit.** The caller picks a variant of
//!   [`Interpolation`] at construction; only `Linear` is implemented today.
//! - **Out-of-range policy is explicit.** [`OutOfRange::ClampToEndpoints`] is
//!   the default to match `numpy.interp` and existing astronomy conventions;
//!   `Zero` and `Error` are the alternatives.
//! - **Integration returns typed quantities.** [`SampledSpectrum::integrate`]
//!   yields `Quantity<Prod<Y, X>>`, preserving units across operations.
//! - **Provenance is first-class.** Datasets pulled from publications, files,
//!   or external services should carry a [`Provenance`] record so downstream
//!   consumers can reason about reproducibility.
//!
//! ## Untyped algorithm helpers
//!
//! For consumers that already have raw `&[f64]` slices (e.g. legacy NSB code,
//! external file readers), the [`algo`] submodule exposes the same numerical
//! kernels operating on plain slices. These are guaranteed to be bit-for-bit
//! identical with the typed `SampledSpectrum` methods.

pub mod algo;
pub mod error;
pub mod interp;
pub mod integrate;
pub mod loaders;
pub mod provenance;
pub mod sampled;

pub use error::SpectrumError;
pub use interp::{Interpolation, OutOfRange};
pub use provenance::{DataSource, Provenance};
pub use sampled::SampledSpectrum;
