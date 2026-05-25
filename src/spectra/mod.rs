// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Astronomical spectral data: photometric passbands and throughput.
//!
//! ## Scientific scope
//!
//! This module provides photometric passbands (wavelength-dependent
//! throughput curves) used for synthetic photometry and sky-brightness
//! calculations. The generic sampled-spectrum infrastructure (interpolation,
//! integration, loaders) lives in [`optica::spectrum`].
//!
//! ## Technical scope
//!
//! `siderust::spectra` retains only astronomical data: calibrated passband
//! tables, provenance metadata, and the [`Throughput`] unit type. All
//! spectrum computation uses [`optica::spectrum::SampledSpectrum`].

pub mod passbands;

pub use passbands::Throughput;
