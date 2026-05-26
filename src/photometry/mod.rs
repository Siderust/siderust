// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Astronomical photometry: photometric systems, passbands, and throughput.
//!
//! ## Scientific scope
//!
//! This module owns astronomical photometric concepts: spectral response curves
//! (passbands), throughput unit, synthetic photometry, and zero-point
//! conventions. The photometric systems currently bundled are Johnson–Cousins
//! UBVRI from Bessell (1990), the de-facto modern reference for synthetic
//! broad-band photometry.
//!
//! Generic sampled-spectrum infrastructure (interpolation, integration,
//! loaders) lives in [`optica::spectrum`].
//!
//! ## Technical scope
//!
//! `siderust::photometry` retains only astronomical data and domain APIs:
//! calibrated passband tables, provenance metadata, and the [`Throughput`]
//! unit type. All spectrum computation uses
//! [`optica::spectrum::SampledSpectrum`].

pub mod passbands;

pub use passbands::Throughput;
