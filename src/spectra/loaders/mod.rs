// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Spectrum loaders
//!
//! ## Scientific scope
//!
//! Most published photometric filter curves, atmospheric transmission
//! tables, and instrument response functions are distributed as plain
//! ASCII files (and, more rarely, as FITS tables). Materialising those
//! files into typed [`crate::spectra::SampledSpectrum`] values is a
//! routine but error-prone step — wavelength columns may be in
//! Ångström, micrometre, or nanometre depending on the publisher;
//! header lines vary by source. This module collects the loaders that
//! perform that materialisation while preserving the typed-axis
//! invariants expected by the rest of the crate.
//!
//! ## Technical scope
//!
//! Currently provides the [`ascii`] submodule (whitespace- /
//! comma-separated two-column numeric files). Future additions
//! (FITS, VOTable) belong here as additional submodules.
//!
//! ## References
//!
//! None — pure I/O glue with no domain-specific algorithm.

pub mod ascii;
