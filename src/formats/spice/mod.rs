// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # NAIF SPICE file-format parsers
//!
//! Low-level parsers for the NAIF SPICE binary formats used by JPL
//! planetary ephemeris files (BSP kernels).
//!
//! ## Modules
//!
//! - [`daf`] — DAF (Double Precision Array File) container parser.
//! - [`spk`] — SPK Type 2 segment reader built on top of [`daf`].
//!
//! ## Error handling
//!
//! Parse errors are reported as [`SpiceError`], which is independent of
//! the dataset catalog error type. Callers that bridge these layers
//! can convert via `impl From<SpiceError> for datasets::DatasetError`.
//!
//! ## References
//!
//! - NAIF (2022). *DAF Required Reading*.
//!   <https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/daf.html>
//! - NAIF (2022). *SPK Required Reading*.
//!   <https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html>

pub mod daf;
pub mod spk;

/// Errors produced by SPICE format parsers.
#[derive(Debug)]
pub enum SpiceError {
    /// Structural parse failure (bad magic, truncated data, unsupported type, …).
    Parse(String),
}

impl std::fmt::Display for SpiceError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SpiceError::Parse(msg) => write!(f, "SPICE parse error: {}", msg),
        }
    }
}

impl std::error::Error for SpiceError {}
