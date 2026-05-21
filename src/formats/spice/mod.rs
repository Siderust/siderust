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
//! - [`spk`] — SPK Type 2/3 segment reader built on top of [`daf`].
//! - [`kernel`] — High-level SPK kernel with segment chain resolution.
//! - [`segment`] — SPK segment evaluation (Type 2 and Type 3 Chebyshev).
//! - [`naif`] — NAIF body ID lookup table.
//!
//! ## Error handling
//!
//! Parse and kernel errors are reported as [`SpiceError`]. The low-level
//! `daf` and `spk` parsers use [`SpiceError::FormatParse`]; the kernel
//! layer uses the richer structured variants (`OutOfCoverage`, `NoChain`,
//! `UnsupportedDataType`, etc.).
//!
//! ## References
//!
//! - NAIF (2022). *DAF Required Reading*.
//!   <https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/daf.html>
//! - NAIF (2022). *SPK Required Reading*.
//!   <https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html>

pub mod daf;
pub mod error;
pub mod kernel;
pub mod naif;
pub mod segment;
pub mod spk;

pub use error::SpiceError;

pub use kernel::{LoadedSegment, SpkKernel};
pub use naif::{naif_id_for_name, well_known};
pub use segment::{segment_for_summary, ChebSegment, SpkSegment};
