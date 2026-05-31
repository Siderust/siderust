// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # NAIF SPICE file-format parsers
//!
//! Low-level parsers for NAIF SPICE text and binary kernels, including SPK
//! ephemerides, CK attitude kernels, frame definitions, leapseconds, and
//! related spacecraft-geometry metadata.
//!
//! ## Modules
//!
//! - [`daf`] — DAF container parsing, including generic [`daf::DafRaw`].
//! - [`spk`] — SPK Type 2/3 segment reader built on top of [`daf`].
//! - [`kernel`] — High-level SPK kernel with segment chain resolution.
//! - [`segment`] — SPK segment evaluation (Type 2/3/9/13).
//! - [`text`] — Generic SPICE text-kernel parser.
//! - [`lsk`] — Leapseconds-kernel parsing.
//! - [`fk`] — Frame-kernel parsing.
//! - [`pck`] — Text-PCK body-orientation parsing.
//! - [`ck`] — Binary CK attitude-kernel parsing.
//! - [`sclk`] — Spacecraft-clock kernel parsing.
//! - [`ik`] — Instrument-kernel parsing.
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

pub mod ck;
pub mod daf;
pub mod error;
pub mod fk;
pub mod ik;
pub mod kernel;
pub mod lsk;
pub mod naif;
pub mod pck;
pub mod sclk;
pub mod segment;
pub mod spk;
pub mod text;

pub use ck::{CkKernel, CkRecord, CkSegment1};
pub use error::SpiceError;
pub use fk::{FrameClass, FrameKernel, FrameSpec, TkSpec};
pub use ik::IkKernel;
pub use kernel::{LoadedSegment, SpkKernel};
pub use lsk::LeapSecondKernel;
pub use naif::{naif_id_for_name, well_known};
pub use pck::{BodyOrientation, PckKernel};
pub use sclk::SclkKernel;
pub use segment::{segment_for_summary, ChebSegment, SpkSegment};
pub use text::TextKernel;
