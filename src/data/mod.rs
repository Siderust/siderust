// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Centralised Siderust data management
//!
//! ## Scientific scope
//!
//! This module is the single public entry point for dataset cataloging,
//! runtime acquisition metadata, provenance, checksums, and the archive
//! registry used by Siderust astronomy, mission geometry, and POD workflows.
//!
//! ## Technical scope
//!
//! Public users interact with [`DatasetId`], [`DatasetMeta`], [`DATASETS`],
//! [`lookup`], [`DataSource`], [`Provenance`], [`checksum`], and, when the
//! `runtime-data` feature is enabled, [`runtime`]. JPL planetary ephemeris
//! binary data is crate-private under `data::compiled::jpl`.
//!
//! ## References
//!
//! - IAU SOFA standards documents.
//! - JPL DE technical memoranda.
//! - NASA NAIF SPICE Toolkit documentation.

pub mod archive;
mod catalog;
pub mod checksum;
pub(crate) mod compiled;
mod provenance;
#[cfg(feature = "runtime-data")]
pub mod runtime;

pub use catalog::*;
pub use provenance::{DataSource, Provenance};

/// Re-export of the embedded Sun-Earth Lagrange Chebyshev archive API.
#[cfg(feature = "lagrange-centers")]
pub mod embedded_lagrange {
    pub use crate::ephemeris::lagrange::*;
}
