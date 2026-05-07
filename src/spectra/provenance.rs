// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Provenance metadata for sampled spectra (re-export)
//!
//! ## Scientific scope
//!
//! Compatibility shim that re-exports the canonical provenance types
//! from [`crate::provenance`]. Provenance lives at the crate root so it
//! can be shared with [`crate::tables`]; this module preserves the
//! older `siderust::spectra::provenance::*` import path used by
//! downstream code.
//!
//! ## Technical scope
//!
//! Re-exports [`DataSource`](crate::provenance::DataSource) and
//! [`Provenance`](crate::provenance::Provenance).
//!
//! ## References
//!
//! - See [`crate::provenance`] for the canonical citations (IVOA
//!   Provenance Data Model; IERS Conventions).

pub use crate::provenance::{DataSource, Provenance};
