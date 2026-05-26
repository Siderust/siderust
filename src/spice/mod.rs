// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! High-level SPICE kernel context for astronomical and spacecraft geometry.
//!
//! ## Scientific scope
//!
//! This module provides the user-facing SPICE context for ephemeris state
//! queries, frame registration, limited frame rotation composition, and basic
//! time-offset support backed by SPICE kernels.
//!
//! ## Technical scope
//!
//! V1 layers `formats::spice` parsers into a loadable kernel set with NAIF-like
//! priority rules. It supports SPK state lookups plus selected FK, PCK text,
//! LSK, SCLK, CK, and IK metadata use cases.
//!
//! ## References
//!
//! - NAIF SPICE Overview: <https://naif.jpl.nasa.gov/naif/aboutspice.html>
//! - NAIF. *Kernel Required Reading*.

#![forbid(unsafe_code)]
#![deny(missing_docs)]

mod context;
mod error;
mod frame_registry;
mod kernel_set;

pub use context::SpiceContext;
pub use error::SpiceContextError;
pub use frame_registry::FrameRegistry;
pub use kernel_set::KernelSet;
