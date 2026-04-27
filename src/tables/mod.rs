// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Generic typed gridded tables.
//!
//! A second-generation, typed counterpart to [`crate::spectra`]. This module
//! provides:
//!
//! - [`Grid1D<X, V, S>`](Grid1D) — 1D linear interpolator over a typed `X`
//!   axis and a typed `V` value, with a configurable [`OutOfRange`] policy.
//! - [`Grid2D<X, Y, V, S>`](Grid2D) — 2D bilinear interpolator over two typed
//!   axes (row-major `[NY][NX]`), with per-axis [`OutOfRange`] policies.
//! - [`algo`] — untyped `f64` numerical kernels ([`algo::linear_1d`],
//!   [`algo::bilinear`]) for callers that need bit-for-bit
//!   `numpy.interp`-style parity with existing pipelines.
//! - [`Provenance`] — re-export of the shared dataset metadata record.
//!
//! Available behind the `tables` feature.

pub mod algo;
mod error;
mod grid1d;
mod grid2d;

pub use error::TableError;
pub use grid1d::Grid1D;
pub use grid2d::Grid2D;

pub use crate::interp::OutOfRange;
pub use crate::provenance::{DataSource, Provenance};
