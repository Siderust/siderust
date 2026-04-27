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
//! - [`Grid3D<X, Y, Z, V, S>`](Grid3D) — 3D trilinear interpolator over three
//!   typed axes (storage convention `(iz·NY + iy)·NX + ix`), with per-axis
//!   [`OutOfRange`] policies.
//! - [`algo`] — untyped `f64` numerical kernels ([`algo::linear_1d`],
//!   [`algo::bilinear`], [`algo::trilinear`]) for callers that need bit-for-bit
//!   `numpy.interp`-style parity with existing pipelines.
//! - [`Provenance`] — re-export of the shared dataset metadata record.
//!
//! Available behind the `tables` feature.

pub mod algo;
mod error;
mod grid1d;
mod grid2d;
mod grid3d;

pub use error::TableError;
pub use grid1d::Grid1D;
pub use grid2d::Grid2D;
pub use grid3d::Grid3D;

pub use crate::interp::OutOfRange;
pub use crate::provenance::{DataSource, Provenance};

/// Direction of a strictly-monotonic axis.
///
/// Returned by [`algo::validate_axis`] and stored inside [`Grid1D`] /
/// [`Grid2D`] / [`Grid3D`] so that [`OutOfRange::ClampToEndpoints`] and the
/// internal [`algo::locate`] kernel know which end is the low-value end and
/// which is the high-value end.
///
/// For `ClampToEndpoints` on a **descending** axis, a query above `xs[0]`
/// (the largest stored value) clamps to `xs[0]` / `ys[0]`, while a query
/// below `xs[n-1]` (the smallest stored value) clamps to `xs[n-1]` /
/// `ys[n-1]`.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum AxisDirection {
    /// Axis values strictly increase with index (`xs[i+1] > xs[i]`).
    Ascending,
    /// Axis values strictly decrease with index (`xs[i+1] < xs[i]`).
    Descending,
}
