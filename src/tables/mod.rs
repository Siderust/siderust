// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Typed multi-dimensional look-up tables
//!
//! ## Scientific scope
//!
//! Many astronomical models are distributed as gridded look-up tables
//! over two or three independent parameters: the Leinert et al. (1998)
//! zodiacal-light intensity as a function of ecliptic latitude and
//! solar elongation, the ESO sky-brightness model as a function of
//! wavelength × airmass × moon separation, refraction and atmospheric
//! transmission tables, and so on. This module extends the 1-D
//! sampled-spectrum concept of [`crate::spectra`] to those higher-rank
//! grids while keeping the same hygiene properties: typed axes via
//! `qtty::Unit` markers, configurable out-of-range policy, and an
//! optional [`Provenance`] record for reproducibility.
//!
//! Validity is bounded by the sampled domain on each axis; the
//! interpolation kernels are piecewise-multilinear and match
//! `numpy.interp`-style boundary semantics so that results agree
//! bit-for-bit with reference pipelines.
//!
//! ## Technical scope
//!
//! Available behind the `tables` feature. Provides:
//!
//! - [`Grid1D<X, V, S>`](Grid1D) — 1-D linear interpolator over a typed
//!   `X` axis with a typed `V` value and a configurable [`OutOfRange`]
//!   policy.
//! - [`Grid2D<X, Y, V, S>`](Grid2D) — 2-D bilinear interpolator over
//!   two typed axes (row-major `[NY][NX]`) with per-axis [`OutOfRange`]
//!   policies; [`ConstantRegion`] short-circuits a half-open
//!   rectangular sub-region to a constant value (the lower-corner
//!   clamp pattern from NSB Leinert lookups).
//! - [`Grid3D<X, Y, Z, V, S>`](Grid3D) — 3-D trilinear interpolator
//!   with storage convention `(iz·NY + iy)·NX + ix`.
//! - [`algo`] — untyped `f64` numerical kernels ([`algo::linear_1d`],
//!   [`algo::bilinear`], [`algo::trilinear`], [`algo::validate_axis`],
//!   [`algo::locate`]) for callers that need numpy-parity with existing
//!   pipelines.
//! - [`AxisDirection`] — direction tag returned by
//!   [`algo::validate_axis`] and stored inside each grid for descending
//!   axes.
//! - [`Provenance`] / [`DataSource`] — re-exports of the shared
//!   provenance metadata.
//! - [`OutOfRange`] — re-export of the shared boundary policy enum.
//! - [`TableError`] — error taxonomy.
//!
//! ## References
//!
//! - Leinert, Ch., Bowyer, S., Haikala, L. K., et al. (1998). "The
//!   1997 reference of diffuse night sky brightness". *Astronomy &
//!   Astrophysics Supplement Series* **127**, 1–99.
//!   doi:10.1051/aas:1998105.
//! - Noll, S., Kausch, W., Barden, M., et al. (2012). "An atmospheric
//!   radiation model for Cerro Paranal. I. The optical spectral
//!   range". *Astronomy & Astrophysics* **543**, A92.
//!   doi:10.1051/0004-6361/201219040.
//! - European Southern Observatory. *SkyCalc — ESO Sky Model
//!   Calculator*. <https://www.eso.org/observing/etc/skycalc/skycalc.htm>.

pub mod algo;
mod error;
mod grid1d;
mod grid2d;
mod grid3d;

pub use error::TableError;
pub use grid1d::Grid1D;
pub use grid2d::{ConstantRegion, Grid2D};
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
