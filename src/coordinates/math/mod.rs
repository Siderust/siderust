//! # Mathematical Kernels for Coordinate Operations
//!
//! This module provides pure mathematical functions separated from the coordinate
//! type system. These functions operate on primitive `f64` values and are designed
//! to be:
//!
//! - **Pure**: No side effects, deterministic outputs for given inputs
//! - **Testable**: Easy to unit test in isolation with primitive values
//! - **Reusable**: Can be used by FFI, benchmarks, or alternative APIs
//! - **Type-agnostic**: No phantom type parameters or unit wrappers
//!
//! ## Modules
//!
//! - [`rotations`]: Frame rotation matrices (ecliptic ↔ equatorial, etc.)
//! - [`conversions`]: Coordinate system conversions (spherical ↔ cartesian)
//! - [`geometry`]: Geometric operations (distance, normalization, angular separation)
//!
//! ## Design Philosophy
//!
//! The coordinate types in [`super::cartesian`] and [`super::spherical`] provide
//! the type-safe API with phantom types for compile-time safety. This module
//! contains the underlying arithmetic that those types delegate to.
//!
//! This separation allows:
//! 1. Clear distinction between API surface and implementation
//! 2. Easier testing of mathematical correctness
//! 3. Potential optimization of hot paths
//! 4. Reuse in contexts where type parameters add overhead

pub mod conversions;
pub mod geometry;
pub mod rotations;
