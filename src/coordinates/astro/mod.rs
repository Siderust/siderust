//! # Astro Module
//!
//! This module contains the physical/astronomical coordinate system implementations
//! with convenience constructors and astronomical conventions (e.g., RA/Dec, lon/lat naming).
//!
//! The pure algebraic coordinate types are in the `algebra` module. This module
//! provides the convenient, domain-specific wrappers and extensions for astronomical use.

pub mod cartesian;
pub mod spherical;
