//! Interpolation, typed data tables, and data provenance primitives.
//!
//! Re-exported from [`siderust`] as `siderust::interp`, `siderust::tables`, and
//! `siderust::provenance`. Downstream crates that want only the core data-layer
//! primitives can depend on this crate directly.

pub mod interp;
pub mod provenance;

#[cfg(feature = "tables")]
pub mod tables;
