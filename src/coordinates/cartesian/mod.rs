//! Cartesian coordinate types.
//!
//! This module re-exports types from [`affn::cartesian`] and provides
//! convenient type aliases for astronomical coordinate systems.

// Re-export all types from affn
pub use affn::cartesian::*;

// Include astronomical type aliases
pub mod astro;

// Re-export specific type aliases at this level for convenience
pub use astro::{direction, position, velocity};
