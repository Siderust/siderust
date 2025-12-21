//! Spherical coordinate types.
//!
//! This module re-exports types from [`affn::spherical`] and provides
//! convenient constructors and extensions for astronomical coordinate systems.

// Re-export all types from affn
pub use affn::spherical::*;

// Include astronomical extensions and type aliases
pub mod astro;
pub mod ext;

// Re-export for convenience at module level
pub use astro::{direction, position};
pub use ext::*;
