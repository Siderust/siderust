//! # Algebra Module
//!
//! This module contains the pure mathematical/geometric concepts for coordinate systems:
//! - Reference frames and centers (orientation and origin)
//! - Cartesian and spherical vector types (Position, Direction, Velocity)
//!
//! The coordinate types are re-exported from the `affn` crate (the pure geometry kernel).
//! Astronomy-specific frames and centers are defined locally in this module.

pub mod centers;
pub mod frames;

// Re-export cartesian and spherical types from affn
// This provides zero-cost access to affn's coordinate algebra
pub mod cartesian {
    //! Cartesian coordinate types.
    //!
    //! This module re-exports types from [`affn::cartesian`].
    pub use affn::cartesian::*;
}

pub mod spherical {
    //! Spherical coordinate types.
    //!
    //! This module re-exports types from [`affn::spherical`].
    pub use affn::spherical::*;
}
