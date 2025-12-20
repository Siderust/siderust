//! # Algebra Module
//!
//! This module contains the pure mathematical/geometric concepts for coordinate systems:
//! - Reference frames and centers (orientation and origin)
//! - Cartesian and spherical vector types (Position, Direction, Velocity)
//!
//! These are the abstract algebraic structures independent of physical/astronomical
//! context. For physical coordinate systems with convenience constructors and
//! astronomical conventions, see the `astro` module.

pub mod cartesian;
pub mod centers;
pub mod frames;
pub mod spherical;
