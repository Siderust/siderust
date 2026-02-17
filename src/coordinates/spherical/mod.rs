// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Spherical coordinate types with astronomical conventions.
//!
//! - [`Position<C, F, U>`]: spherical **position** (center + frame + distance)
//! - [`Direction<F>`]: spherical **direction** (frame-only, no center)
//!
//! The core spherical coordinate functionality lives in the `affn` crate
//! (domain-agnostic geometry kernel). This module provides:
//!
//! - Astronomical reference frame types (ICRS, EclipticMeanJ2000, Horizontal, etc.)
//! - Inherent named constructors and getters (`ra`/`dec`, `lon`/`lat`, `alt`/`az`)
//! - Convenient type aliases for common coordinate combinations
//!
//! # Usage
//!
//! Named constructors and getters are inherent — no trait imports needed:
//!
//! ```rust
//! use siderust::coordinates::spherical::direction;
//! use qtty::*;
//!
//! let dir = direction::ICRS::new(120.0 * DEG, 45.0 * DEG);
//! assert_eq!(dir.ra(), 120.0 * DEG);
//! assert_eq!(dir.dec(), 45.0 * DEG);
//! ```
//!
//! Or use generic types directly with `polar`/`azimuth`:
//!
//! ```rust
//! use siderust::coordinates::spherical::Direction;
//! use siderust::coordinates::frames::EquatorialMeanJ2000;
//! use qtty::*;
//!
//! let dir = Direction::<EquatorialMeanJ2000>::new(120.0 * DEG, 45.0 * DEG);
//! assert_eq!(dir.ra(), 120.0 * DEG);
//! ```

// Public alias modules (re-exported for public API)
pub mod direction;
pub mod position;

// Re-export affn types directly — no wrappers
pub use affn::spherical::Direction;
pub use affn::spherical::Position;
