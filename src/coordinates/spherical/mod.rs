// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Spherical coordinate types with astronomical conventions.
//!
//! - [`Position<C, F, U>`]: spherical **position** (center + frame + distance)
//! - [`Direction<F>`]: spherical **direction** (frame-only, no center)
//!
//! This module provides siderust-owned wrapper types around `affn` with
//! frame-specific inherent constructors following IAU conventions.
//!
//! # Architecture Note
//!
//! **Core spherical coordinate functionality is in the `affn` crate** (domain-agnostic
//! geometry kernel). This module provides thin wrappers with:
//! - Astronomical-specific frame types (ICRS, Ecliptic, Horizontal, etc.)
//! - Frame-specific constructors with IAU naming (e.g., `ra`/`dec`, `lon`/`lat`, `az`/`alt`)
//! - Getter methods with astronomical field names (instead of generic `polar`/`azimuth`)
//!
//! # Usage
//!
//! No trait imports needed—just use the types directly:
//!
//! ```rust
//! use siderust::coordinates::spherical::direction;
//! use qtty::*;
//!
//! // Inherent constructor with astronomical argument order
//! let dir = direction::ICRS::new(120.0 * DEG, 45.0 * DEG);  // (ra, dec)
//! assert_eq!(dir.ra(), 120.0 * DEG);
//! assert_eq!(dir.dec(), 45.0 * DEG);
//! ```
//!
//! Or use generic types directly:
//!
//! ```rust
//! use siderust::coordinates::spherical::Direction;
//! use siderust::coordinates::frames::EquatorialMeanJ2000;
//! use qtty::*;
//!
//! let dir = Direction::<EquatorialMeanJ2000>::new(120.0 * DEG, 45.0 * DEG);
//! ```

use qtty::Degrees;

// Internal helper functions
#[inline]
fn normalize_azimuth(az: Degrees) -> Degrees {
    az.normalize()
}

#[inline]
fn clamp_polar(polar: Degrees) -> Degrees {
    polar.wrap_quarter_fold()
}

// Core type modules
mod direction_core;
mod direction_impls;
mod position_core;
mod position_impls;

#[cfg(feature = "serde")]
mod feature_serde;

// Public alias modules (re-exported for public API)
pub mod direction;
pub mod position;

// Public exports
pub use direction_core::Direction;
pub use position_core::Position;
