// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Observation Module
//!
//! This module provides types and transformations for observational astronomy,
//! separating **geometric** coordinates from **observed** (apparent) coordinates.
//!
//! ## Core Concepts
//!
//! ### Observational State
//!
//! Celestial directions can be in different observational states:
//!
//! - **Astrometric** (`Astrometric`): The geometric direction to a target, corrected for
//!   proper motion and parallax, but *without* corrections for observer motion effects.
//!   This is the "true" direction to the object in space.
//!
//! - **Apparent** (`Apparent`): The direction as it would be observed, including corrections
//!   for aberration (and potentially other effects like deflection). This is where the
//!   object *appears* to be due to the finite speed of light and observer motion.
//!
//! ### Separation of Concerns
//!
//! The observation module maintains a clear separation:
//!
//! - **Center transforms** (in `transform/centers`): Pure geometric translations.
//!   Moving the origin from Sun to Earth is just vector subtraction. No aberration.
//!
//! - **Frame transforms** (in `transform/frames`): Pure geometric rotations.
//!   Changing from ecliptic to equatorial is just a rotation matrix. No aberration.
//!
//! - **Observation transforms** (this module): Observer-dependent effects.
//!   Aberration depends on observer velocity. Requires explicit `ObserverState`.
//!
//! ## Usage
//!
//! ```rust
//! use siderust::coordinates::observation::{Astrometric, Apparent, ObserverState};
//! use siderust::coordinates::spherical::direction::EquatorialMeanJ2000;
//! use siderust::astro::JulianDate;
//! use qtty::*;
//!
//! // A geometric direction to a star
//! let astrometric_dir: Astrometric<EquatorialMeanJ2000> = Astrometric::new(
//!     EquatorialMeanJ2000::new(45.0 * DEG, 20.0 * DEG)
//! );
//!
//! // Get observer state (Earth-bound observer at J2000)
//! let obs = ObserverState::geocentric(JulianDate::J2000);
//!
//! // Convert to apparent direction (applies aberration)
//! let apparent_dir: Apparent<EquatorialMeanJ2000> = astrometric_dir.to_apparent(&obs);
//!
//! // The transformation is explicit - you can't confuse the two types
//! ```
//!
//! ## Design Rationale
//!
//! This design follows the IAU-style coordinate pipeline:
//!
//! 1. Catalog positions (ICRS, typically astrometric)
//! 2. Apply proper motion, parallax → geometric position at date
//! 3. Center transform (barycentric → geocentric) → geometric position from Earth
//! 4. Frame transform (ICRS → equatorial of date) → geometric direction in observing frame
//! 5. Apply aberration → **apparent** direction
//! 6. (Optional) Apply refraction → observed direction
//!
//! Each step is explicit. The type system tracks observational state.

mod observational_direction;
mod observer_state;

pub use observational_direction::{Apparent, Astrometric};
pub use observer_state::ObserverState;

// Re-export aberration utilities for convenience
pub use crate::astro::aberration::{
    apply_aberration_to_direction, remove_aberration_from_direction,
};
