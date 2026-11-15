// SPDX-License-Identifier: AGPL-3.0-only
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
//! ## Geometric and aberration-corrected horizon directions
//!
//! ```rust
//! use siderust::coordinates::{
//!     centers::Geodetic,
//!     frames::{EquatorialMeanJ2000, EquatorialTrueOfDate, ECEF},
//!     observation::{Astrometric, ObserverState},
//!     spherical::Direction,
//!     transform::{AstroContext, SphericalDirectionAstroExt, horizontal::project_direction_with},
//! };
//! use siderust::qtty::{DEG, M};
//!
//! let jd = siderust::J2000;
//! let context: AstroContext = AstroContext::default();
//! let site = Geodetic::<ECEF>::new(-17.9 * DEG, 28.8 * DEG, 2300.0 * M);
//! let geometric = Direction::<EquatorialMeanJ2000>::new(45.0 * DEG, 20.0 * DEG);
//! let observer = ObserverState::topocentric(&site, jd);
//! let apparent = Astrometric::new(geometric).to_apparent(&observer);
//!
//! let geometric_of_date = geometric.to_frame_with::<EquatorialTrueOfDate, _>(&jd, &context);
//! let apparent_of_date = apparent.direction().to_frame_with::<EquatorialTrueOfDate, _>(&jd, &context);
//! let geometric_horizon = project_direction_with(&geometric_of_date, &site, jd, &context);
//! let apparent_horizon = project_direction_with(&apparent_of_date, &site, jd, &context);
//! assert_ne!(geometric_horizon.az(), apparent_horizon.az());
//! ```
//!
//! This direction-only example applies annual and diurnal aberration explicitly.
//! For finite distances, form the observer-relative line of sight before
//! applying aberration. The observer state and frame rotation use the same TT
//! epoch and default models.
//!
//! ## Design Rationale
//!
//! Carry the target's epoch and frame through geometric transformations.
//! Aberration uses an observer velocity at the reception epoch, expressed in
//! the same axes as the sightline. The wrappers perform that frame alignment
//! for their supported inertial frames. Light-time, deflection, and atmospheric
//! refraction are separate operations.

mod observational_direction;
mod observer_state;

pub use observational_direction::{Apparent, Astrometric};
pub use observer_state::ObserverState;

// Re-export aberration utilities for convenience
pub use crate::astro::aberration::{
    apply_aberration_to_direction, remove_aberration_from_direction,
};
