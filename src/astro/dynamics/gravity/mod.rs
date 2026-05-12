// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Geopotential gravity models for astrodynamics.
//!
//! ## Quick reference
//!
//! | Type | Purpose |
//! |------|---------|
//! | [`GravityFieldProvider`] | Trait for normalised Stokes-coefficient sources |
//! | [`GravityConstants`]     | Legacy helper: GM, R, and max_degree |
//! | [`TwoBodyEarth`]         | Point-mass Earth (degree 0) |
//! | [`LowDegreeEarth`]       | EGM2008 Earth through degree/order 4 |
//! | `geopotential_acceleration` | Spherical-harmonic acceleration kernel |

pub mod acceleration;
pub mod egm_low_degree;
pub mod provider;

pub use acceleration::geopotential_acceleration;
pub use egm_low_degree::{LowDegreeEarth, TwoBodyEarth};
pub use provider::{GravityConstants, GravityFieldProvider};
