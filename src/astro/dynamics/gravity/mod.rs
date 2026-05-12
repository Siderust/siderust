// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Geopotential gravity models and acceleration kernel.
//!
//! ## Scope
//!
//! Provides the [`GravityFieldProvider`] trait abstraction, built-in providers
//! ([`TwoBodyEarth`], [`LowDegreeEarth`]), and the low-level spherical-harmonic
//! acceleration kernel ([`geopotential_acceleration`]).
//!
//! ## Equations
//!
//! The gravitational acceleration is derived from the potential
//!
//! ```text
//! U = (GM/r) Σ_{n=0}^{N} (R/r)^n  Σ_{m=0}^{n}
//!       [ C̄_{nm} P̄_{nm}(sin φ) cos(mλ)
//!       + S̄_{nm} P̄_{nm}(sin φ) sin(mλ) ]
//! a = ∇U
//! ```
//!
//! where C̄, S̄ are fully-normalised Stokes coefficients.
//!
//! ## Units & frames
//!
//! Position km (geocentric).  Acceleration km/s².
//! No frame restrictions (the gradient is frame-independent).
//!
//! ## Quick reference
//!
//! | Type | Purpose |
//! |------|---------|
//! | [`GravityFieldProvider`] | Trait for normalised Stokes-coefficient sources |
//! | [`GravityConstants`]     | Legacy helper: GM, R, and max_degree |
//! | [`TwoBodyEarth`]         | Point-mass Earth (degree 0) |
//! | [`LowDegreeEarth`]       | EGM2008 Earth through degree/order 4 |
//! | [`geopotential_acceleration`] | Spherical-harmonic acceleration kernel |
//!
//! ## References
//!
//! * Montenbruck & Gill, *Satellite Orbits*, §3.2.
//! * Vallado, *Fundamentals of Astrodynamics and Applications*, §8.6.
//! * IERS Conventions (2010).

pub mod acceleration;
pub mod egm_low_degree;
pub mod provider;

pub use acceleration::geopotential_acceleration;
pub use egm_low_degree::{LowDegreeEarth, TwoBodyEarth};
pub use provider::{GravityConstants, GravityFieldProvider};
