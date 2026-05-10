// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Spacecraft dynamics
//!
//! Cartesian propagated state and (in later phases) force models, integrators,
//! and orbit-relative local frames for satellite mechanics.
//!
//! ## Scientific scope
//!
//! While [`crate::astro::orbit`] models orbits as Keplerian element sets,
//! many high-fidelity applications — orbit propagation under non-Keplerian
//! perturbations (J2, third-body, drag, SRP), state-transition matrices,
//! covariance transport, precise orbit determination — work directly with a
//! Cartesian inertial state `(r, v)`. This module provides the typed primitive
//! for that representation, [`OrbitState`], plus the surrounding spacecraft
//! description (mass, drag and SRP cross-sections).
//!
//! ## Technical scope
//!
//! All position, velocity, and acceleration values are [`affn`] vectors with
//! [`qtty`] units, so frame and unit constraints are checked at compile time.
//! Time derivatives are carried as the typed [`StateDerivative`] aggregate.
//!
//! Subsequent submodules (added in later phases) will host:
//!
//! - `forces` — `ForceModel` trait and standard models (two-body, J2, …),
//! - `integrators` — RK4 / DOPRI5,
//! - `covariance` — frame-tagged 6×6 state covariance with transport between
//!   inertial and local orbital frames.

pub mod covariance;
pub mod frames;
pub mod state;

pub use state::{
    OrbitState, Position, SpacecraftProperties, SpacecraftState, StateDerivative, Velocity,
};
