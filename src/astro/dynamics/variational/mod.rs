// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Analytic state-transition-matrix propagation via variational equations.
//!
//! This module provides a self-contained implementation of the variational
//! equations of motion, integrating the orbit ODE and the STM equation
//! simultaneously with a classical fixed-step RK4 scheme.
//!
//! ## Quick start
//!
//! ```rust
//! use siderust::astro::dynamics::variational::{propagate_stm, VariationalConfig};
//! use siderust::astro::dynamics::forces::TwoBody;
//! use siderust::astro::dynamics::context::DynamicsContext;
//! use siderust::astro::dynamics::{OrbitState, Position, Velocity};
//! use siderust::coordinates::frames::GCRS;
//! use siderust::time::JulianDate;
//! use siderust::qtty::Second;
//!
//! let s0 = OrbitState::new_at_jd(
//!     JulianDate::from_raw_unchecked(qtty::Day::new(2_451_545.0)),
//!     Position::<GCRS>::new(7_000.0, 0.0, 0.0),
//!     Velocity::<GCRS>::new(0.0, 7.5, 0.0),
//! );
//! let ctx = DynamicsContext::empty();
//! let force = TwoBody::earth();
//! let dt = Second::new(600.0);
//!
//! let (s_final, phi) = propagate_stm(&force, s0, dt, &ctx).unwrap();
//! // phi is the 6×6 state-transition matrix Φ(t, t₀) tagged to GCRS.
//! ```
//!
//! ## See also
//!
//! - [`equations`] — building the A matrix and the variational derivative.
//! - [`propagator`] — the propagator entry points and `VariationalConfig`.
//! - The finite-difference STM in `super::stm` is preserved side-by-side.

pub mod equations;
pub mod propagator;

pub use propagator::{propagate_stm, propagate_stm_with, StateTransitionMatrix, VariationalConfig};
