// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! # POD propagation and STM
//!
//! ## Scientific scope
//! Integrator adapters, variational propagation, and STM validation for POD.

pub mod adapter;
pub mod error;
pub mod integrators;
pub mod pod_error;
pub mod validation;
pub mod variational;

pub use adapter::{propagate_orbit, propagate_spacecraft};
pub use error::DynamicsError;
pub use integrators::{Dop853Integrator, Dopri5Integrator, Integrator, Rk4Integrator};
pub use pod_error::PodDynamicsError;
#[doc(no_inline)]
pub use principia::{propagate_stm, StateTransitionMatrix};
pub use variational::{
    param_partials_central_diff, ParamColumn, ParamStmReport, PropagatedArc, VariationalPropagator,
};
