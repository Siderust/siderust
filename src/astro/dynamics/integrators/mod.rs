// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Numerical integrators for [`OrbitState`].
//!
//! All integrators consume a [`ForceModel`] and integrate the typed
//! [`StateDerivative`] returned by [`StateDerivative::new`]; no raw
//! `[f64; N]` plumbing escapes the public API.

pub mod dopri5;
pub mod rk4;

pub use dopri5::{dopri5_propagate, dopri5_step, Tolerance};
pub use rk4::{rk4_propagate, rk4_propagate_series, rk4_step};

#[allow(unused_imports)]
use super::forces::ForceModel;
#[allow(unused_imports)]
use crate::astro::dynamics::{state::StateDerivative, OrbitState};
