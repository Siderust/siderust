// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Numerical integrators for [`OrbitState`].
//!
//! ## Scope
//!
//! Provides fixed-step (RK4) and adaptive-step (DOPRI5, DOP853) integrators
//! for the Cartesian state derivative `dy/dt = [v; a(t, y)]`.
//!
//! ## Algorithm overview
//!
//! All integrators consume a [`ForceModel`] and integrate the typed
//! [`StateDerivative`] returned by [`StateDerivative::new`]; no raw
//! `[f64; N]` plumbing escapes the public API.
//!
//! | Integrator | Type | Order | Stages | Use case |
//! |------------|------|-------|--------|----------|
//! | RK4 | Fixed | 4th | 4 | Fast, simple propagation; requires careful step choice |
//! | DOPRI5 | Adaptive | 5th | 7 | General-purpose; sufficient for most LEO applications |
//! | DOP853 | Adaptive | 8th | 12 | High-precision; POD, maneuver design |
//!
//! ## Units & frames
//!
//! Position km, velocity km/s, acceleration km/s² (GCRS by default).
//! Time steps in seconds.
//!
//! ## Trait abstraction
//!
//! Two traits define the integrator contract:
//! - [`FixedStepper`] — fixed-step integrators (RK4)
//! - [`AdaptiveStepper`] — adaptive-step integrators (DOPRI5, DOP853)
//!
//! Both are generic over center `C` and frame `F`.
//!
//! ## References
//!
//! * Hairer, Norsett & Wanner, *Solving Ordinary Differential Equations I*,
//!   2nd ed., Springer (1993).
//! * Vallado, *Fundamentals of Astrodynamics and Applications* (2013), §4.4.
//! * Montenbruck & Gill, *Satellite Orbits* (2001), §4.4.

pub mod dop853;
pub mod dopri5;
pub mod rk4;

pub use dop853::{dop853_propagate, dop853_step, Dop853, Dop853Step};
pub use dopri5::{dopri5_propagate, dopri5_step, Dopri5};
pub use rk4::{rk4_propagate, rk4_propagate_series, rk4_step, Rk4};

#[allow(unused_imports)]
use super::forces::ForceModel;
#[allow(unused_imports)]
use crate::astro::dynamics::{state::StateDerivative, OrbitState};

use crate::astro::dynamics::context::DynamicsContext;
use crate::astro::dynamics::errors::DynamicsError;
use crate::coordinates::centers::{Geocentric, ReferenceCenter};
use crate::coordinates::frames::{ReferenceFrame, GCRS};
use crate::qtty::Second;

/// Contract for variable-step Runge-Kutta integrators used by the generic
/// [`crate::astro::dynamics::propagation::propagate`] driver.
///
/// `C` is the reference center, `F` the reference frame.  The method
/// [`AdaptiveStepper::step`] attempts a single step of size `h_try` and
/// returns `(accepted_state, h_used, h_next, steps_rejected)` where
/// `steps_rejected` is the number of internal trial steps the controller
/// discarded before accepting this step.
pub trait AdaptiveStepper<C = Geocentric, F = GCRS>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
{
    fn step<FM: ForceModel<C, F>>(
        &self,
        force: &FM,
        state: &OrbitState<C, F>,
        h_try: Second,
        ctx: &DynamicsContext,
    ) -> Result<(OrbitState<C, F>, Second, Second, u32), DynamicsError>;
}

/// Contract for fixed-step integrators (e.g. classical RK4).
pub trait FixedStepper<C = Geocentric, F = GCRS>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
{
    fn step<FM: ForceModel<C, F>>(
        &self,
        force: &FM,
        state: &OrbitState<C, F>,
        h: Second,
        ctx: &DynamicsContext,
    ) -> Result<OrbitState<C, F>, DynamicsError>;
}

#[inline]
pub(super) fn state_component(s: &OrbitState, i: usize) -> f64 {
    match i {
        0 => s.position.x().value(),
        1 => s.position.y().value(),
        2 => s.position.z().value(),
        3 => s.velocity.x().value(),
        4 => s.velocity.y().value(),
        5 => s.velocity.z().value(),
        _ => panic!("index out of range"),
    }
}

#[inline]
pub(super) fn deriv_component(d: &StateDerivative, i: usize) -> f64 {
    match i {
        0 => d.vel.x().value(),
        1 => d.vel.y().value(),
        2 => d.vel.z().value(),
        3 => d.acc.x().value(),
        4 => d.acc.y().value(),
        5 => d.acc.z().value(),
        _ => panic!("index out of range"),
    }
}

#[inline]
pub(super) fn rhs<FM: ForceModel>(
    force: &FM,
    s: &OrbitState,
    ctx: &DynamicsContext,
) -> Result<StateDerivative, DynamicsError> {
    Ok(StateDerivative::new(
        s.velocity,
        force.acceleration(s, ctx)?,
    ))
}

#[inline]
pub(super) fn state_at(s: &OrbitState, d: &StateDerivative, h: f64, dt: f64) -> OrbitState {
    let new_epoch = s.epoch + Second::new(dt);
    let advanced = s.advance(d, Second::new(h));
    OrbitState {
        epoch: new_epoch,
        position: advanced.position,
        velocity: advanced.velocity,
    }
}
