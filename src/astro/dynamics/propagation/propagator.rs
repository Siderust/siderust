// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! High-level propagator facade for astronomy callers.
//!
//! ## Scientific scope
//!
//! Wraps `principia`'s numerical integrators in an astronomy-typed facade:
//! accepts and returns [`OrbitState`] (TT, geocentric, GCRS) so callers never
//! need to touch raw `principia` types for standard orbit propagation.
//!
//! ## Technical scope
//!
//! [`Propagator`] is generic over the integrator `I: AdaptiveStepper`.  Use
//! the provided type aliases [`Dop853Propagator`], [`Dopri5Propagator`], or
//! [`Rk4Propagator`] (wraps [`FixedRk4Adapter`]) to select the algorithm at
//! compile time, consistent with the codebase's type-level model-selection
//! convention.
//!
//! `PropagatorConfig` carries tolerances, output cadence, and step bounds
//! independent of the integrator choice.
//!
//! ## References
//!
//! * Hairer, Nørsett, Wanner, *Solving Ordinary Differential Equations I*.

use std::marker::PhantomData;

use principia::{
    self, propagate, rk4_step, AccelerationModel, AdaptiveStepper, Dop853, Dopri5, EventDetector,
    PropagationConfig, PropagationResult, StateTransitionMatrix,
};

use crate::astro::dynamics::context::DynamicsContext;
use crate::astro::dynamics::errors::DynamicsError;
use crate::astro::dynamics::state::OrbitState;
use crate::coordinates::centers::{Geocentric, ReferenceCenter};
use crate::coordinates::frames::{ReferenceFrame, GCRS};
use crate::qtty::Second;
use crate::time::TT;

// ─────────────────────────────────────────────────────────────────────────────
// FixedRk4Adapter
// ─────────────────────────────────────────────────────────────────────────────

/// Adapts the fixed-step [`principia::Rk4`] stepper to the [`AdaptiveStepper`]
/// interface expected by [`Propagator`].
///
/// Construct with the desired sub-step magnitude; the adapter caps each
/// requested `h_try` at `step` so the effective step never exceeds `step`.
///
/// Use [`Rk4Propagator`] as a convenience alias for
/// `Propagator<FixedRk4Adapter, …>`.
#[derive(Debug, Clone, Copy)]
pub struct FixedRk4Adapter {
    /// Fixed sub-step size used for every RK4 stage.
    pub step: Second,
}

impl FixedRk4Adapter {
    /// Create an adapter with the given sub-step magnitude.
    pub fn new(step: Second) -> Self {
        Self { step }
    }
}

impl<C, F> AdaptiveStepper<DynamicsContext, TT, C, F> for FixedRk4Adapter
where
    C: ReferenceCenter,
    F: ReferenceFrame,
{
    fn step<M: AccelerationModel<DynamicsContext, TT, C, F>>(
        &self,
        model: &M,
        state: &OrbitState<C, F>,
        h_try: Second,
        ctx: &DynamicsContext,
    ) -> Result<(OrbitState<C, F>, Second, Second, u32), principia::PrincipiaError> {
        if !h_try.value().is_finite() {
            return Err(principia::PrincipiaError::InvalidStepRequest {
                reason: "h_try must be finite",
            });
        }
        let max_h = self.step.value().abs();
        if !max_h.is_finite() || max_h == 0.0 {
            return Err(principia::PrincipiaError::InvalidStepRequest {
                reason: "FixedRk4Adapter.step must be finite and non-zero",
            });
        }
        let h = if h_try.value().abs() < max_h {
            h_try
        } else {
            Second::new(h_try.value().signum() * max_h)
        };
        let new_state = rk4_step(model, state, h, ctx)?;
        Ok((new_state, h, h, 0))
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PropagatorConfig
// ─────────────────────────────────────────────────────────────────────────────

/// Configuration for [`Propagator`], independent of the integrator choice.
///
/// The integrator algorithm is encoded in the type parameter `I` of
/// [`Propagator`], not here.  This config carries only tunable knobs that are
/// common to all integrators.
#[derive(Debug, Clone)]
pub struct PropagatorConfig {
    /// Adaptive-step tolerances (used by DOPRI5 and DOP853; ignored by RK4).
    pub tolerances: principia::IntegratorTolerances,
    /// Maximum number of accepted steps before the propagation is aborted.
    pub max_steps: usize,
    /// Optional regular output sampling interval for event-recording runs.
    pub event_dt: Option<Second>,
}

impl Default for PropagatorConfig {
    fn default() -> Self {
        Self {
            tolerances: principia::IntegratorTolerances::uniform(1e-9, 1e-6, 1e-9),
            max_steps: 1_000_000,
            event_dt: None,
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Propagator
// ─────────────────────────────────────────────────────────────────────────────

/// Orbit propagator facade, typed over the integrator algorithm.
///
/// `I` must implement [`AdaptiveStepper`] for the
/// `(DynamicsContext, TT, C, F)` tuple. Use the provided type aliases to
/// select the algorithm at compile time:
///
/// ```rust
/// use siderust::astro::dynamics::propagation::{Dop853Propagator, PropagatorConfig};
/// use siderust::astro::dynamics::context::DynamicsContext;
/// use siderust::astro::dynamics::forces::TwoBody;
/// use siderust::astro::dynamics::GM_EARTH;
///
/// let _p: Dop853Propagator = Dop853Propagator::new(
///     principia::Dop853::new(principia::IntegratorTolerances::uniform(1e-9, 1e-6, 1e-9)),
///     TwoBody::new(GM_EARTH),
///     PropagatorConfig::default(),
///     DynamicsContext::empty(),
/// );
/// ```
pub struct Propagator<I, C = Geocentric, F = GCRS, FM = principia::TwoBody>
where
    I: AdaptiveStepper<DynamicsContext, TT, C, F>,
    C: ReferenceCenter,
    F: ReferenceFrame,
    FM: AccelerationModel<DynamicsContext, TT, C, F>,
{
    integrator: I,
    model: FM,
    config: PropagatorConfig,
    ctx: DynamicsContext,
    _marker: PhantomData<(C, F)>,
}

impl<I, C, F, FM> Propagator<I, C, F, FM>
where
    I: AdaptiveStepper<DynamicsContext, TT, C, F>,
    C: ReferenceCenter,
    F: ReferenceFrame,
    FM: AccelerationModel<DynamicsContext, TT, C, F>,
{
    /// Construct a propagator with the given integrator, force model,
    /// configuration, and runtime context.
    pub fn new(integrator: I, model: FM, config: PropagatorConfig, ctx: DynamicsContext) -> Self {
        Self {
            integrator,
            model,
            config,
            ctx,
            _marker: PhantomData,
        }
    }
}

impl<I, FM> Propagator<I, Geocentric, GCRS, FM>
where
    I: AdaptiveStepper<DynamicsContext, TT, Geocentric, GCRS>,
    FM: AccelerationModel<DynamicsContext, TT, Geocentric, GCRS>,
{
    /// Propagate `state` forward by `dt` seconds.
    pub fn propagate(
        &self,
        state: OrbitState<Geocentric, GCRS>,
        dt: Second,
    ) -> Result<OrbitState<Geocentric, GCRS>, DynamicsError> {
        let t_start = state.epoch;
        let t_end = state.epoch + dt;
        let h0 = Second::new(30.0);
        let cfg = PropagationConfig::<DynamicsContext, TT, Geocentric, GCRS>::new(t_start, t_end)
            .with_initial_step(h0)
            .with_max_steps(self.config.max_steps)
            .with_tolerances(self.config.tolerances);
        let result = propagate(&self.integrator, &self.model, state, &cfg, &self.ctx)
            .map_err(|e: principia::PropagationError| DynamicsError::from(e))?;
        Ok(result
            .samples
            .into_iter()
            .last()
            .expect("propagation always produces at least one sample"))
    }

    /// Propagate `state` forward by `dt` seconds, recording crossings of `event`.
    pub fn propagate_with_events<E>(
        &self,
        state: OrbitState<Geocentric, GCRS>,
        dt: Second,
        event: E,
    ) -> Result<PropagationResult<TT, Geocentric, GCRS>, DynamicsError>
    where
        E: EventDetector<DynamicsContext, TT, Geocentric, GCRS> + 'static,
    {
        let t_start = state.epoch;
        let t_end = state.epoch + dt;
        let mut cfg =
            PropagationConfig::<DynamicsContext, TT, Geocentric, GCRS>::new(t_start, t_end)
                .with_initial_step(Second::new(30.0))
                .with_max_steps(self.config.max_steps)
                .with_tolerances(self.config.tolerances)
                .with_event(Box::new(event));
        if let Some(dt_ev) = self.config.event_dt {
            cfg = cfg.with_output_every(dt_ev);
        }
        propagate(&self.integrator, &self.model, state, &cfg, &self.ctx).map_err(Into::into)
    }

    /// Propagate `state` forward by `dt` seconds and return the final state
    /// together with its State Transition Matrix.
    pub fn propagate_with_stm(
        &self,
        state: OrbitState<Geocentric, GCRS>,
        dt: Second,
    ) -> Result<(OrbitState<Geocentric, GCRS>, StateTransitionMatrix<GCRS>), DynamicsError> {
        principia::propagate_stm(&self.model, state, dt, &self.ctx).map_err(Into::into)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Convenience type aliases
// ─────────────────────────────────────────────────────────────────────────────

/// [`Propagator`] using the DOP853 adaptive integrator.
pub type Dop853Propagator<C = Geocentric, F = GCRS, FM = principia::TwoBody> =
    Propagator<Dop853, C, F, FM>;

/// [`Propagator`] using the Dormand–Prince DOPRI5 adaptive integrator.
pub type Dopri5Propagator<C = Geocentric, F = GCRS, FM = principia::TwoBody> =
    Propagator<Dopri5, C, F, FM>;

/// [`Propagator`] using the fixed-step classical RK4 integrator via
/// [`FixedRk4Adapter`].
pub type Rk4Propagator<C = Geocentric, F = GCRS, FM = principia::TwoBody> =
    Propagator<FixedRk4Adapter, C, F, FM>;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::astro::dynamics::context::DynamicsContext;
    use crate::astro::dynamics::forces::{TwoBody, GM_EARTH};
    use crate::astro::dynamics::{Position, Velocity};
    use crate::qtty::{Kilometers, Second};
    use crate::time::JulianDate;

    const R: f64 = 7_000.0;

    fn circular_state() -> (OrbitState, f64) {
        let v = (GM_EARTH.value() / R).sqrt();
        let s0 = OrbitState::new(
            JulianDate::new(2_451_545.0).to_j2000s(),
            Position::<GCRS>::new(R, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, v, 0.0),
        );
        let period = 2.0 * std::f64::consts::PI * (R.powi(3) / GM_EARTH.value()).sqrt();
        (s0, period)
    }

    fn eccentric_state() -> OrbitState {
        let rp = 7_000.0;
        let ra = 8_000.0;
        let a = 0.5 * (rp + ra);
        let vp = (GM_EARTH.value() * (2.0 / rp - 1.0 / a)).sqrt();
        OrbitState::new(
            JulianDate::new(2_451_545.0).to_j2000s(),
            Position::<GCRS>::new(rp, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, vp, 0.0),
        )
    }

    #[test]
    fn propagator_dop853_orbit_closes() {
        let (s0, period) = circular_state();
        let p = Dop853Propagator::new(
            Dop853::new(principia::IntegratorTolerances::uniform(1e-9, 1e-6, 1e-9)),
            TwoBody::new(GM_EARTH),
            PropagatorConfig::default(),
            DynamicsContext::empty(),
        );
        let s = p.propagate(s0, Second::new(period)).unwrap();
        let dr = ((s.position.x().value() - R).powi(2) + s.position.y().value().powi(2)).sqrt();
        assert!(dr < 1e-3);
    }

    #[test]
    fn propagator_with_events_records_crossing() {
        let cfg = PropagatorConfig {
            event_dt: Some(Second::new(10.0)),
            ..PropagatorConfig::default()
        };
        let p = Dop853Propagator::new(
            Dop853::new(cfg.tolerances),
            TwoBody::new(GM_EARTH),
            cfg,
            DynamicsContext::empty(),
        );
        let event = principia::RadialThresholdEvent::new(Kilometers::new(7_500.0));
        let result = p
            .propagate_with_events(eccentric_state(), Second::new(4_000.0), event)
            .unwrap();
        assert!(!result.events.is_empty());
    }
}
