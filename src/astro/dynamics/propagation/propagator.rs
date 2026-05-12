// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! High-level [`Propagator`] facade that dispatches over the three integrators.
//!
//! ## Quick start
//!
//! ```rust
//! use siderust::astro::dynamics::propagation::{Propagator, PropagatorConfig, IntegratorChoice};
//! use siderust::astro::dynamics::forces::TwoBody;
//! use siderust::astro::dynamics::context::DynamicsContext;
//! use siderust::astro::dynamics::{OrbitState, Position, Velocity};
//! use siderust::coordinates::frames::GCRS;
//! use siderust::time::JulianDate;
//! use siderust::qtty::{Second, IntegratorTolerances};
//!
//! let s0 = OrbitState::new_at_jd(
//!     JulianDate::new(2_451_545.0),
//!     Position::<GCRS>::new(7_000.0, 0.0, 0.0),
//!     Velocity::<GCRS>::new(0.0, 7.5450, 0.0),
//! );
//! let cfg = PropagatorConfig {
//!     integrator: IntegratorChoice::Dop853,
//!     tolerances: IntegratorTolerances::uniform(1e-9, 1e-6, 1e-9),
//!     max_steps: 1_000_000,
//!     event_dt: None,
//! };
//! let p = Propagator::new(TwoBody::earth(), cfg, DynamicsContext::empty());
//! let s1 = p.propagate(s0, Second::new(600.0)).unwrap();
//! assert!(s1.epoch != s0.epoch);
//! ```

use std::marker::PhantomData;

use crate::astro::dynamics::context::DynamicsContext;
use crate::astro::dynamics::errors::DynamicsError;
use crate::astro::dynamics::forces::ForceModel;
use crate::astro::dynamics::integrators::rk4::rk4_step;
use crate::astro::dynamics::integrators::{
    dop853_propagate, dopri5_propagate, rk4_propagate, AdaptiveStepper, Dop853, Dopri5,
};
use crate::astro::dynamics::state::OrbitState;
use crate::astro::dynamics::variational;
use crate::astro::dynamics::variational::StateTransitionMatrix;
use crate::coordinates::centers::{Geocentric, ReferenceCenter};
use crate::coordinates::frames::{ReferenceFrame, GCRS};
use crate::ext_qtty::tolerances::IntegratorTolerances;
use crate::qtty::Second;

use super::config::PropagationConfig;
use super::events::EventDetector;
use super::result::PropagationResult;

// =============================================================================
// Public types
// =============================================================================

/// Integrator selection for [`Propagator`].
#[derive(Debug, Clone, Copy)]
pub enum IntegratorChoice {
    /// Classical fixed-step RK4. `step` is the sub-step size in seconds.
    Rk4 { step: Second },
    /// Dormand–Prince adaptive RK4(5).
    Dopri5,
    /// Hairer DOP853 adaptive 8th-order Runge–Kutta.
    Dop853,
}

/// High-level configuration for [`Propagator`].
///
/// This type is intentionally minimal: it captures the integrator choice,
/// tolerance, step budget, and optional dense-output cadence.  Time windows
/// and initial conditions are supplied at call time via
/// [`Propagator::propagate`] / [`Propagator::propagate_with_events`].
#[derive(Debug, Clone)]
pub struct PropagatorConfig {
    /// Which integrator to use.
    pub integrator: IntegratorChoice,
    /// Tolerances forwarded to adaptive integrators (`Dopri5`, `Dop853`).
    /// Ignored for fixed-step `Rk4`.
    pub tolerances: IntegratorTolerances,
    /// Maximum number of steps before the propagation is aborted.
    pub max_steps: usize,
    /// Dense-output sampling cadence forwarded to the event driver.
    pub event_dt: Option<Second>,
}

impl Default for PropagatorConfig {
    fn default() -> Self {
        Self {
            integrator: IntegratorChoice::Dop853,
            tolerances: IntegratorTolerances::uniform(1e-9, 1e-6, 1e-9),
            max_steps: 1_000_000,
            event_dt: None,
        }
    }
}

// =============================================================================
// Propagator
// =============================================================================

/// Orbit propagator facade.
///
/// Wraps a [`ForceModel`], a [`PropagatorConfig`], and a [`DynamicsContext`].
/// Dispatches to RK4, DOPRI5, or DOP853 depending on the configured
/// [`IntegratorChoice`].
///
/// # Type parameters
///
/// - `C` — reference center (defaults to [`Geocentric`]).
/// - `F` — reference frame (defaults to [`GCRS`]).
/// - `FM` — force model, implementing [`ForceModel<C, F>`].
///
/// # Current dispatch support
///
/// The propagation methods ([`propagate`][Self::propagate],
/// [`propagate_with_stm`][Self::propagate_with_stm],
/// [`propagate_with_events`][Self::propagate_with_events]) are currently
/// implemented only for the default `Geocentric`/`GCRS` frame pair because
/// the underlying concrete integrators ([`Dop853`], [`Dopri5`]) are
/// specialized to that frame.  The struct itself is generic so that
/// higher-order code can store typed propagators without loss of information.
pub struct Propagator<C = Geocentric, F = GCRS, FM = crate::astro::dynamics::forces::TwoBody>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
    FM: ForceModel<C, F>,
{
    force: FM,
    config: PropagatorConfig,
    ctx: DynamicsContext,
    _marker: PhantomData<(C, F)>,
}

impl<C, F, FM> Propagator<C, F, FM>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
    FM: ForceModel<C, F>,
{
    /// Create a new propagator.
    pub fn new(force: FM, config: PropagatorConfig, ctx: DynamicsContext) -> Self {
        Self {
            force,
            config,
            ctx,
            _marker: PhantomData,
        }
    }
}

// =============================================================================
// Geocentric/GCRS-specific dispatch methods
//
// The underlying integrators (Dop853, Dopri5) implement AdaptiveStepper only
// for the Geocentric/GCRS pair.  FixedRk4Adapter is generic but for
// consistency all dispatch methods are placed in the same impl block.
// =============================================================================

impl<FM> Propagator<Geocentric, GCRS, FM>
where
    FM: ForceModel<Geocentric, GCRS>,
{
    /// Propagate `state` forward (or backward) by `dt` seconds.
    ///
    /// Dispatches to the integrator specified in [`PropagatorConfig::integrator`].
    pub fn propagate(
        &self,
        state: OrbitState<Geocentric, GCRS>,
        dt: Second,
    ) -> Result<OrbitState<Geocentric, GCRS>, DynamicsError> {
        match self.config.integrator {
            IntegratorChoice::Rk4 { step } => {
                let total_s = dt.value();
                let step_s = step.value().abs();
                let n = ((total_s.abs() / step_s).ceil() as usize).max(1);
                let h = Second::new(total_s / n as f64);
                rk4_propagate(&self.force, state, h, n, &self.ctx)
            }
            IntegratorChoice::Dopri5 => {
                dopri5_propagate(&self.force, state, dt, self.config.tolerances, &self.ctx)
            }
            IntegratorChoice::Dop853 => {
                dop853_propagate(&self.force, state, dt, self.config.tolerances, &self.ctx)
            }
        }
    }

    /// Propagate with an event detector.
    ///
    /// Returns a [`PropagationResult`] that includes samples and any event
    /// occurrences detected during propagation.  The event is passed by
    /// value and boxed internally so it can be held by the driver config.
    pub fn propagate_with_events<E>(
        &self,
        state: OrbitState<Geocentric, GCRS>,
        dt: Second,
        event: E,
    ) -> Result<PropagationResult<Geocentric, GCRS>, DynamicsError>
    where
        E: EventDetector<Geocentric, GCRS> + 'static,
    {
        let t_start = state.epoch;
        let t_end = state.epoch + dt;

        let h0 = match self.config.integrator {
            IntegratorChoice::Rk4 { step } => step,
            _ => Second::new(30.0),
        };
        let mut cfg = PropagationConfig::<Geocentric, GCRS>::new(t_start, t_end)
            .with_initial_step(h0)
            .with_max_steps(self.config.max_steps.min(u32::MAX as usize) as u32)
            .with_event(Box::new(event));
        if let Some(dt_ev) = self.config.event_dt {
            cfg = cfg.with_output_every(dt_ev);
        }

        match self.config.integrator {
            IntegratorChoice::Rk4 { step } => {
                let adapter = FixedRk4Adapter { step };
                super::driver::propagate(&adapter, &self.force, state, &cfg, &self.ctx)
            }
            IntegratorChoice::Dopri5 => {
                let stepper = Dopri5::new(self.config.tolerances);
                super::driver::propagate(&stepper, &self.force, state, &cfg, &self.ctx)
            }
            IntegratorChoice::Dop853 => {
                let stepper = Dop853::new(self.config.tolerances);
                super::driver::propagate(&stepper, &self.force, state, &cfg, &self.ctx)
            }
        }
    }

    /// Propagate the orbit state and the 6×6 state-transition matrix Φ(t, t₀)
    /// jointly using the variational-equations RK4 integrator.
    ///
    /// Returns `(final_state, Φ)` where Φ is a [`StateTransitionMatrix`] (`FrameMatrix6`)
    /// tagged to [`GCRS`].
    pub fn propagate_with_stm(
        &self,
        state: OrbitState<Geocentric, GCRS>,
        dt: Second,
    ) -> Result<(OrbitState<Geocentric, GCRS>, StateTransitionMatrix<GCRS>), DynamicsError> {
        variational::propagate_stm(&self.force, state, dt, &self.ctx)
    }
}

// =============================================================================
// Internal: thin RK4 adapter so the fixed-step integrator works with the
// generic adaptive driver (used by propagate_with_events + Rk4).
// =============================================================================

struct FixedRk4Adapter {
    step: Second,
}

impl AdaptiveStepper<Geocentric, GCRS> for FixedRk4Adapter {
    fn step<FM: ForceModel<Geocentric, GCRS>>(
        &self,
        force: &FM,
        state: &OrbitState<Geocentric, GCRS>,
        h_try: Second,
        ctx: &DynamicsContext,
    ) -> Result<(OrbitState<Geocentric, GCRS>, Second, Second, u32), DynamicsError> {
        // Respect the driver's step clip: use h_try if it is smaller in magnitude.
        let h = if h_try.value().abs() < self.step.value().abs() {
            h_try
        } else {
            self.step
        };
        let new_state = rk4_step(force, state, h, ctx)?;
        Ok((new_state, h, h, 0))
    }
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::astro::dynamics::context::DynamicsContext;
    use crate::astro::dynamics::forces::TwoBody;
    use crate::astro::dynamics::propagation::events::AltitudeEvent;
    use crate::astro::dynamics::{Position, Velocity};
    use crate::coordinates::frames::GCRS;
    use crate::ext_qtty::tolerances::IntegratorTolerances;
    use crate::qtty::{Kilometers, Second};
    use crate::time::JulianDate;

    const MU: f64 = 398_600.441_8;
    const R: f64 = 7_000.0;

    fn circular_state() -> (OrbitState, f64) {
        let v = (MU / R).sqrt();
        let s0 = OrbitState::new_at_jd(
            JulianDate::new(2_451_545.0),
            Position::<GCRS>::new(R, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, v, 0.0),
        );
        let period = 2.0 * std::f64::consts::PI * (R.powi(3) / MU).sqrt();
        (s0, period)
    }

    fn dop853_config() -> PropagatorConfig {
        PropagatorConfig {
            integrator: IntegratorChoice::Dop853,
            tolerances: IntegratorTolerances::uniform(1e-9, 1e-6, 1e-9),
            max_steps: 1_000_000,
            event_dt: None,
        }
    }

    fn dopri5_config() -> PropagatorConfig {
        PropagatorConfig {
            integrator: IntegratorChoice::Dopri5,
            tolerances: IntegratorTolerances::uniform(1e-9, 1e-6, 1e-9),
            max_steps: 1_000_000,
            event_dt: None,
        }
    }

    fn rk4_config(step_s: f64) -> PropagatorConfig {
        PropagatorConfig {
            integrator: IntegratorChoice::Rk4 {
                step: Second::new(step_s),
            },
            tolerances: IntegratorTolerances::uniform(1e-9, 1e-6, 1e-9),
            max_steps: 1_000_000,
            event_dt: None,
        }
    }

    fn orbit_closure_error(s_end: &OrbitState) -> f64 {
        let dx = s_end.position.x().value() - R;
        let dy = s_end.position.y().value();
        let dz = s_end.position.z().value();
        (dx * dx + dy * dy + dz * dz).sqrt()
    }

    // --- Orbit closure tests -----------------------------------------------

    #[test]
    fn propagator_dop853_orbit_closes() {
        let (s0, period) = circular_state();
        let p = Propagator::new(TwoBody::earth(), dop853_config(), DynamicsContext::empty());
        let s = p.propagate(s0, Second::new(period)).unwrap();
        let dr = orbit_closure_error(&s);
        assert!(dr < 1e-3, "DOP853 orbit closure {dr} km exceeds 1e-3 km");
    }

    #[test]
    fn propagator_dopri5_orbit_closes() {
        let (s0, period) = circular_state();
        let p = Propagator::new(TwoBody::earth(), dopri5_config(), DynamicsContext::empty());
        let s = p.propagate(s0, Second::new(period)).unwrap();
        let dr = orbit_closure_error(&s);
        assert!(dr < 1.0, "DOPRI5 orbit closure {dr} km exceeds 1 km");
    }

    #[test]
    fn propagator_rk4_orbit_closes() {
        let (s0, period) = circular_state();
        // Use a small enough step for reasonable accuracy.
        let p = Propagator::new(TwoBody::earth(), rk4_config(10.0), DynamicsContext::empty());
        let s = p.propagate(s0, Second::new(period)).unwrap();
        let dr = orbit_closure_error(&s);
        assert!(dr < 5.0, "RK4 orbit closure {dr} km exceeds 5 km");
    }

    // --- STM identity test --------------------------------------------------

    #[test]
    fn propagate_with_stm_identity_at_zero_dt() {
        let (s0, _) = circular_state();
        let p = Propagator::new(TwoBody::earth(), dop853_config(), DynamicsContext::empty());
        let (_, phi) = p.propagate_with_stm(s0, Second::new(0.0)).unwrap();
        let arr = phi.as_array();
        for i in 0..6 {
            for j in 0..6 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (arr[i][j] - expected).abs() < 1e-12,
                    "phi[{i}][{j}] = {} (expected {expected})",
                    arr[i][j]
                );
            }
        }
    }

    // --- Event detection test -----------------------------------------------

    #[test]
    fn propagate_with_events_altitude_event_fires() {
        let (_, period) = circular_state();
        let mu = MU;
        let r = R;
        let v_circ = (mu / r).sqrt();
        // Add 1% perturbation to velocity to create eccentricity so altitude
        // oscillates and crosses the trigger threshold.
        let s0_ecc = OrbitState::new_at_jd(
            JulianDate::new(2_451_545.0),
            Position::<GCRS>::new(r, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, v_circ * 1.01, 0.0),
        );
        let earth_radius = Kilometers::new(6_371.0);
        // Periapsis altitude ≈ 629 km; apoapsis altitude ≈ 913 km.
        // Trigger at 750 km (r = 7121 km) sits between them so the spacecraft
        // crosses the threshold twice per orbit.
        let trigger_alt = Kilometers::new(750.0);
        let event = AltitudeEvent::new(trigger_alt, earth_radius);

        let p = Propagator::new(TwoBody::earth(), dop853_config(), DynamicsContext::empty());
        let result = p
            .propagate_with_events(s0_ecc, Second::new(period), event)
            .unwrap();
        assert!(
            !result.events.is_empty(),
            "expected at least one altitude event, got 0"
        );
    }

    #[test]
    fn propagator_config_default_uses_dop853() {
        let cfg = PropagatorConfig::default();
        assert!(matches!(cfg.integrator, IntegratorChoice::Dop853));
        assert_eq!(cfg.max_steps, 1_000_000);
        assert!(cfg.event_dt.is_none());
    }

    #[test]
    fn propagator_propagate_with_events_dopri5_fires() {
        let v_circ = (MU / R).sqrt();
        let s0_ecc = OrbitState::new_at_jd(
            JulianDate::new(2_451_545.0),
            Position::<GCRS>::new(R, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, v_circ * 1.01, 0.0),
        );
        let period = 2.0 * std::f64::consts::PI * (R.powi(3) / MU).sqrt();
        let event = AltitudeEvent::new(Kilometers::new(750.0), Kilometers::new(6_371.0));
        let p = Propagator::new(TwoBody::earth(), dopri5_config(), DynamicsContext::empty());
        let result = p
            .propagate_with_events(s0_ecc, Second::new(period), event)
            .unwrap();
        assert!(!result.events.is_empty());
    }

    #[test]
    fn propagator_event_dt_produces_dense_output() {
        let (s0, period) = circular_state();
        let mut cfg = dop853_config();
        cfg.event_dt = Some(Second::new(period / 10.0));
        let event = AltitudeEvent::new(Kilometers::new(100.0), Kilometers::new(6_371.0));
        let p = Propagator::new(TwoBody::earth(), cfg, DynamicsContext::empty());
        let result = p
            .propagate_with_events(s0, Second::new(period), event)
            .unwrap();
        assert!(
            result.samples.len() > 3,
            "expected dense samples, got {}",
            result.samples.len()
        );
    }

    #[test]
    fn propagator_propagate_with_events_rk4_runs() {
        let (s0, period) = circular_state();
        let event = AltitudeEvent::new(Kilometers::new(100.0), Kilometers::new(6_371.0));
        let p = Propagator::new(TwoBody::earth(), rk4_config(60.0), DynamicsContext::empty());
        let result = p
            .propagate_with_events(s0, Second::new(period), event)
            .unwrap();
        assert!(!result.samples.is_empty());
    }
}
