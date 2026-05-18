// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Generic adaptive propagation driver.
//!
//! ## Scope
//!
//! Provides [`propagate`], the main entry point for orbit propagation.  It
//! composes any [`super::super::integrators::AdaptiveStepper`] with a
//! [`super::config::PropagationConfig`] to manage event detection, output
//! sampling, step control, and termination logic.
//!
//! ## Workflow
//!
//! 1. Supply an adaptive integrator (DOPRI5, DOP853), force model, initial
//!    state, configuration, and dynamics context.
//! 2. The driver steps the integrator while:
//!    - Checking event detector switching functions for zero crossings.
//!    - Recording output at requested times (`output_every`, `output_at`).
//!    - Honoring the step budget (`max_steps`).
//!    - Terminating on terminal events or reaching `t_end`.
//! 3. Return a [`PropagationResult`] with samples, event times, and statistics.
//!
//! ## Event detection
//!
//! Events are zero-crossing detectors: if `g(t) changes sign between steps,
//! a zero-crossing is recorded.  Terminal events immediately halt propagation.
//!
//! ## Failure modes
//!
//! Returns errors on:
//! - Step control failure (integrator error)
//! - Step shrinking below `h_min`
//! - Step budget exhaustion
//! - Event evaluation error
//!
//! ## References
//!
//! * Vallado, *Fundamentals of Astrodynamics and Applications* (2013), §4.4.

use super::config::PropagationConfig;
use super::result::{EventOccurrence, PropagationResult};
use crate::astro::dynamics::context::DynamicsContext;
use crate::astro::dynamics::errors::DynamicsError;
use crate::astro::dynamics::forces::ForceModel;
use crate::astro::dynamics::integrators::AdaptiveStepper;
use crate::astro::dynamics::state::OrbitState;
use crate::coordinates::centers::{Geocentric, ReferenceCenter};
use crate::coordinates::frames::{ReferenceFrame, GCRS};
use crate::qtty::Second;

/// Linear interpolation between two [`OrbitState`] values at fraction `theta`
/// in `[0, 1]`.  Used to estimate the state at an event or output crossing
/// between two integration steps.
fn lerp_state<C: ReferenceCenter, F: ReferenceFrame>(
    a: &OrbitState<C, F>,
    b: &OrbitState<C, F>,
    theta: f64,
) -> OrbitState<C, F>
where
    C::Params: Clone,
{
    use crate::astro::dynamics::state::VelocityUnit;
    use crate::coordinates::cartesian::Velocity as CartVelocity;
    use crate::qtty::unit::Kilometer;
    let px = a.position.x().value() + theta * (b.position.x().value() - a.position.x().value());
    let py = a.position.y().value() + theta * (b.position.y().value() - a.position.y().value());
    let pz = a.position.z().value() + theta * (b.position.z().value() - a.position.z().value());
    let vx = a.velocity.x().value() + theta * (b.velocity.x().value() - a.velocity.x().value());
    let vy = a.velocity.y().value() + theta * (b.velocity.y().value() - a.velocity.y().value());
    let vz = a.velocity.z().value() + theta * (b.velocity.z().value() - a.velocity.z().value());
    let epoch = a.epoch + Second::new(theta * (b.epoch - a.epoch).value());
    OrbitState {
        epoch,
        position: crate::coordinates::cartesian::Position::<C, F, Kilometer>::new_with_params(
            a.position.center_params().clone(),
            px,
            py,
            pz,
        ),
        velocity: CartVelocity::<F, VelocityUnit>::new(vx, vy, vz),
    }
}

/// Drive an [`AdaptiveStepper`] from `cfg.t_start` to `cfg.t_end`.
///
/// Honors `output_every`, `output_at`, terminal events, and a hard
/// `max_steps` budget.
#[allow(clippy::too_many_lines)]
pub fn propagate<I, C, F, FM>(
    integrator: &I,
    force: &FM,
    initial: OrbitState<C, F>,
    cfg: &PropagationConfig<C, F>,
    ctx: &DynamicsContext,
) -> Result<PropagationResult<C, F>, DynamicsError>
where
    I: AdaptiveStepper<C, F>,
    FM: ForceModel<C, F>,
    C: ReferenceCenter,
    F: ReferenceFrame,
    C::Params: Clone,
{
    let total = cfg.total_duration_s();
    let direction = if total >= 0.0 { 1.0 } else { -1.0 };

    let mut samples: Vec<OrbitState<C, F>> = Vec::new();
    let mut events: Vec<EventOccurrence<C, F>> = Vec::new();
    let mut steps_taken: u32 = 0;
    let mut steps_rejected: u32 = 0;

    samples.push(initial.clone());

    if total == 0.0 {
        return Ok(PropagationResult {
            samples,
            events,
            steps_taken,
            steps_rejected: 0,
        });
    }

    // Pre-evaluate the switching function for each event.
    let mut g_prev: Vec<f64> = Vec::with_capacity(cfg.events.len());
    for ev in &cfg.events {
        g_prev.push(ev.evaluate(&initial, ctx)?);
    }

    let mut state = initial;
    let mut h = direction * cfg.h0.value().abs().min(cfg.h_max.value().abs());

    let mut next_output_t = cfg
        .output_every
        .map(|dt| cfg.t_start + Second::new(direction * dt.value().abs()));
    let mut output_at_iter = cfg.output_at.iter().peekable();
    let mut terminated = false;

    while !terminated {
        let remaining = (cfg.t_end - state.epoch).value();
        if remaining * direction <= 1e-12 {
            break;
        }
        if h.abs() > remaining.abs() {
            h = remaining;
        }
        if h.abs() > cfg.h_max.value().abs() {
            h = direction * cfg.h_max.value().abs();
        }
        if h.abs() < cfg.h_min.value().abs() && h.abs() < remaining.abs() {
            // Underflow guard: still advance with the minimum allowed step.
            h = direction * cfg.h_min.value().abs();
        }

        // Clip h to not overshoot upcoming output times so the stored state
        // is at the exact requested epoch (no interpolation needed).
        if let Some(t_target) = next_output_t {
            let to_target = (t_target - state.epoch).value();
            if to_target.abs() > 1e-12 && to_target * direction >= 0.0 && to_target.abs() < h.abs()
            {
                h = to_target;
            }
        }
        if let Some(&t_target) = output_at_iter.peek() {
            let to_target = (*t_target - state.epoch).value();
            if to_target.abs() > 1e-12 && to_target * direction >= 0.0 && to_target.abs() < h.abs()
            {
                h = to_target;
            }
        }

        let (new_state, _h_used, h_next, rejected) =
            integrator.step(force, &state, Second::new(h), ctx)?;
        steps_taken += 1;
        steps_rejected += rejected;
        if steps_taken > cfg.max_steps {
            return Err(DynamicsError::InvalidStepRequest {
                reason: "propagator max_steps exceeded",
            });
        }

        // Event detection: sign change between previous and new state.
        for (i, ev) in cfg.events.iter().enumerate() {
            let g_new = ev.evaluate(&new_state, ctx)?;
            if g_prev[i] == 0.0 || (g_prev[i] * g_new) < 0.0 {
                // Linear interpolation to estimate the state at the crossing.
                let theta = if (g_new - g_prev[i]).abs() > 1e-300 {
                    g_prev[i].abs() / (g_prev[i].abs() + g_new.abs())
                } else {
                    0.5
                };
                let crossing_state = lerp_state(&state, &new_state, theta);
                events.push(EventOccurrence {
                    event_name: ev.name(),
                    state: crossing_state,
                });
                if ev.terminal() {
                    terminated = true;
                }
            }
            g_prev[i] = g_new;
        }

        // Output cadence: emit a sample if a requested output time was
        // crossed by this step.
        if let Some(t_target) = next_output_t {
            let crossed = (t_target - state.epoch).value() * direction >= 0.0
                && (new_state.epoch - t_target).value() * direction >= 0.0;
            if crossed {
                samples.push(new_state.clone());
                if let Some(dt) = cfg.output_every {
                    next_output_t = Some(t_target + Second::new(direction * dt.value().abs()));
                }
            }
        }
        while let Some(&t_target) = output_at_iter.peek() {
            let crossed = (*t_target - state.epoch).value() * direction >= 0.0
                && (new_state.epoch - *t_target).value() * direction >= 0.0;
            if crossed {
                samples.push(new_state.clone());
                output_at_iter.next();
            } else {
                break;
            }
        }

        state = new_state;
        h = h_next.value();
    }

    // Always include the final state.
    if samples
        .last()
        .map(|s| (s.epoch - state.epoch).value().abs() > 0.0)
        .unwrap_or(true)
    {
        samples.push(state);
    }

    Ok(PropagationResult {
        samples,
        events,
        steps_taken,
        steps_rejected,
    })
}

// Re-export marker types so doc links resolve.
#[allow(unused_imports)]
use crate as _siderust;
#[allow(unused_imports)]
use Geocentric as _Geocentric;
#[allow(unused_imports)]
use GCRS as _Gcrs;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::astro::dynamics::context::DynamicsContext;
    use crate::astro::dynamics::forces::two_body::TwoBody;
    use crate::astro::dynamics::integrators::Dopri5;
    use crate::astro::dynamics::propagation::events::AltitudeEvent;
    use crate::astro::dynamics::state::OrbitState;
    use crate::astro::dynamics::{Position, Velocity};
    use crate::coordinates::frames::GCRS;
    use crate::ext_qtty::tolerances::IntegratorTolerances;
    use crate::qtty::{Kilometers, Second};
    use crate::time::JulianDate;

    const MU: f64 = 398_600.441_8;
    const R: f64 = 7_000.0;

    fn circ_state() -> OrbitState {
        let v = (MU / R).sqrt();
        OrbitState::new_at_jd(
            JulianDate::new(2_451_545.0),
            Position::<GCRS>::new(R, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, v, 0.0),
        )
    }

    fn tol() -> IntegratorTolerances {
        IntegratorTolerances::uniform(1e-9, 1e-6, 1e-9)
    }

    #[test]
    fn zero_duration_returns_initial_state_only() {
        let s0 = circ_state();
        let cfg = PropagationConfig::new(s0.epoch, s0.epoch);
        let stepper = Dopri5::new(tol());
        let result = propagate(
            &stepper,
            &TwoBody::earth(),
            s0,
            &cfg,
            &DynamicsContext::empty(),
        )
        .unwrap();
        assert_eq!(result.samples.len(), 1);
        assert_eq!(result.steps_taken, 0);
    }

    #[test]
    fn propagate_forward_one_orbit_arrives() {
        let s0 = circ_state();
        let period = 2.0 * std::f64::consts::PI * (R.powi(3) / MU).sqrt();
        let t_end = s0.epoch + Second::new(period);
        let cfg = PropagationConfig::new(s0.epoch, t_end);
        let stepper = Dopri5::new(tol());
        let result = propagate(
            &stepper,
            &TwoBody::earth(),
            s0,
            &cfg,
            &DynamicsContext::empty(),
        )
        .unwrap();
        let final_state = result.samples.last().unwrap();
        let dr = ((final_state.position.x().value() - R).powi(2)
            + final_state.position.y().value().powi(2)
            + final_state.position.z().value().powi(2))
        .sqrt();
        assert!(dr < 1.0, "orbit closure error {dr} km exceeds 1 km");
    }

    #[test]
    fn output_every_produces_samples() {
        let s0 = circ_state();
        let period = 2.0 * std::f64::consts::PI * (R.powi(3) / MU).sqrt();
        let t_end = s0.epoch + Second::new(period);
        let cfg =
            PropagationConfig::new(s0.epoch, t_end).with_output_every(Second::new(period / 10.0));
        let stepper = Dopri5::new(tol());
        let result = propagate(
            &stepper,
            &TwoBody::earth(),
            s0,
            &cfg,
            &DynamicsContext::empty(),
        )
        .unwrap();
        assert!(
            result.samples.len() > 3,
            "expected dense output samples, got {}",
            result.samples.len()
        );
    }

    #[test]
    fn output_at_appends_samples() {
        let s0 = circ_state();
        let period = 2.0 * std::f64::consts::PI * (R.powi(3) / MU).sqrt();
        let t_target = s0.epoch + Second::new(period * 0.25);
        let t_end = s0.epoch + Second::new(period);
        let cfg = PropagationConfig::new(s0.epoch, t_end).with_output_at(vec![t_target]);
        let stepper = Dopri5::new(tol());
        let result = propagate(
            &stepper,
            &TwoBody::earth(),
            s0,
            &cfg,
            &DynamicsContext::empty(),
        )
        .unwrap();
        assert!(result.samples.len() >= 3, "expected at least 3 samples");
    }

    #[test]
    fn event_detection_fires() {
        let v_circ = (MU / R).sqrt();
        let s0_ecc = OrbitState::new_at_jd(
            JulianDate::new(2_451_545.0),
            Position::<GCRS>::new(R, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, v_circ * 1.01, 0.0),
        );
        let period = 2.0 * std::f64::consts::PI * (R.powi(3) / MU).sqrt();
        let t_end = s0_ecc.epoch + Second::new(period);
        let event = AltitudeEvent::new(Kilometers::new(750.0), Kilometers::new(6_371.0));
        let cfg = PropagationConfig::new(s0_ecc.epoch, t_end).with_event(Box::new(event));
        let stepper = Dopri5::new(tol());
        let result = propagate(
            &stepper,
            &TwoBody::earth(),
            s0_ecc,
            &cfg,
            &DynamicsContext::empty(),
        )
        .unwrap();
        assert!(
            !result.events.is_empty(),
            "altitude event must fire for eccentric orbit"
        );
    }
}
