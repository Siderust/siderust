// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Generic adaptive propagation driver.
//!
//! Composes any [`AdaptiveStepper`] with a [`PropagatorConfig`] (event
//! detectors, output cadence, step budget) into a single high-level call.

use super::config::PropagatorConfig;
use super::result::{EventOccurrence, PropagationResult};
use crate::astro::dynamics::context::DynamicsContext;
use crate::astro::dynamics::errors::DynamicsError;
use crate::astro::dynamics::forces::ForceModel;
use crate::astro::dynamics::integrators::AdaptiveStepper;
use crate::astro::dynamics::state::OrbitState;
use crate::coordinates::centers::{Geocentric, ReferenceCenter};
use crate::coordinates::frames::{ReferenceFrame, GCRS};
use crate::qtty::Second;

/// Drive an [`AdaptiveStepper`] from `cfg.t_start` to `cfg.t_end`.
///
/// Honors `output_every`, `output_at`, terminal events, and a hard
/// `max_steps` budget.
#[allow(clippy::too_many_lines)]
pub fn propagate<I, C, F, FM>(
    integrator: &I,
    force: &FM,
    initial: OrbitState<C, F>,
    cfg: &PropagatorConfig<C, F>,
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

        let (new_state, h_used, h_next) =
            integrator.step(force, &state, Second::new(h), ctx)?;
        steps_taken += 1;
        if steps_taken > cfg.max_steps {
            return Err(DynamicsError::InvalidStepRequest {
                reason: "propagator max_steps exceeded",
            });
        }

        // Event detection: sign change between previous and new state.
        for (i, ev) in cfg.events.iter().enumerate() {
            let g_new = ev.evaluate(&new_state, ctx)?;
            if g_prev[i] == 0.0 || (g_prev[i] * g_new) < 0.0 {
                events.push(EventOccurrence {
                    event_name: ev.name(),
                    state: new_state.clone(),
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
                    next_output_t =
                        Some(t_target + Second::new(direction * dt.value().abs()));
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
        let _ = h_used;
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
        steps_rejected: 0,
    })
}

// Re-export marker types so doc links resolve.
#[allow(unused_imports)]
use crate as _siderust;
#[allow(unused_imports)]
use Geocentric as _Geocentric;
#[allow(unused_imports)]
use GCRS as _Gcrs;
