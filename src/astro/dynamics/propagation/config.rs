// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Low-level driver configuration for the [`super::driver::propagate`] function.
//!
//! ## Scope
//!
//! Provides [`PropagationConfig<C, F>`], the typed configuration bundle for the
//! adaptive propagation driver: time bounds, step hints, tolerances, event
//! detectors, output cadence, and step budgeting.
//!
//! See [`super::propagator::Propagator`] and the high-level API for the
//! user-facing wrapper.
//!
//! ## Typical usage
//!
//! ```rust,ignore
//! use siderust::astro::dynamics::propagation::PropagationConfig;
//! use siderust::time::{Time, TT, JulianDate};
//! use siderust::qtty::Second;
//!
//! let t_start = Time::<TT>::from_jd(JulianDate::new(2_451_545.0));
//! let t_end = t_start + Second::new(86_400.0);
//!
//! let cfg = PropagationConfig::new(t_start, t_end)
//!     .with_initial_step(Second::new(30.0))
//!     .with_max_step(Second::new(600.0))
//!     .with_output_every(Second::new(60.0));
//! ```

use super::events::EventDetector;
use crate::coordinates::centers::{Geocentric, ReferenceCenter};
use crate::coordinates::frames::{ReferenceFrame, GCRS};
use crate::ext_qtty::tolerances::IntegratorTolerances;
use crate::qtty::Second;
use crate::time::{Time, TT};

/// Low-level driver configuration.  Use [`super::propagator::PropagatorConfig`]
/// + [`super::propagator::Propagator`] for the ergonomic high-level API.
pub struct PropagationConfig<C = Geocentric, F = GCRS>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
{
    pub t_start: Time<TT>,
    pub t_end: Time<TT>,
    pub h0: Second,
    pub h_max: Second,
    pub h_min: Second,
    pub tolerances: IntegratorTolerances,
    pub output_every: Option<Second>,
    pub output_at: Vec<Time<TT>>,
    pub events: Vec<Box<dyn EventDetector<C, F>>>,
    pub max_steps: u32,
}

impl<C, F> PropagationConfig<C, F>
where
    C: ReferenceCenter,
    F: ReferenceFrame,
{
    pub fn new(t_start: Time<TT>, t_end: Time<TT>) -> Self {
        Self {
            t_start,
            t_end,
            h0: Second::new(30.0),
            h_max: Second::new(86_400.0),
            h_min: Second::new(1e-6),
            tolerances: IntegratorTolerances::uniform(1e-9, 1e-3, 1e-6),
            output_every: None,
            output_at: Vec::new(),
            events: Vec::new(),
            max_steps: 1_000_000,
        }
    }

    pub fn with_initial_step(mut self, h0: Second) -> Self {
        self.h0 = h0;
        self
    }
    pub fn with_max_step(mut self, h_max: Second) -> Self {
        self.h_max = h_max;
        self
    }
    pub fn with_min_step(mut self, h_min: Second) -> Self {
        self.h_min = h_min;
        self
    }
    pub fn with_tolerances(mut self, tols: IntegratorTolerances) -> Self {
        self.tolerances = tols;
        self
    }
    pub fn with_output_every(mut self, dt: Second) -> Self {
        self.output_every = Some(dt);
        self
    }
    pub fn with_output_at(mut self, times: Vec<Time<TT>>) -> Self {
        self.output_at = times;
        self
    }
    pub fn with_event(mut self, ev: Box<dyn EventDetector<C, F>>) -> Self {
        self.events.push(ev);
        self
    }
    pub fn with_max_steps(mut self, max_steps: u32) -> Self {
        self.max_steps = max_steps;
        self
    }

    /// Returns the total signed duration `t_end - t_start` in seconds.
    pub fn total_duration_s(&self) -> f64 {
        (self.t_end - self.t_start).value()
    }
}
